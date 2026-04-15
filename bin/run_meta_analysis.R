#!/usr/bin/env Rscript
# ================================================================
# run_meta_analysis.R

# Author: Nadeem Khan, INRS-CAFSB
# ================================================================

suppressPackageStartupMessages({
    library(data.table)
    library(metafor)
    library(optparse)
    library(glue)
    library(R.utils)
})

# ── Arguments ──────────────────────────────────────────────────
opt_list <- list(
    make_option("--input",       type="character"),
    make_option("--trait",       type="character"),
    make_option("--type",        type="character", default="qt"),
    make_option("--method",      type="character", default="ivw_re"),
    make_option("--p_sig",       type="double",    default=5e-8),
    make_option("--p_sug",       type="double",    default=1e-5),
    make_option("--i2_thresh",   type="double",    default=75),
    make_option("--min_cohorts", type="integer",   default=2),
    make_option("--threads",     type="integer",   default=4),
    make_option("--out_results", type="character"),
    make_option("--out_het",     type="character"),
    make_option("--out_qq",      type="character"),
    make_option("--out_summary", type="character")
)
opt <- parse_args(OptionParser(option_list=opt_list))

cat(glue("\n[meta] Trait: {opt$trait}  |  Method: {opt$method}  |  Type: {opt$type}\n\n"))

# ── Load harmonised data ───────────────────────────────────────
cat("[meta] Loading harmonised data...\n")
dt <- fread(opt$input, showProgress=FALSE)

# Expected columns: SNP, CHR, POS, EA, OA, EAF,
#                   BETA_{cohort}, SE_{cohort}, N_{cohort}
beta_cols <- grep("^BETA_", names(dt), value=TRUE)
se_cols   <- grep("^SE_",   names(dt), value=TRUE)
n_cols    <- grep("^N_",    names(dt), value=TRUE)
cohorts   <- sub("^BETA_", "", beta_cols)

cat(glue("[meta] Cohorts detected: {paste(cohorts, collapse=', ')}\n"))
cat(glue("[meta] Variants loaded: {nrow(dt)}\n\n"))

if (nrow(dt) == 0) {
    cat("[meta] ERROR: harmonised file contains 0 variants.\n")
    cat("[meta] This usually means the cohort SNP identifiers do not overlap.\n")
    cat("[meta] Check that both cohorts use the same SNP ID format (e.g. rsID vs CHR:POS:REF:ALT).\n")
    quit(status=1)
}

# ── Core meta-analysis function ───────────────────────────────
run_meta_snp <- function(betas, ses, method="DL") {
    # Filter out cohorts with missing data
    ok  <- !is.na(betas) & !is.na(ses) & ses > 0
    nb  <- sum(ok)
    if (nb < 2) return(list(beta=NA, se=NA, p=NA, q=NA, i2=NA, n_cohorts=nb))

    b <- betas[ok]
    s <- ses[ok]
    w <- 1 / s^2

    # Fixed-effects IVW
    beta_fe <- sum(w * b) / sum(w)
    se_fe   <- sqrt(1 / sum(w))
    z_fe    <- beta_fe / se_fe
    p_fe    <- 2 * pnorm(-abs(z_fe))

    # Cochran Q
    q_stat  <- sum(w * (b - beta_fe)^2)
    df      <- nb - 1
    p_q     <- pchisq(q_stat, df=df, lower.tail=FALSE)
    i2      <- max(0, (q_stat - df) / q_stat * 100)

    # DerSimonian-Laird tau²
    tau2    <- max(0, (q_stat - df) / (sum(w) - sum(w^2)/sum(w)))
    w_re    <- 1 / (s^2 + tau2)
    beta_re <- sum(w_re * b) / sum(w_re)
    se_re   <- sqrt(1 / sum(w_re))
    z_re    <- beta_re / se_re
    p_re    <- 2 * pnorm(-abs(z_re))

    list(
        beta_fe   = beta_fe, se_fe  = se_fe,  p_fe  = p_fe,
        beta_re   = beta_re, se_re  = se_re,  p_re  = p_re,
        q_stat    = q_stat,  p_q    = p_q,    i2    = i2,
        tau2      = tau2,    n_cohorts = nb
    )
}

# ── Apply per SNP ─────────────────────────────────────────────
cat("[meta] Running meta-analysis across", nrow(dt), "variants...\n")

beta_mat <- as.matrix(dt[, ..beta_cols])
se_mat   <- as.matrix(dt[, ..se_cols])

results <- lapply(seq_len(nrow(dt)), function(i) {
    run_meta_snp(beta_mat[i,], se_mat[i,])
})

res_dt <- rbindlist(results)
dt_out <- cbind(dt[, .(SNP, CHR, POS, EA, OA, EAF)], res_dt)

# Genomic control lambda
p_use <- if ("p_re" %in% names(dt_out)) dt_out$p_re else dt_out$p_fe
p_use <- p_use[!is.na(p_use) & p_use > 0 & p_use < 1]
lambda_gc <- round(median(qchisq(1 - p_use, 1)) / qchisq(0.5, 1), 4)
cat(glue("[meta] Genomic lambda (GC): {lambda_gc}\n"))

# Apply GC correction to RE p-value
if (isTRUE(lambda_gc > 1.05)) {
    cat("[meta] Applying GC correction...\n")
    z_corrected     <- qnorm(1 - dt_out$p_re / 2) / sqrt(lambda_gc)
    dt_out$p_re_gc  <- 2 * pnorm(-abs(z_corrected))
} else {
    dt_out$p_re_gc  <- dt_out$p_re   # lambda <= 1.05 or NA: no correction needed
}

# Significance flags
dt_out[, sig    := p_re_gc < opt$p_sig]
dt_out[, sug    := p_re_gc < opt$p_sug & !sig]
dt_out[, hi_het := i2 > opt$i2_thresh]

n_sig <- sum(dt_out$sig, na.rm=TRUE)
n_sug <- sum(dt_out$sug, na.rm=TRUE)
cat(glue("[meta] Genome-wide significant SNPs: {n_sig}\n"))
cat(glue("[meta] Suggestive SNPs: {n_sug}\n\n"))

# ── Write results ─────────────────────────────────────────────
cat("[meta] Writing results...\n")
fwrite(dt_out, opt$out_results, sep="\t", compress="gzip")

# Heterogeneity stats table
het_dt <- dt_out[!is.na(q_stat), .(SNP, CHR, POS, EA, OA,
                                    q_stat, p_q, i2, tau2, n_cohorts)]
fwrite(het_dt, opt$out_het, sep="\t", compress="gzip")

# QQ data (subsample for speed)
set.seed(42)
p_col     <- "p_re_gc"
qq_p      <- dt_out[[p_col]][!is.na(dt_out[[p_col]]) & dt_out[[p_col]] > 0]
n_qq      <- length(qq_p)
# Keep all top signals; subsample remainder
top_idx   <- which(qq_p < opt$p_sug)
other_idx <- setdiff(seq_len(n_qq), top_idx)
samp_idx  <- c(top_idx, sample(other_idx, min(200000, length(other_idx))))
qq_dt     <- data.table(
    observed = sort(-log10(qq_p[samp_idx]), decreasing=TRUE),
    expected = -log10(ppoints(length(samp_idx)))
)
fwrite(qq_dt, opt$out_qq, sep="\t", compress="gzip")

# Summary table
summ <- data.table(
    trait         = opt$trait,
    type          = opt$type,
    n_variants    = nrow(dt_out),
    n_cohorts     = length(cohorts),
    cohorts       = paste(cohorts, collapse=";"),
    n_sig         = n_sig,
    n_sug         = n_sug,
    lambda_gc     = lambda_gc,
    method        = opt$method,
    top_snp       = if (n_sig > 0) dt_out[sig==TRUE][which.min(p_re_gc)]$SNP else NA_character_,
    top_p         = if (n_sig > 0) min(dt_out$p_re_gc, na.rm=TRUE) else NA_real_,
    top_beta      = if (n_sig > 0) dt_out[sig==TRUE][which.min(p_re_gc)]$beta_re else NA_real_
)
fwrite(summ, opt$out_summary, sep="\t")

cat(glue("[meta] Done: {opt$trait}\n\n"))
