#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// ================================================================
// GWAS Meta-Analysis Pipeline — CDC 1.0.0 · Phase B
// Author  : Nadeem Khan, Bioinformatician, INRS-CAFSB
// GitHub  : github.com/nkhan119
// ================================================================

include { VALIDATE_SUMSTATS          } from './modules/validate_sumstats.nf'
include { HARMONISE_META             } from './modules/harmonise_meta.nf'
include { META_ANALYSIS              } from './modules/meta_analysis.nf'
include { HETEROGENEITY_FILTER       } from './modules/heterogeneity_filter.nf'
include { ANNOTATION                 } from './modules/annotation.nf'
include { META_REPORT                } from './modules/meta_report.nf'
include { META_FIGURES               } from './modules/meta_figures.nf'

// ── Banner ─────────────────────────────────────────────────────
log.info """
╔══════════════════════════════════════════════════════════════════╗
║   GWAS Meta-Analysis Pipeline · CDC 1.0.0 · Phase B             ║
║   Author    : ${params.author.padRight(48)}║
║   Institute : ${(params.institute ?: "").take(48).padRight(48)}║
╠══════════════════════════════════════════════════════════════════╣
║   sumstats_dir  : ${params.sumstats_dir.take(46).padRight(46)}║
║   out_dir       : ${params.out_dir.take(46).padRight(46)}║
║   traits (QT)   : ${params.quant_traits.take(46).padRight(46)}║
║   traits (BT)   : ${params.binary_traits.take(46).padRight(46)}║
╠══════════════════════════════════════════════════════════════════╣
║   method        : ${params.meta_method.padRight(46)}║
║   skip_validate : ${params.skip_validate.toString().padRight(46)}║
║   skip_annotate : ${params.skip_annotate.toString().padRight(46)}║
╚══════════════════════════════════════════════════════════════════╝
""".stripIndent()

// ── Helper: parse trait list ────────────────────────────────────
def parseTraits(raw) {
    raw instanceof List ? raw : raw.toString().split(",").collect { it.trim() }
}

workflow {

    def qt_traits = parseTraits(params.quant_traits)
    def bt_traits = parseTraits(params.binary_traits)
    def all_traits = qt_traits + bt_traits

    // ── B0  Discover per-cohort sumstat files ──────────────────
    // Expects: {sumstats_dir}/{Cohort_*}/sumstats/{Cohort_X}_{trait}_sumstats.tsv.gz
    // The trait key must be cohort-agnostic so files from Cohort_A and Cohort_B
    // group together.  Strip the leading "{Cohort_X}_" prefix from the filename.
    ch_raw_ss = Channel
        .fromPath("${params.sumstats_dir}/**/sumstats/*_sumstats.tsv.gz")
        .map { f ->
            def parts  = f.toString().tokenize('/')
            def cohort = parts.find { it.startsWith('Cohort_') } ?: 'unknown'
            // Remove compression + _sumstats suffix
            def base   = f.baseName.replaceAll('_sumstats\\.tsv', '').replaceAll('\\.gz','')
            // Strip the leading cohort prefix (e.g. "Cohort_A_") so the key is
            // the bare trait name shared across cohorts (e.g. "bmi", "cad")
            def trait  = base.replaceFirst(/^${cohort}_/, '')
            def type   = trait in qt_traits ? 'qt' : (trait in bt_traits ? 'bt' : 'qt')
            tuple(trait, cohort, type, f)
        }

    // Group by trait → list of (cohort, file) tuples per trait
    ch_by_trait = ch_raw_ss
        .groupTuple(by: 0)
        .map { trait, cohorts, types, files ->
            def cohort_files = [cohorts, files].transpose()
            def type = types[0]
            tuple(trait, type, cohort_files.collect { c, f -> tuple(c, f) })
        }

    // ── B1  Validate & QC sumstats ────────────────────────────
    // ch_by_trait shape: (trait, type, [[cohort,file],...])
    // Reshape to (trait, type, [names], [paths]) so VALIDATE_SUMSTATS
    // can stage files via path().  The validated output is re-grouped
    // and names are carried through for HARMONISE_META.
    if (!params.skip_validate) {
        ch_for_validate = ch_by_trait.map { trait, type, cohort_files ->
            def names = cohort_files.collect { pair -> pair[0] }
            def paths = cohort_files.collect { pair -> pair[1] }
            tuple(trait, type, names, paths)
        }
        VALIDATE_SUMSTATS(ch_for_validate)

        // VALIDATE_SUMSTATS emits one file per cohort per trait.
        // Re-group by [trait, type]; infer cohort name from filename.
        // Output shape: (trait, type, [cohort_names], [staged_paths])
        ch_validated = VALIDATE_SUMSTATS.out.validated_files
            .groupTuple(by: [0, 1])
            .map { trait, type, files ->
                def flist = files instanceof List ? files.flatten() : [files]
                // Derive cohort name from validated filename:
                // e.g. cad_validated_Cohort_A.tsv.gz → Cohort_A
                def names = flist.collect { f ->
                    (f.name =~ /validated_(.+)\.tsv\.gz/)[0][1]
                }
                tuple(trait, type, names, flist)
            }
    } else {
        // Skip validation — pass raw files directly, names from ch_by_trait
        ch_validated = ch_by_trait.map { trait, type, cohort_files ->
            def names = cohort_files.collect { pair -> pair[0] }
            def paths = cohort_files.collect { pair -> pair[1] }
            tuple(trait, type, names, paths)
        }
    }

    // ── B2  Harmonise to common variant set ───────────────────
    // HARMONISE_META expects (trait, type, val(names), path(files))
    // so Nextflow stages files properly into the work directory.
    ch_harmonised = HARMONISE_META(ch_validated)

    // ── B3  Meta-analysis (IVW + RE + FE) ────────────────────
    ch_meta = META_ANALYSIS(ch_harmonised.harmonised)

    // ── B4  Heterogeneity filtering + conditional signals ─────
    ch_filtered = HETEROGENEITY_FILTER(ch_meta.results)

    // ── B5  Functional annotation ─────────────────────────────
    if (!params.skip_annotate) {
        ch_annotated = ANNOTATION(ch_filtered.top_loci)
        ch_top_loci  = ch_annotated.annotated.map { trait, f -> f }.collect()
    } else {
        ch_top_loci  = ch_filtered.top_loci.map { trait, type, f -> f }.collect()
    }

    // ── B6  Collect all outputs ────────────────────────────────
    ch_all_meta   = ch_meta.results.map   { t, typ, f -> f }.collect()
    ch_all_filter = ch_filtered.summary.collect()
    ch_het_stats  = ch_meta.het_stats.map { t, typ, f -> f }.collect()
    ch_qq_data    = ch_meta.qq_data.map   { t, typ, f -> f }.collect()

    // ── B7  Publication-quality figures ───────────────────────
    META_FIGURES(
        ch_all_meta,
        ch_het_stats,
        ch_qq_data,
        ch_top_loci,
        params.author,
        params.institute ?: "INRS-CAFSB"
    )

    // ── B8  Interactive HTML report ───────────────────────────
    META_REPORT(
        ch_all_meta,
        ch_het_stats,
        ch_qq_data,
        ch_top_loci,
        ch_all_filter,
        params.author,
        params.affiliation,
        params.institute ?: "INRS-CAFSB"
    )
}

// ── Completion banner ──────────────────────────────────────────
workflow.onComplete {
    def status = workflow.success ? "SUCCESS ✓" : "FAILED ✗"
    log.info """
╔══════════════════════════════════════════════════════════════════╗
║   Meta-Analysis ${status.padRight(47)}║
║   Duration  : ${workflow.duration.toString().padRight(50)}║
║   HTML      : ${(params.out_dir + '/Meta_Report.html').padRight(50)}║
║   Figures   : ${(params.out_dir + '/figures/').padRight(50)}║
╚══════════════════════════════════════════════════════════════════╝
""".stripIndent()
}
