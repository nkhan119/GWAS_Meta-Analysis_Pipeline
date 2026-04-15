# GWAS Meta-Analysis Pipeline — CDC 1.0.0 · Phase B

**Author:** Nadeem Khan, Bioinformatician, INRS-Centre Armand-Frappier Santé-Biotechnologie  
**GitHub:** [github.com/nkhan119](https://github.com/nkhan119)

---

## What this pipeline does

Picks up where `gwas-pipeline-v3` (Phase A) leaves off. It reads the harmonised
summary statistics produced by REGENIE → `HARMONISE_SUMSTATS` and runs a full
cross-cohort meta-analysis:

```
sumstats (Cohort_A + Cohort_B)
  │
  ├── B1  VALIDATE_SUMSTATS    — column checks, SE/P filters, palindromic flags
  ├── B2  HARMONISE_META       — allele alignment, common-variant merge, GC lambda
  ├── B3  META_ANALYSIS        — IVW fixed-effects + DerSimonian-Laird RE, Cochran Q, I²
  ├── B4  HETEROGENEITY_FILTER — I² flagging, positional LD clumping, independent loci
  ├── B5  ANNOTATION           — nearest gene, consequence class (GRCh37/38)
  ├── B7  META_FIGURES         — Manhattan, QQ, Forest, Heatmap, Effect scatter (PNG+PDF)
  └── B8  META_REPORT          — interactive HTML with embedded Plotly
```

---

## Directory structure expected

```
{sumstats_dir}/
  Cohort_A/
    sumstats/
      ldl_cholesterol_sumstats.tsv.gz
      bmi_sumstats.tsv.gz
      cad_sumstats.tsv.gz
      ...
  Cohort_B/
    sumstats/
      ldl_cholesterol_sumstats.tsv.gz
      ...
```

This matches the output of `gwas-pipeline-v3` exactly — no reformatting needed.

---

## Quick start

### 1. Build Docker image (once)

```bash
cd gwas-meta-pipeline/
docker build -t nkhan119/gwas-meta:1.0.0 .
```

### 2. Run locally (Docker)

```bash
nextflow run main.nf \
    -profile local,docker \
    --sumstats_dir ~/Downloads/1000G_phase3_common_norel/GWAS_analysis \
    --out_dir      ~/Downloads/1000G_phase3_common_norel/META_results \
    -resume
```

### 3. Run on Narval (Singularity + SLURM)

```bash
nextflow run main.nf \
    -profile slurm \
    --base_dir   ~/projects/1000G \
    --sumstats_dir ~/projects/1000G/GWAS_analysis \
    --out_dir    ~/projects/1000G/META_results \
    -resume
```

### 4. Skip validation (if sumstats are already clean)

```bash
nextflow run main.nf \
    -profile local,docker \
    --sumstats_dir ... \
    --skip_validate true \
    --skip_annotate false \
    -resume
```

---

## Outputs

```
results/
  Meta_Report.html
  meta/{trait}/
    {trait}_meta_results.tsv.gz
    {trait}_het_stats.tsv.gz
    {trait}_qq_data.tsv.gz
    {trait}_meta_summary.tsv
  filtered/{trait}/
    {trait}_top_loci.tsv
    {trait}_filter_summary.tsv
    {trait}_high_het_loci.tsv
  annotated/{trait}/
    {trait}_annotated_loci.tsv
  figures/
    {trait}_manhattan.png
    {trait}_qqplot.png
    {trait}_forest.png
    {trait}_effect_scatter.png
    heterogeneity_heatmap.png
    lambda_summary.png
  validation/{trait}/
    {trait}_validated_{cohort}.tsv.gz
    {trait}_validation_report.tsv
  harmonised/{trait}/
    {trait}_harmonised.tsv.gz
    {trait}_lambda_report.tsv
  logs/
    trace_*.txt
    nf_report_*.html
    timeline_*.html
    dag.svg
```

---

## Key parameters

| Parameter | Default | Description |
|---|---|---|
| `--sumstats_dir` | — | Path containing Cohort_A/ Cohort_B/ |
| `--out_dir` | `./results` | Output directory |
| `--quant_traits` | `ldl_cholesterol,bmi,crp_log,height_cm` | Quantitative traits |
| `--binary_traits` | `cad,t2d,hypertension` | Binary traits |
| `--meta_method` | `ivw_re` | `ivw_re` \| `ivw_fe` \| `both` |
| `--p_sig` | `5e-8` | GW significance threshold |
| `--p_suggestive` | `1e-5` | Suggestive threshold |
| `--i2_threshold` | `75` | I² % to flag high heterogeneity |
| `--min_cohorts` | `2` | Minimum cohorts per variant |
| `--clump_kb` | `500` | LD clumping window (kb) |
| `--ref_build` | `GRCh37` | Reference build |
| `--skip_validate` | `false` | Skip B1 validation |
| `--skip_annotate` | `false` | Skip B5 annotation |

---

## Conda fallback (no Docker)

```bash
conda env create -f envs/meta_env.yaml
conda activate cdc_meta

nextflow run main.nf \
    -profile local,conda \
    --sumstats_dir ... \
    -resume
```

---

## Dependencies

All bundled in the Docker image (`nkhan119/gwas-meta:1.0.0`):

- **R:** `metafor`, `meta`, `data.table`, `TwoSampleMR`, `coloc`
- **Python:** `pandas`, `numpy`, `scipy`, `plotly`, `statsmodels`
- **System:** `plink2`, Ubuntu 22.04

---

## Citation

> Khan N. (2026). CDC 1.0.0 — Causal Deep Consensus GWAS Meta-Analysis Pipeline.
> INRS-Centre Armand-Frappier Santé-Biotechnologie.
> https://github.com/nkhan119

---

*Phase A (QC → PCA → REGENIE) lives in [gwas-pipeline-v3](https://github.com/nkhan119/gwas-pipeline-v3).*
