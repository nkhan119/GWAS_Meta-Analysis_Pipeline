// ================================================================
// Module: META_FIGURES
// Publication-quality figures (300 dpi PNG + PDF):
//   - Miami plots (bidirectional Manhattan)
//   - Forest plots per lead SNP
//   - QQ plots with genomic lambda annotation
//   - Heterogeneity heatmap
//   - Effect-size comparison scatter (Cohort A vs B)
//   - Per-trait LD clump bubble plot
// ================================================================

process META_FIGURES {

    label 'medium'
    publishDir "${params.out_dir}/figures", mode: 'copy', overwrite: true, saveAs: { fn -> fn }

    input:
    path meta_files,    stageAs: 'meta/*'
    path het_files,     stageAs: 'het/*'
    path qq_files,      stageAs: 'qq/*'
    path top_loci_files,stageAs: 'loci/*'
    val  author
    val  institute

    output:
    // Figures are written to typed sub-directories by generate_meta_figures.py
    path "manhattan/*.png", optional: true
    path "qqplot/*.png",    optional: true
    path "forest/*.png",    optional: true
    path "scatter/*.png",   optional: true
    path "global/*.png",    optional: true
    path "**/*.pdf",        optional: true

    script:
    """
    python3 /opt/pipeline/bin/generate_meta_figures.py \\
        --meta_dir   meta/ \\
        --het_dir    het/ \\
        --qq_dir     qq/ \\
        --loci_dir   loci/ \\
        --author     "${author}" \\
        --institute  "${institute}" \\
        --p_sig      ${params.p_sig} \\
        --p_sug      ${params.p_suggestive} \\
        --i2_thresh  ${params.i2_threshold} \\
        --out_dir    ./
    """
}
