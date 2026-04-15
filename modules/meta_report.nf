// ================================================================
// Module: META_REPORT
// Generates the interactive HTML report with embedded Plotly charts
// ================================================================

process META_REPORT {

    label 'small'
    publishDir "${params.out_dir}", mode: 'copy', overwrite: true

    input:
    path meta_files,    stageAs: 'meta/*'
    path het_files,     stageAs: 'het/*'
    path qq_files,      stageAs: 'qq/*'
    path top_loci_files,stageAs: 'loci/*'
    path filter_files,  stageAs: 'filter/*'
    val  author
    val  affiliation
    val  institute

    output:
    path "Meta_Report.html"

    script:
    """
    python3 /opt/pipeline/bin/generate_meta_report.py \\
        --meta_dir    meta/ \\
        --het_dir     het/ \\
        --qq_dir      qq/ \\
        --loci_dir    loci/ \\
        --filter_dir  filter/ \\
        --author      "${author}" \\
        --affiliation "${affiliation}" \\
        --institute   "${institute}" \\
        --p_sig       ${params.p_sig} \\
        --p_sug       ${params.p_suggestive} \\
        --i2_thresh   ${params.i2_threshold} \\
        --out          Meta_Report.html
    """
}
