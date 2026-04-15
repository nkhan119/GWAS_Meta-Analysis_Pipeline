// ================================================================
// Module: META_ANALYSIS
// Fixed-effects IVW  (weight = 1/SE²)
// Random-effects IVW (DerSimonian–Laird)
// Cochran Q  + I²  per SNP
// Genomic-control corrected Z-scores
// ================================================================

process META_ANALYSIS {

    tag "${trait}"
    label 'medium'
    publishDir "${params.out_dir}/meta/${trait}", mode: 'copy', overwrite: true

    input:
    tuple val(trait), val(type), path(harmonised)

    output:
    tuple val(trait), val(type), path("${trait}_meta_results.tsv.gz"), emit: results
    tuple val(trait), val(type), path("${trait}_het_stats.tsv.gz"),    emit: het_stats
    tuple val(trait), val(type), path("${trait}_qq_data.tsv.gz"),      emit: qq_data
    path "${trait}_meta_summary.tsv",                                   emit: summary

    script:
    """
    run_meta_analysis.R \\
        --input       ${harmonised} \\
        --trait       ${trait} \\
        --type        ${type} \\
        --method      ${params.meta_method} \\
        --p_sig       ${params.p_sig} \\
        --p_sug       ${params.p_suggestive} \\
        --i2_thresh   ${params.i2_threshold} \\
        --min_cohorts ${params.min_cohorts} \\
        --threads     ${task.cpus} \\
        --out_results ${trait}_meta_results.tsv.gz \\
        --out_het     ${trait}_het_stats.tsv.gz \\
        --out_qq      ${trait}_qq_data.tsv.gz \\
        --out_summary ${trait}_meta_summary.tsv
    """
}
