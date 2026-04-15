// ================================================================
// Module: HETEROGENEITY_FILTER
// 1. Flag loci with I² > threshold
// 2. Positional LD clumping for independent signals
// 3. Produce top-loci table for annotation
// Input : (trait, type, path(meta_results))
// Output: top_loci → (trait, type, path)
//         summary  → path
//         high_het → path
// ================================================================

process HETEROGENEITY_FILTER {

    tag "${trait}"
    label 'medium'
    publishDir "${params.out_dir}/filtered/${trait}", mode: 'copy', overwrite: true

    input:
    tuple val(trait), val(type), path(meta_results)

    output:
    tuple val(trait), val(type), path("${trait}_top_loci.tsv"), emit: top_loci
    path "${trait}_filter_summary.tsv",                          emit: summary
    path "${trait}_high_het_loci.tsv",                          emit: high_het

    script:
    """
    heterogeneity_filter.py \\
        --input       ${meta_results} \\
        --trait       ${trait} \\
        --type        ${type} \\
        --p_sig       ${params.p_sig} \\
        --p_sug       ${params.p_suggestive} \\
        --i2_thresh   ${params.i2_threshold} \\
        --clump_p1    ${params.clump_p1} \\
        --clump_p2    ${params.clump_p2} \\
        --clump_r2    ${params.clump_r2} \\
        --clump_kb    ${params.clump_kb} \\
        --out_loci    ${trait}_top_loci.tsv \\
        --out_summary ${trait}_filter_summary.tsv \\
        --out_het     ${trait}_high_het_loci.tsv
    """
}
