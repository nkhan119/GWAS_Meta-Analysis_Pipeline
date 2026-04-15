// ================================================================
// Module: HARMONISE_META
// Aligns alleles across cohorts, filters to common variants,
// applies GC correction if lambda > 1.1.
//
// Input channel shape (from main.nf B2):
//   tuple val(trait), val(type), val(cohort_names), path(input_files)
//
// cohort_names : list of cohort name strings  [Cohort_A, Cohort_B]
// input_files  : staged list of Path objects  (work-dir local copies)
//
// Separating names (val) from files (path) is required so Nextflow
// can stage every file into the work directory before execution.
//
// Output:
//   harmonised → tuple(trait, type, path)
//   lambda     → path
// ================================================================

process HARMONISE_META {

    tag "${trait}"
    label 'small'
    publishDir "${params.out_dir}/harmonised/${trait}", mode: 'copy', overwrite: true

    input:
    tuple val(trait), val(type), val(cohort_names), path(input_files)

    output:
    tuple val(trait), val(type), path("${trait}_harmonised.tsv.gz"), emit: harmonised
    path "${trait}_lambda_report.tsv",                                emit: lambda

    script:
    def names = cohort_names instanceof List ? cohort_names : [cohort_names]
    def files = input_files instanceof List  ? input_files  : [input_files]
    def file_args = [names, files].transpose()
        .collect { name, f -> "--input ${name}:${f}" }
        .join(" \\\n        ")
    """
    harmonise_meta.py \\
        --trait       ${trait} \\
        --type        ${type} \\
        ${file_args} \\
        --min_cohorts ${params.min_cohorts} \\
        --ref_build   ${params.ref_build} \\
        --out         ${trait}_harmonised.tsv.gz \\
        --lambda_out  ${trait}_lambda_report.tsv
    """
}
