// ================================================================
// Module: VALIDATE_SUMSTATS
// Validates and QC-filters per-cohort summary statistics.
//
// Input channel shape (reshaped in main.nf):
//   tuple val(trait), val(type), val(cohort_names), path(cohort_files)
//
// Using path() for cohort_files ensures Nextflow stages the files
// into the process work directory before execution — this is what
// caused the FileNotFoundError when files were passed via val().
//
// Output:
//   validated_files → tuple(trait, type, path("*_validated_*.tsv.gz"))
//   report          → path
// ================================================================

process VALIDATE_SUMSTATS {

    tag "${trait}"
    label 'small'
    publishDir "${params.out_dir}/validation/${trait}", mode: 'copy', overwrite: true

    input:
    tuple val(trait), val(type), val(cohort_names), path(cohort_files)

    output:
    tuple val(trait), val(type), path("${trait}_validated_*.tsv.gz"), emit: validated_files
    path "${trait}_validation_report.tsv",                             emit: report

    script:
    // cohort_names: list of cohort name strings
    // cohort_files: staged list of Path objects (work-dir local copies)
    def names = cohort_names instanceof List ? cohort_names : [cohort_names]
    def files = cohort_files instanceof List ? cohort_files : [cohort_files]
    def cf_args = [names, files].transpose()
        .collect { name, f -> "--cohort ${name}:${f}" }
        .join(' \\\n        ')
    """
    validate_sumstats.py \\
        --trait      ${trait} \\
        --type       ${type} \\
        ${cf_args} \\
        --min_n      ${params.min_cohorts} \\
        --out_prefix ${trait}_validated \\
        --report     ${trait}_validation_report.tsv
    """
}
