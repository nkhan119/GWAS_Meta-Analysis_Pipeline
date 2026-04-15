// ================================================================
// Module: ANNOTATION
// Positional annotation of top loci:
//   - Nearest gene, consequence class, GWAS Catalog lookup
// Input : (trait, type, path(top_loci))
// Output: annotated → (trait, path)
// ================================================================

process ANNOTATION {

    tag "${trait}"
    label 'small'
    publishDir "${params.out_dir}/annotated/${trait}", mode: 'copy', overwrite: true

    input:
    tuple val(trait), val(type), path(top_loci)

    output:
    tuple val(trait), path("${trait}_annotated_loci.tsv"), emit: annotated

    script:
    """
    annotate_loci.py \\
        --input       ${top_loci} \\
        --trait       ${trait} \\
        --build       ${params.ref_build} \\
        --out         ${trait}_annotated_loci.tsv
    """
}
