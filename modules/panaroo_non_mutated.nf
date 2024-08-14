process PANAROO_NON_MUTATED {
    tag "PANAROO on non-mutated genomes"
    publishDir "${params.outputDir}/panaroo_non_mutated", mode: 'copy'

    input:
    path gff_files

    output:
    path "panaroo_output", emit: panaroo_results

    script:
    """
    echo "Input GFF files:"
    ls -l ${gff_files}
    panaroo -i ${gff_files.join(' ')} \
            -o panaroo_output \
            --clean-mode ${params.panaroo_clean_mode} \
            --merge_paralogs \
            --threads ${task.cpus}
    """
}