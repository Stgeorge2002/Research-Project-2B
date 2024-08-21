process SEARCH_TRUNCATED_GFF {
    publishDir "${projectDir}/${params.outputDir}_further_analysis", mode: 'copy'

    input:
    path truncated_genes
    path cleaned_fasta_files

    output:
    path "truncated_*.fa", emit: truncated_fa_files, optional: true
    path "search_log.txt", emit: search_log

    script:
    """
    chmod +x ${projectDir}/bin/search_truncated_gff.sh
    ${projectDir}/bin/search_truncated_gff.sh ${truncated_genes} ${cleaned_fasta_files}
    """
}