process CREATE_CUSTOM_ACETYLTRANSFERASE_DB {
    errorStrategy 'retry'
    maxRetries 1
    publishDir "${params.outputDir}/Custom_BLAST_DB", mode: 'copy'

    output:
    path "spneumo_custom_db*", emit: db
    path "gene_info.tsv", emit: gene_info
    path "filtered_sequences.fasta", emit: sequences
    
    beforeScript "chmod +x $projectDir/bin/create_custom_acetyltransferase_db.sh"

    script:
    """
    chmod +x $projectDir/bin/create_custom_acetyltransferase_db.sh
    $projectDir/bin/create_custom_acetyltransferase_db.sh "${params.search_term}" "${params.acetyltransferase_search_term}" "${params.acetyltransferase_search_retmax}"
    """
}