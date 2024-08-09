process CREATE_CUSTOM_ACETYLTRANSFERASE_DB {
    errorStrategy 'retry'
    maxRetries 3
    
    output:
    path "spneumo_acetyltransferase_db*", emit: db
    
    beforeScript "chmod +x $projectDir/bin/create_custom_acetyltransferase_db.sh"

    script:
    """
    chmod +x $projectDir/bin/create_custom_acetyltransferase_db.sh
    $projectDir/bin/create_custom_acetyltransferase_db.sh "${params.acetyltransferase_search_term}"
    """
}