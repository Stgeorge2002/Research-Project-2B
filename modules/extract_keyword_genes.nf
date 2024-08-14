process EXTRACT_KEYWORD_GENES {
    tag "Extracting genes matching ${params.keyword}"
    publishDir "${params.outputDir}/extracted_genes", mode: 'copy'

    input:
    path panaroo_output
    val keyword

    output:
    path "matching_genes.fasta", emit: extracted_genes

    script:
    """
    python3 $projectDir/bin/extract_keyword_genes_real.py ${panaroo_output} ${keyword} matching_genes.fasta
    """
}