process EXTRACT_ALL_KEYWORD_GENES_MUTATED {
    tag "Extracting all genes with keyword from mutated results"
    publishDir "${params.outputDir}/extracted_genes_mutated", mode: 'copy'

    input:
    path gene_data_csv
    val keyword

    output:
    path "all_mutated_extracted_genes.fasta", emit: extracted_genes_mutated

    script:
    """
    chmod +x $projectDir/bin/extract_keyword_genes.py
    $projectDir/bin/extract_keyword_genes.py ${gene_data_csv} all_mutated_extracted_genes.fasta "${keyword}"
    """
}