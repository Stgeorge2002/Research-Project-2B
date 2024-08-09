process EXTRACT_ALL_KEYWORD_GENES_NON_MUTATED {
    tag "Extracting all genes with keyword from non-mutated results"
    publishDir "${params.outputDir}/extracted_genes", mode: 'copy'

    input:
    path gene_data_csv
    val keyword

    output:
    path "all_non_mutated_extracted_genes.fasta", emit: extracted_genes

    script:
    """
    chmod +x $projectDir/bin/extract_keyword_genes.py
    $projectDir/bin/extract_keyword_genes.py ${gene_data_csv} all_non_mutated_extracted_genes.fasta "${keyword}"
    """
}