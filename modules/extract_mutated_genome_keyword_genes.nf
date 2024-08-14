process EXTRACT_MUTATED_GENOME_KEYWORD_GENES {
    tag "Extracting genes with keyword from mutated results"
    publishDir "${params.outputDir}/extracted_genes", mode: 'copy'

    input:
    tuple val(sampleName), path(gene_data_csv)
    val keyword

    output:
    path "${sampleName}_mutated_extracted_genes.fasta", emit: extracted_genes

    beforeScript "chmod +x $projectDir/bin/extract_keyword_genes.py"

    script:
    """
    chmod +x $projectDir/bin/extract_keyword_genes_test.py
    $projectDir/bin/extract_keyword_genes_test.py ${gene_data_csv} ${sampleName}_mutated_extracted_genes.fasta ${sampleName} "${keyword}"
    """
}