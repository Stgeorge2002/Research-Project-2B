process EXTRACT_NON_MUTATED_GENOME_KEYWORD_GENES {
    tag "Extracting genes with keyword from non-mutated results"
    publishDir "${params.outputDir}/extracted_genes", mode: 'copy'

    input:
    tuple val(sampleName), path(gene_data_csv)
    val keyword

    output:
    path "${sampleName}_non_mutated_extracted_genes.fasta", emit: extracted_genes

    script:
    """
    chmod +x $projectDir/bin/extract_keyword_genes_test.py
    $projectDir/bin/extract_keyword_genes_test.py ${gene_data_csv} ${sampleName}_non_mutated_extracted_genes.fasta ${sampleName} "${keyword}"
    """
}