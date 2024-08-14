process MUTATE_GENES {
    tag "Mutating genes in ${sampleName}"
    publishDir "${params.outputDir}/mutated_genes", mode: 'copy'
    
    container 'mutation_holder:latest'

    input:
    tuple val(sampleName), path(gff), path(fna)

    output:
    tuple val("${sampleName}_mutated"), path("${sampleName}_mutated.gff"), path("${sampleName}_mutated.fna"), emit: mutated_genome
    path "${sampleName}_combined_genes.fasta", emit: combined_genes
    path "${sampleName}_mutation_info.txt", emit: mutation_info

    script:
    """
    Rscript /usr/local/bin/mutate_genes.R \
        ${gff} \
        ${fna} \
        ${sampleName}_mutated.fna \
        ${sampleName}_combined_genes.fasta \
        ${sampleName}_mutation_info.txt \
        ${params.num_genes_to_mutate} \
        ${params.mutation_keyword} \
        ${params.mutation_mode}  # Pass the mode parameter

    cp ${gff} ${sampleName}_mutated.gff
    """
}