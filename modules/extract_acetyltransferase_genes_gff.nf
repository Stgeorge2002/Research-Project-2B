process EXTRACT_ACETYLTRANSFERASE_GENES_GFF {
    tag "Extract genes from ${sampleName}"
    publishDir "${params.outputDir}/extracted_genes_${task.process.tokenize(':')[-1].toLowerCase()}", mode: 'copy'

    input:
    tuple val(sampleName), path(gff), path(fna)
    val(keyword)

    output:
    tuple val(sampleName), path("${sampleName}_extracted_genes.fasta"), emit: extracted_genes

    script:
    """
    chmod +x $projectDir/bin/extract_genes_from_gff.py
    python3 $projectDir/bin/extract_genes_from_gff.py ${gff} ${fna} ${keyword} ${sampleName}_extracted_genes.fasta
    """
}