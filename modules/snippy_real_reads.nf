process SNIPPY_REAL_READS {
    tag "Snippy on ${sampleName}"
    publishDir "${params.outputDir}/snippy_results/${sampleName}", mode: 'copy'

    input:
    tuple val(sampleName), path(read1), path(read2), path(reference_fasta)

    output:
    tuple val(sampleName), path("${sampleName}_snippy_results"), emit: snippy_results

    script:
    """
    snippy --cpus ${task.cpus} \
           --outdir ${sampleName}_snippy_results \
           --ref ${reference_fasta} \
           --R1 ${read1} \
           --R2 ${read2} \
           --prefix ${sampleName}
    """
}