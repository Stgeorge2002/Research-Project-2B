process SNIPPY_REAL_READS {
    tag "Snippy on ${sampleName}"
    publishDir "${params.outputDir}/snippy_real_reads", mode: 'copy'

    input:
    tuple val(sampleName), path(read1), path(read2)
    path reference_genome

    output:
    tuple val(sampleName), path("${sampleName}_snippy"), path(read1), path(read2), emit: snippy_results

    script:
    """
    snippy --cpus ${task.cpus} --outdir ${sampleName}_snippy --ref ${reference_genome} --R1 ${read1} --R2 ${read2}
    """
}