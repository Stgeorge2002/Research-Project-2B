process DWGSIM_MUTATED {
    tag "DWGSIM on ${sampleName}"
    publishDir "${params.outputDir}/simulated_reads", mode: 'copy'

    input:
    tuple val(sampleName), path(fna_file)

    output:
    tuple val(sampleName), path("${sampleName}.bwa.read1.fastq.gz"), path("${sampleName}.bwa.read2.fastq.gz"), emit: simulated_reads

    script:
    """
    dwgsim \
        -e ${params.dwgsimErrorRate} \
        -E ${params.dwgsimErrorRate} \
        -N ${params.dwgsimNumReads} \
        -1 ${params.dwgsimReadLength} \
        -2 ${params.dwgsimReadLength} \
        -d ${params.dwgsimOuterDistance} \
        -s ${params.dwgsimStdDev} \
        -r 0 \
        -R 0 \
        -X 0 \
        "${fna_file}" \
        "${sampleName}"
    """
}