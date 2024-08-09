process SPADES_ASSEMBLY {
    tag "SPAdes on ${sampleName}"
    publishDir "${params.outputDir}/spades_assembly", mode: 'copy'

    input:
    tuple val(sampleName), path(snippy_dir), path(read1), path(read2)

    output:
    tuple val(sampleName), path("${sampleName}_assembly.fasta"), emit: assembled_genome
    path "${sampleName}_spades_log.txt"

    script:
    """
    spades.py -1 ${read1} -2 ${read2} \
              -o ./${sampleName}_spades \
              --isolate \
              --cov-cutoff auto \
              -t ${task.cpus}
    
    mv ./${sampleName}_spades/scaffolds.fasta ${sampleName}_assembly.fasta
    mv ./${sampleName}_spades/spades.log ${sampleName}_spades_log.txt
    """
}