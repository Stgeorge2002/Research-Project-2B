process UNICYCLER_ASSEMBLY {
    tag "Unicycler on ${sampleName}"
    container 'staphb/unicycler:latest'
    publishDir "${params.outputDir}/Unicycler_assembly", mode: 'copy'

    input:
    tuple val(sampleName), path(read1), path(read2)
    
    output:
    tuple val(sampleName), path("output/assembly.fasta"), emit: assembly
    path "output/unicycler.log", emit: log

    script:
    """
    unicycler \
        -1 ${read1} \
        -2 ${read2} \
        -o ./output \
        --threads ${task.cpus} \
        --mode ${params.unicycler_mode ?: 'bold'} \
        --min_fasta_length ${params.unicycler_min_fasta_length ?: 300} \
        --kmer_count ${params.unicycler_kmer_count ?: 4} \
        --min_component_size ${params.unicycler_min_component_size ?: 300} \
        --min_dead_end_size ${params.unicycler_min_dead_end_size ?: 300} \
        --verbosity ${params.unicycler_verbosity ?: 2} \
        --min_bridge_qual ${params.unicycler_min_bridge_qual ?: 3} \
        ${params.unicycler_keep_temp ? '--keep_temp' : ''} \
        ${params.unicycler_no_correct ? '--no_correct' : ''}
    """
}