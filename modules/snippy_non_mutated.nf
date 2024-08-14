process SNIPPY_NON_MUTATED {
    tag "SNIPPY_NON_MUTATED on ${sampleName}"
    label 'process_high'
    publishDir "${params.outputDir}/snippy_output_non_mutated", mode: 'copy'

    input:
    tuple val(sampleName), path(forwardReads), path(reverseReads)
    path mutated_reference_fasta

    output:
    tuple val(sampleName), path("snippy_output"), emit: snippy_results

    script:
    """
    set -e
    echo "Running Snippy command:"
    snippy --cpus ${task.cpus} \
           --outdir snippy_output \
           --ref ${mutated_reference_fasta} \
           --R1 ${forwardReads} \
           --R2 ${reverseReads} \
           --force
           
    echo "SNIPPY_NON_MUTATED process completed for ${sampleName}"
    """
}