process SNIPPY_MUTATED {
    tag "SNIPPY_MUTATED on ${sampleName}"
    label 'process_high'
    publishDir "${params.outputDir}/snippy_output_mutated", mode: 'copy'

    input:
    tuple val(sampleName), path(forwardReads), path(reverseReads)
    path non_mutated_reference_fasta

    output:
    path "snippy_output/*", emit: snippy_results

    script:
    """
    set -e
    echo "Running Snippy command:"
    snippy --cpus ${task.cpus} \
           --outdir snippy_output \
           --ref ${non_mutated_reference_fasta} \
           --R1 ${forwardReads} \
           --R2 ${reverseReads} \
           --force

    echo "SNIPPY_MUTATED process completed for ${sampleName}"
    """
}