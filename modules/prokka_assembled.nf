process PROKKA_ASSEMBLED {
    tag "PROKKA on ${sampleName}"
    publishDir "${params.outputDir}/prokka_assembled", mode: 'copy'

    input:
    tuple val(sampleName), path(gap_filled_scaffolds)

    output:
    tuple val(sampleName), path("${sampleName}/${sampleName}.gff"), emit: gff_file
    tuple val(sampleName), path("${sampleName}/${sampleName}.fna"), emit: fna_file
    path "${sampleName}/*"

    script:
    """
    prokka --cpus ${task.cpus} \
           --prefix ${sampleName} \
           --outdir ${sampleName} \
           --genus ${params.genus} \
           --species ${params.species} \
           --strain ${params.strain} \
           ${gap_filled_scaffolds}
    """
}