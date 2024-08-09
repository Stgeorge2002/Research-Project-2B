process PROKKA_INITIAL {
    tag "PROKKA on ${sampleName}"
    publishDir "${params.outputDir}/prokka_initial", mode: 'copy'

    input:
    tuple val(sampleName), path(fasta)

    output:
    tuple val(sampleName), path("${sampleName}/${sampleName}.gff"), path("${sampleName}/${sampleName}.fna"), emit: prokka_results
    path "${sampleName}/${sampleName}.gff", emit: gff_files
    path "${sampleName}/${sampleName}.fna", emit: fna_files

    script:
    """
    prokka --outdir ${sampleName} \
           --prefix ${sampleName} \
           --cpus ${task.cpus} \
           --genus ${params.prokka_genus} \
           --species ${params.prokka_species} \
           --strain ${params.prokka_strain} \
           ${fasta}
    """
}