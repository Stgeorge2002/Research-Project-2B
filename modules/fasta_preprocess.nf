process FASTA_PREPROCESS {
    tag "Preprocess ${sampleName}"
    publishDir "${params.outputDir}/fasta_preprocess", mode: 'copy'
    
    container 'staphb/seqkit:latest'

    input:
    tuple val(sampleName), path(fasta)

    output:
    tuple val(sampleName), path("${sampleName}_cleaned.fa"), emit: cleaned_fasta
    tuple val(sampleName), path("${sampleName}_qc_report.txt"), emit: qc_report

    script:
    def remove_duplicates = params.fasta_preprocess.remove_duplicates ? "| seqkit rmdup -s" : ""
    """
    # QC report
    seqkit stats ${fasta} > ${sampleName}_qc_report.txt

    # Clean and preprocess
    seqkit seq -m ${params.fasta_preprocess.min_length} ${fasta} \
    ${remove_duplicates} \
    -o ${sampleName}_cleaned.fa

    # Update QC report with cleaned data
    seqkit stats ${sampleName}_cleaned.fa >> ${sampleName}_qc_report.txt
    """
}