process FASTP {
    tag "Preprocess ${sampleName}"
    publishDir "${params.outputDir}/fastp", mode: 'copy'

    container 'staphb/fastp:latest'

    input:
    tuple val(sampleName), path(read1), path(read2)

    output:
    tuple val(sampleName), path("${sampleName}_cleaned_R1.fastq.gz"), path("${sampleName}_cleaned_R2.fastq.gz"), emit: cleaned_reads
    path "${sampleName}_fastp.json", emit: json_report
    path "${sampleName}_fastp.html", emit: html_report

    script:
    """
    fastp -i ${read1} \
          -I ${read2} \
          -o ${sampleName}_cleaned_R1.fastq.gz \
          -O ${sampleName}_cleaned_R2.fastq.gz \
          -j ${sampleName}_fastp.json -h ${sampleName}_fastp.html \
          -q ${params.fastp.qualified_quality_phred} \
          -u ${params.fastp.unqualified_percent_limit} \
          -5 ${params.fastp.cut_mean_quality} \
          -W ${params.fastp.cut_window_size} \
          ${params.fastp.cut_front ? '--cut_front' : ''} \
          ${params.fastp.cut_tail ? '--cut_tail' : ''} \
          --cut_front_window_size ${params.fastp.cut_front_window_size} \
          --cut_front_mean_quality ${params.fastp.cut_front_mean_quality} \
          --cut_tail_window_size ${params.fastp.cut_tail_window_size} \
          --cut_tail_mean_quality ${params.fastp.cut_tail_mean_quality} \
          -l ${params.fastp.length_required} \
          -w ${task.cpus}
    """
}