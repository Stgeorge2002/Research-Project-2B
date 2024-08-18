process EXTRACT_TOP_BLAST_HITS {
    tag "Extracting top BLAST hits for ${sampleName}"
    publishDir "${params.outputDir}/top_blast_hits", mode: 'copy'

    input:
    tuple val(sampleName), path(blast_results_with_seq)
    tuple val(sampleName), path(query_fasta)

    output:
    tuple val(sampleName), path("${sampleName}_top_hits_query.fasta"), emit: query_fasta
    tuple val(sampleName), path("${sampleName}_top_hits_subject.fasta"), emit: subject_fasta

    script:
    """
    python3 ${projectDir}/bin/extract_top_blast_hits.py \
        --blast_results ${blast_results_with_seq} \
        --query_fasta ${query_fasta} \
        --output_query ${sampleName}_top_hits_query.fasta \
        --output_subject ${sampleName}_top_hits_subject.fasta
    """
}