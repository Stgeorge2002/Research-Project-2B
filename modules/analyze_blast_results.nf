// modules/analyze_blast_results.nf
process ANALYZE_BLAST_RESULTS {
    tag "Analyzing BLAST results for ${sampleName}"
    publishDir "${params.outputDir}/blast_analysis", mode: 'copy'

    input:
    tuple val(sampleName), path(blast_results), path(query_fasta)
    path blast_db

    output:
    tuple val(sampleName), path("${sampleName}_blast_analysis.txt"), emit: analysis_results

    script:
    """
    python3 ${workflow.projectDir}/bin/analyze_blast_results.py ${blast_results} ${query_fasta} ${sampleName}
    """
}