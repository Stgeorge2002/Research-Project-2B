process ANALYZE_NON_MUTATED_GENES {
    tag "Analyzing non-mutated genes"
    publishDir "${params.outputDir}/non_mutated_analysis", mode: 'copy'

    input:
    path fasta_file

    output:
    path "${fasta_file.simpleName}_non_mutated_analysis.txt", emit: analysis_results

    script:
    """
    python3 /app/analyze_premature_stops.py ${fasta_file} ${fasta_file.simpleName}_non_mutated_analysis.txt
    """
}