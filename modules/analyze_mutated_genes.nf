process ANALYZE_MUTATED_GENES {
    tag "Analyzing mutated genes"
    publishDir "${params.outputDir}/mutated_analysis", mode: 'copy'

    input:
    path fasta_file

    output:
    path "${fasta_file.simpleName}_mutated_analysis.txt", emit: analysis_results

    script:
    """
    python3 $projectDir/bin/analyze_premature_stops.py ${fasta_file} ${fasta_file.simpleName}_mutated_analysis.txt
    """
}