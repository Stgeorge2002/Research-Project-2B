process ANALYZE_PREMATURE_STOPS_ASSEMBLED {
    tag "Analyzing premature stops in ${sampleName}"
    publishDir "${params.outputDir}/Analyze_Premature_Stops", mode: 'copy'

    input:
    tuple val(sampleName), path(assembled_genome)

    output:
    path "${sampleName}_premature_stops_analysis.txt", emit: analysis_results

    script:
    """
    python3 $projectDir/bin/analyze_premature_stops.py ${assembled_genome} ${sampleName}_premature_stops_analysis.txt
    """
}