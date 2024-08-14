process COMBINE_REAL_RESULTS {
    tag "Combining results for ${sampleName}"
    publishDir "${params.outputDir}/Final_results_Combined_Reads", mode: 'copy'

    input:
    tuple val(sampleName), path(snippy_analysis)
    tuple val(sampleName), path(blast_results)
    tuple val(sampleName), path(extracted_genes)

    output:
    tuple val(sampleName), path("${sampleName}_combined_analysis.txt"), emit: combined_results

    script:
    """
    python3 $projectDir/bin/combine_real_results.py \
        --snippy ${snippy_analysis} \
        --blast ${blast_results} \
        --fasta ${extracted_genes} \
        --output ${sampleName}_combined_analysis.txt
    """
}