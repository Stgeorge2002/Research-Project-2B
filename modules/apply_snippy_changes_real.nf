process APPLY_SNIPPY_CHANGES_REAL {
    tag "Applying Snippy changes and analyzing premature stops for ${sampleName}"
    publishDir "${params.outputDir}/modified_genes_real/${sampleName}", mode: 'copy'

    input:
    tuple val(sampleName), path(snippy_output)

    output:
    tuple val(sampleName), path("modified_reference.fasta"), emit: modified_fasta
    tuple val(sampleName), path("modified_reference_analysis.txt"), emit: analysis_results

    script:
    """
    # Use the assumed file structure
    SNPS_FILE="${snippy_output}/${sampleName}.csv"
    REF_FILE="${snippy_output}/reference/ref.fa"

    python ${projectDir}/bin/apply_snippy_changes.py \
        \$SNPS_FILE \
        \$REF_FILE \
        modified_reference.fasta

    # Analyze premature stops in the modified reference
    python ${projectDir}/bin/analyze_premature_stops.py \
        modified_reference.fasta \
        modified_reference_analysis.txt
    """
}