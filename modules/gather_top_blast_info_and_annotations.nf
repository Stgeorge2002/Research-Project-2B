process GATHER_TOP_BLAST_INFO_AND_ANNOTATIONS {
    tag "Gathering top BLAST info and annotations for ${sampleName}"
    publishDir "${projectDir}/${params.outputDir}_XtraDetect_Output", mode: 'copy'

    input:
    tuple val(sampleName), path(blast_results_with_seq)
    tuple val(sampleName), path(extracted_genes)
    path pangenome_alignment_results
    path snippy_analysis
    path query_gff
    path pangenome_reference
    path gene_info

    output:
    tuple val(sampleName), path("${sampleName}_top_blast_info_and_annotations.tsv"), emit: top_blast_info_and_annotations

    script:
    """
    python3 ${projectDir}/bin/gather_top_blast_info_and_annotations.py ${sampleName} ${blast_results_with_seq} ${extracted_genes} ${pangenome_alignment_results} ${snippy_analysis} ${query_gff} ${params.email} ${pangenome_reference} ${gene_info}
    """
}