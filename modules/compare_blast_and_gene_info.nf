// modules/compare_blast_and_gene_info.nf

process COMPARE_BLAST_AND_GENE_INFO {
    tag "Comparing BLAST results with gene info and reference alignments"
    publishDir "${projectDir}/PanG_DB_Output", mode: 'copy'

    input:
    path blast_results
    path gene_info
    path reference_alignments
    path gene_data_csv
    path pangenome_reference

    output:
    path "${blast_results.baseName}_with_gene_info_and_reference.txt", emit: comparison_results
    path "truncated_genes.txt", emit: truncated_genes

    script:
    """
    python3 ${projectDir}/bin/compare_blast_and_gene_info.py ${blast_results} ${gene_info} ${reference_alignments} ${gene_data_csv} ${pangenome_reference}
    """
}