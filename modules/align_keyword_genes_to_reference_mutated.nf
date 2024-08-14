process ALIGN_KEYWORD_GENES_TO_REFERENCE_MUTATED {
    tag "Aligning individual mutated acetyltransferases to pangenome reference"
    publishDir "${params.outputDir}/aligned_acetyltransferases_mutated", mode: 'copy'

    input:
    path extracted_genes
    path pangenome_reference

    output:
    path "acetyltransferases_alignments_mutated.txt", emit: alignments_mutated

    script:
    """
    # Create a BLAST database from the pangenome reference
    makeblastdb -in ${pangenome_reference} -dbtype nucl -out pangenome_db

    # Align each extracted gene to the pangenome reference using BLAST
    blastn -query ${extracted_genes} \
           -db pangenome_db \
           -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq" \
           -evalue ${params.blast_evalue} \
           -word_size ${params.blast_word_size} \
           -perc_identity ${params.blast_perc_identity} \
           -max_target_seqs ${params.blast_max_target_seqs} \
           -out acetyltransferases_alignments_mutated.txt

    # Clean up
    rm pangenome_db*
    """
}