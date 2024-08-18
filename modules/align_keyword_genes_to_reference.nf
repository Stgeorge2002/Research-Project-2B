process ALIGN_KEYWORD_GENES_TO_REFERENCE {
    tag "Aligning individual acetyltransferases to pangenome reference"
    publishDir "${params.outputDir}/aligned_acetyltransferases", mode: 'copy'

    input:
    path extracted_genes
    path pangenome_reference

    output:
    path "acetyltransferases_alignments_with_seq.txt", emit: blast_results

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
           -out acetyltransferases_alignments.txt

    # Process BLAST results and add corresponding sequences
    awk -F'\\t' 'NR==FNR {
        if (\$0 ~ /^>/) {
            header = substr(\$0, 2);
        } else {
            seq[header] = seq[header] \$0;
        }
        next
    }
    {
        if (\$1 in seq) {
            print \$0 "\\t" seq[\$1]
        } else {
            print \$0 "\\tSequence not found"
        }
    }' ${extracted_genes} acetyltransferases_alignments.txt > acetyltransferases_alignments_with_seq.txt

    # Clean up
    rm pangenome_db*
    """
}