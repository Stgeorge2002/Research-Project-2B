// modules/blast_acetyltransferases_non_mutated.nf
process BLAST_ACETYLTRANSFERASES_NON_MUTATED {
    tag "BLASTing non-mutated acetyltransferases"
    publishDir "${params.outputDir}/blast_results/non_mutated", mode: 'copy'

    input:
    path fasta_file
    path blast_db_files

    output:
    path "${fasta_file.baseName}_blast_results.tsv", emit: blast_results

    script:
    def db_name = blast_db_files.first().baseName
    """
    blastn -task blastn-short \
           -query ${fasta_file} \
           -db ${db_name} \
           -outfmt "6 qseqid sseqid pident length qlen slen mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq" \
           -max_target_seqs ${params.blast_max_target_seqs} \
           -word_size ${params.blast_word_size} \
           -evalue ${params.blast_evalue} \
           -perc_identity ${params.blast_perc_identity} \
           -dust ${params.blast_dust} \
           -soft_masking ${params.blast_soft_masking} \
           -xdrop_ungap ${params.blast_xdrop_ungap} \
           -xdrop_gap ${params.blast_xdrop_gap} \
           -ungapped \
           -num_threads ${task.cpus} \
           -out ${fasta_file.baseName}_blast_results.tsv
    """
}