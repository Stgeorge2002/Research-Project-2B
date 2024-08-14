process BLAST_ACETYLTRANSFERASES {
    tag "BLASTing acetyltransferases for ${sampleName}"
    publishDir "${params.outputDir}/blast_results/${task.process.tokenize(':')[-1].toLowerCase()}", mode: 'copy'

    input:
    tuple val(sampleName), path(fasta_file)
    path blast_db_files

    output:
    tuple val(sampleName), path("${sampleName}_blast_results.tsv"), emit: blast_results

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
           -out ${sampleName}_blast_results.tsv
    """
}