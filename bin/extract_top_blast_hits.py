#!/usr/bin/env python3
import argparse
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parse_blast_results(blast_results_file):
    top_hits = {}
    with open(blast_results_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            qseqid, sseqid = row[0], row[1]
            qstart, qend, sstart, send = map(int, row[8:12])
            evalue = float(row[12])
            qseq, sseq = row[14], row[16]  # Updated column indices
            
            if qseqid not in top_hits or evalue < top_hits[qseqid]['evalue']:
                top_hits[qseqid] = {
                    'sseqid': sseqid,
                    'qseq': qseq,
                    'sseq': sseq,
                    'qstart': qstart,
                    'qend': qend,
                    'sstart': sstart,
                    'send': send,
                    'evalue': evalue
                }
    return top_hits

def extract_top_hits(blast_results_file, query_fasta_file, output_query_file, output_subject_file):
    top_hits = parse_blast_results(blast_results_file)
    
    # Write query sequences
    query_records = []
    original_queries = SeqIO.to_dict(SeqIO.parse(query_fasta_file, "fasta"))
    
    for i, (qseqid, hit_info) in enumerate(top_hits.items(), 1):
        full_seq = original_queries[qseqid].seq
        record = SeqRecord(full_seq, 
                           id=f"query_{i}", 
                           description=f"{qseqid} | E-value: {hit_info['evalue']} | Alignment: {hit_info['qstart']}-{hit_info['qend']}")
        query_records.append(record)
    
    SeqIO.write(query_records, output_query_file, "fasta")

    # Write subject sequences
    subject_records = []
    for i, (qseqid, hit_info) in enumerate(top_hits.items(), 1):
        record = SeqRecord(Seq(hit_info['sseq']), 
                           id=f"subject_{i}", 
                           description=f"{hit_info['sseqid']} | E-value: {hit_info['evalue']} | Alignment: {hit_info['sstart']}-{hit_info['send']}")
        subject_records.append(record)
    
    SeqIO.write(subject_records, output_subject_file, "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract top BLAST hits")
    parser.add_argument("--blast_results", required=True, help="BLAST results file")
    parser.add_argument("--query_fasta", required=True, help="Query FASTA file")
    parser.add_argument("--output_query", required=True, help="Output query FASTA file")
    parser.add_argument("--output_subject", required=True, help="Output subject FASTA file")
    args = parser.parse_args()

    extract_top_hits(args.blast_results, args.query_fasta, args.output_query, args.output_subject)