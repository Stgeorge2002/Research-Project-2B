#!/usr/bin/env python3

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import hashlib

def shorten_header(header, max_length=36):
    if len(header) <= max_length:
        return header
    # Use MD5 hash to create a unique, shortened identifier
    hash_object = hashlib.md5(header.encode())
    return f"seq_{hash_object.hexdigest()[:8]}"

def parse_blast_results(blast_file):
    results = {}
    with open(blast_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            query_id = shorten_header(fields[0])
            subject_id = shorten_header(fields[1])
            query_seq = fields[12]
            subject_seq = fields[14]
            results[query_id] = {'subject_id': subject_id, 'query_seq': query_seq, 'subject_seq': subject_seq}
    return results

def main(blast_results, output_query, output_subject):
    results = parse_blast_results(blast_results)

    query_records = []
    subject_records = []

    for query_id, data in results.items():
        query_records.append(SeqRecord(Seq(data['query_seq']), id=query_id, description=""))
        subject_records.append(SeqRecord(Seq(data['subject_seq']), id=data['subject_id'], description=""))

    SeqIO.write(query_records, output_query, "fasta")
    SeqIO.write(subject_records, output_subject, "fasta")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract sequences from BLAST results")
    parser.add_argument("--blast_results", required=True, help="BLAST results file")
    parser.add_argument("--output_query", required=True, help="Output file for query sequences")
    parser.add_argument("--output_subject", required=True, help="Output file for subject sequences")
    args = parser.parse_args()

    main(args.blast_results, args.output_query, args.output_subject)