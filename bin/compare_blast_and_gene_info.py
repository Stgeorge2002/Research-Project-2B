#!/usr/bin/env python3

import argparse
import csv

def parse_blast_results(blast_results_file):
    blast_data = []
    with open(blast_results_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            blast_data.append(row)
    return blast_data

def parse_gene_info(gene_info_file):
    gene_info = {}
    with open(gene_info_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            gene_info[row[0]] = {
                'header': row[1],
                'length': len(row[1].split()[-1])  # Assuming the sequence is the last field
            }
    return gene_info

def check_truncation(query_length, subject_length, alignment_length):
    if int(query_length) < int(subject_length):
        return "Query may be truncated (shorter than subject)"
    elif int(alignment_length) < int(subject_length):
        return "Possible partial alignment"
    else:
        return "No truncation detected"

def create_output(blast_data, gene_info, output_file):
    with open(output_file, 'w') as f:
        for row in blast_data:
            subject_id = row[1]
            query_length = row[4]
            alignment_length = row[3]
            if subject_id in gene_info:
                subject_length = gene_info[subject_id]['length']
                truncation_status = check_truncation(query_length, subject_length, alignment_length)
                
                f.write("BLAST Alignment:\n")
                f.write("\t".join(row) + "\n")
                f.write("Subject Sequence Header:\n")
                f.write(gene_info[subject_id]['header'] + "\n")
                f.write(f"Truncation Status: {truncation_status}\n")
                f.write(f"Query Length: {query_length}, Subject Length: {subject_length}, Alignment Length: {alignment_length}\n\n")
            else:
                f.write(f"Warning: No gene info found for subject ID {subject_id}\n\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare BLAST results with gene info and create output")
    parser.add_argument("--blast_results", required=True, help="BLAST results file")
    parser.add_argument("--gene_info", required=True, help="Gene info TSV file")
    parser.add_argument("--output", required=True, help="Output file")
    args = parser.parse_args()

    blast_data = parse_blast_results(args.blast_results)
    gene_info = parse_gene_info(args.gene_info)
    create_output(blast_data, gene_info, args.output)