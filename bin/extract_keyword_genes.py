#!/usr/bin/env python3
import csv
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def clean_string(s):
    return ''.join(c if c.isalnum() else '_' for c in s)

def extract_genes(input_file, output_file, keyword):
    records = []
    with open(input_file, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if keyword.lower() in row['description'].lower():
                seq = Seq(row['dna_sequence'])
                sample_id = row['scaffold_name'].split('.')[0]  # Extract sample identifier
                unique_id = f"{sample_id}_{clean_string(row['annotation_id'])}_{clean_string(row['gene_name'])}_{clean_string(row['description'][:50])}"
                record = SeqRecord(seq, id=unique_id, description="")
                records.append(record)
    
    SeqIO.write(records, output_file, "fasta")
    print(f"Extracted {len(records)} genes with keyword '{keyword}' to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python extract_keyword_genes.py <input_file> <output_file> <keyword>")
        sys.exit(1)
    
    input_file, output_file, keyword = sys.argv[1:]
    extract_genes(input_file, output_file, keyword)