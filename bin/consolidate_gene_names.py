#!/usr/bin/env python3

from Bio import SeqIO
from collections import defaultdict
import sys

def consolidate_gene_names(input_file, output_file):
    gene_families = defaultdict(list)
    
    # First pass: group sequences by gene name (without suffixes)
    for record in SeqIO.parse(input_file, "fasta"):
        gene_name = record.id.split('_')[0]  # Assume suffix is after first underscore
        gene_families[gene_name].append(record)
    
    # Second pass: write consolidated sequences
    with open(output_file, 'w') as output_handle:
        for gene_name, records in gene_families.items():
            # Sort records by length, descending
            records.sort(key=lambda r: len(r.seq), reverse=True)
            
            # Write the longest sequence with the consolidated name
            consolidated_record = records[0]
            consolidated_record.id = gene_name
            consolidated_record.description = ""  # Remove the description
            SeqIO.write(consolidated_record, output_handle, "fasta")

    print(f"Consolidated gene names from {input_file} to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python consolidate_gene_names.py <input_file> <output_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    consolidate_gene_names(input_file, output_file)