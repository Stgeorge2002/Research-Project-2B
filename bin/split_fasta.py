#!/usr/bin/env python3

import os
import sys
from Bio import SeqIO

def split_fasta(input_file, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    
    # Dictionary to store sequences for each prefix
    prefix_sequences = {}
    
    for record in SeqIO.parse(input_file, "fasta"):
        # Extract the prefix from the sequence identifier
        prefix = record.id.split('_')[1]
        
        if prefix not in prefix_sequences:
            prefix_sequences[prefix] = []
        
        prefix_sequences[prefix].append(record)
    
    # Write sequences for each prefix to separate files
    for prefix, sequences in prefix_sequences.items():
        output_file = os.path.join(output_dir, f"{prefix}_sequences.fasta")
        SeqIO.write(sequences, output_file, "fasta")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python split_fasta.py <input_fasta> <output_dir>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_dir = sys.argv[2]

    split_fasta(input_file, output_dir)
    print(f"FASTA file has been split into multiple files in the directory: {output_dir}")