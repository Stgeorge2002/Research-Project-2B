#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
from Bio import SeqIO
from Bio.Seq import Seq

def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())

def split_into_codons(seq):
    return [seq[i:i+3] for i in range(0, len(seq) - 2, 3)]

def is_stop_codon(codon):
    return codon in ['TAA', 'TAG', 'TGA']

def find_premature_stops(seq, strand):
    seq = seq.upper()
    if strand == '-':
        seq = reverse_complement(seq)
    
    codons = split_into_codons(seq)
    premature_stops = [i * 3 for i, codon in enumerate(codons[:-1]) if is_stop_codon(codon)]
    
    if strand == '-':
        seq_length = len(seq)
        premature_stops = [seq_length - pos - 3 for pos in premature_stops]
    
    return premature_stops

def analyze_fasta(file_path, output_file):
    with open(output_file, 'w') as out:
        for record in SeqIO.parse(file_path, "fasta"):
            seq = str(record.seq)
            gene_name = record.id.split('_')[0]
            
            strand = '+' if '_original' in record.id else '-'
            
            premature_stops = find_premature_stops(seq, strand)
            
            out.write(f"Gene: {gene_name}\n")
            out.write(f"Strand: {strand}\n")
            out.write(f"Sequence length: {len(seq)}\n")
            out.write(f"Premature stop codons: {', '.join(map(str, premature_stops))}\n")
            out.write("\n")

def main():
    if len(sys.argv) != 3:
        print("Usage: python analyze_premature_stops.py <input_fasta> <output_txt>", file=sys.stderr)
        sys.exit(1)
    
    input_fasta = sys.argv[1]
    output_txt = sys.argv[2]
    
    try:
        analyze_fasta(input_fasta, output_txt)
        print(f"Analysis complete. Results written to {output_txt}", file=sys.stderr)
    except Exception as e:
        print(f"Error: {str(e)}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()