#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq

def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())

def is_start_codon(codon):
    return codon in ['ATG', 'GTG', 'TTG']

def is_stop_codon(codon):
    return codon in ['TAA', 'TAG', 'TGA']

def find_orf_stop_codons(seq):
    orf_codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
    stop_positions = [i * 3 for i, codon in enumerate(orf_codons) if is_stop_codon(codon)]
    return stop_positions, len(stop_positions)

def determine_orientation(seq):
    if is_start_codon(seq[:3]) and is_stop_codon(seq[-3:]):
        return '+', seq
    elif is_stop_codon(seq[:3]) and is_start_codon(reverse_complement(seq)[:3]):
        return '-', reverse_complement(seq)
    elif is_start_codon(seq[:3]):
        return '+', seq
    elif is_start_codon(reverse_complement(seq)[:3]):
        return '-', reverse_complement(seq)
    else:
        # If orientation can't be determined, assume forward
        return '+', seq

def parse_fasta_header(header):
    parts = header.split('_')
    if len(parts) < 2:
        return header, "Unknown", "Unknown product"
    gene_id = parts[0]
    gene_number = parts[1]
    gene_name = parts[2] if len(parts) > 2 else "Unknown"
    product = ' '.join(parts[3:]).replace('_', ' ') if len(parts) > 3 else "Unknown product"
    return f"{gene_id}_{gene_number}", gene_name, product

def analyze_fasta(file_path, output_file):
    with open(output_file, 'w') as out:
        for record in SeqIO.parse(file_path, "fasta"):
            seq = str(record.seq).upper()
            try:
                gene_id, gene_name, product = parse_fasta_header(record.id)
            except Exception as e:
                print(f"Error parsing header for sequence: {record.id}", file=sys.stderr)
                print(f"Error details: {str(e)}", file=sys.stderr)
                continue
            
            strand, oriented_seq = determine_orientation(seq)
            
            stop_positions, stop_count = find_orf_stop_codons(oriented_seq)
            
            # Check if the last stop codon is at the end of the sequence
            if stop_positions and stop_positions[-1] >= len(oriented_seq) - 3:
                premature_stops = stop_positions[:-1]
                premature_stop_count = len(premature_stops)
            else:
                premature_stops = stop_positions
                premature_stop_count = stop_count
            
            out.write(f"Gene ID: {gene_id}\n")
            out.write(f"Gene Name: {gene_name}\n")
            out.write(f"Product: {product}\n")
            out.write(f"Strand: {strand}\n")
            out.write(f"Sequence length: {len(seq)}\n")
            out.write(f"Total ORF stop codons: {stop_count}\n")
            out.write(f"Premature stop codons: {premature_stop_count}\n")
            out.write(f"ORF stop codon positions: {', '.join(map(str, stop_positions))}\n")
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