#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import csv
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os

def parse_snippy_output(snippy_file):
    changes = {}
    with open(snippy_file, 'r') as csvfile:
        reader = csv.reader(csvfile)
        header = next(reader)  # Skip the header row
        for row in reader:
            if len(row) < 5:  # Ensure the row has at least 5 columns
                continue
            chrom = row[0]
            pos = int(row[1])
            ref = row[3]
            alt = row[4]
            if chrom not in changes:
                changes[chrom] = []
            changes[chrom].append((pos, ref, alt))
    return changes

def apply_changes(seq, changes):
    seq_list = list(str(seq))
    for pos, ref, alt in sorted(changes, key=lambda x: x[0], reverse=True):
        pos = pos - 1  # Convert to 0-based index
        if seq_list[pos:pos+len(ref)] == list(ref):
            seq_list[pos:pos+len(ref)] = alt
        else:
            print(f"Warning: Expected '{ref}' at position {pos+1}, but found '{''.join(seq_list[pos:pos+len(ref)])}'")
    return ''.join(seq_list)

def main():
    if len(sys.argv) != 4:
        print("Usage: python apply_snippy_changes.py <snippy_snps.csv> <reference.fa> <output.fasta>")
        sys.exit(1)

    snippy_file = sys.argv[1]
    reference_fasta = sys.argv[2]
    output_fasta = sys.argv[3]

    changes = parse_snippy_output(snippy_file)
    
    modified_records = []
    for record in SeqIO.parse(reference_fasta, "fasta"):
        if record.id in changes:
            modified_seq = apply_changes(record.seq, changes[record.id])
            modified_record = SeqRecord(Seq(modified_seq), id=record.id, description=record.description)
            modified_records.append(modified_record)
        else:
            modified_records.append(record)
    
    # Ensure the output directory exists (if there's a directory component)
    output_dir = os.path.dirname(output_fasta)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    
    SeqIO.write(modified_records, output_fasta, "fasta")
    
    print(f"Modified sequences written to {output_fasta}")

if __name__ == "__main__":
    main()