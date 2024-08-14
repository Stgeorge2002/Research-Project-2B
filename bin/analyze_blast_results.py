#!/usr/bin/env python3
import csv
from Bio import SeqIO, Entrez
from collections import defaultdict
import sys
import os
import re

# Set your email for NCBI queries
Entrez.email = "your_email@example.com"

def log_error(message):
    print(f"ERROR: {message}", file=sys.stderr)

def extract_info_from_header(header):
    gene_name = "Unknown"
    product = "Unknown"
    parts = header.split()
    for part in parts:
        if part.startswith("gene:"):
            gene_name = part.split(":")[1]
        elif part.startswith("product:"):
            product = " ".join(parts[parts.index(part)+1:])
            break
    return gene_name, product

def extract_genbank_id(blast_hit):
    match = re.search(r'cds_(\w+\.\d+)', blast_hit)
    return match.group(1) if match else None

def fetch_genbank_data(genbank_id):
    handle = Entrez.efetch(db="protein", id=genbank_id, rettype="gb", retmode="text")
    record = SeqIO.read(handle, "genbank")
    handle.close()
    return record

def get_gene_name(record):
    for feature in record.features:
        if feature.type == "CDS":
            gene_qualifiers = feature.qualifiers.get("gene", [])
            if gene_qualifiers:
                return gene_qualifiers[0]
            locus_tag = feature.qualifiers.get("locus_tag", [])
            if locus_tag:
                return locus_tag[0]
    return "Unknown"

try:
    # Get command-line arguments
    blast_results_file = sys.argv[1]
    query_fasta_file = sys.argv[2]
    sample_id = sys.argv[3]
    output_file = f"{sample_id}_blast_analysis.txt"

    with open(output_file, 'w') as outfile:
        outfile.write(f"Processing sample: {sample_id}\n")
        outfile.write(f"BLAST results file: {blast_results_file}\n")
        outfile.write(f"Query FASTA file: {query_fasta_file}\n\n")

        # Check if files exist
        if not os.path.exists(blast_results_file):
            log_error(f"BLAST results file not found: {blast_results_file}")
            sys.exit(1)
        if not os.path.exists(query_fasta_file):
            log_error(f"Query FASTA file not found: {query_fasta_file}")
            sys.exit(1)

        # Read query sequences
        query_sequences = SeqIO.to_dict(SeqIO.parse(query_fasta_file, "fasta"))

        # Read BLAST results
        blast_results = defaultdict(list)
        with open(blast_results_file, 'r') as f:
            reader = csv.reader(f, delimiter='\t')
            for row in reader:
                qseqid, sseqid, pident, length, qlen, slen = row[:6]
                blast_results[qseqid].append((sseqid, float(pident), int(length), int(qlen), int(slen)))

        # Analyze results
        for qseqid, matches in blast_results.items():
            outfile.write(f"Reference: {qseqid}\n")
            
            # Extract reference information
            if qseqid in query_sequences:
                ref_seq = query_sequences[qseqid]
                ref_gene, ref_product = extract_info_from_header(ref_seq.description)
                outfile.write(f"Reference Gene: {ref_gene}\n")
                outfile.write(f"Reference Product: {ref_product}\n")
            
            # Select the best match based on percent identity and coverage
            best_match = max(matches, key=lambda x: (x[1], x[2]/x[3]))
            sseqid, pident, length, qlen, slen = best_match
            
            outfile.write(f"  Best Database Match: {sseqid}\n")
            outfile.write(f"  Percent identity: {pident:.2f}%\n")
            outfile.write(f"  Alignment length: {length}\n")
            outfile.write(f"  Reference length: {qlen}, Database Match length: {slen}\n")
            
            # Check for truncation
            if qlen < slen:
                outfile.write(f"  Truncation: Yes (Reference is shorter by {slen - qlen} bp)\n")
            else:
                outfile.write("  Truncation: No\n")
            
            # Analyze alignment coverage
            ref_coverage = (length / qlen) * 100
            db_coverage = (length / slen) * 100
            outfile.write(f"  Reference coverage: {ref_coverage:.2f}%\n")
            outfile.write(f"  Database Match coverage: {db_coverage:.2f}%\n")
            
            # Fetch GenBank data for the best database match
            genbank_id = extract_genbank_id(sseqid)
            if genbank_id:
                try:
                    record = fetch_genbank_data(genbank_id)
                    db_gene_name = get_gene_name(record)
                    db_description = record.description
                    outfile.write(f"  Database Gene Name: {db_gene_name}\n")
                    outfile.write(f"  Database Description: {db_description}\n")
                except Exception as e:
                    outfile.write(f"  Error fetching GenBank data: {str(e)}\n")
            
            outfile.write("\n")

        outfile.write(f"Analysis complete for sample: {sample_id}\n")

    print(f"Analysis complete. Results written to {output_file}")

except Exception as e:
    log_error(f"An error occurred: {str(e)}")
    sys.exit(1)