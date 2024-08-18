#!/usr/bin/env python3
import csv
import sys
from Bio import SeqIO

def search_genes(keyword, gene_descriptions):
    matching_genes = []
    for gene, description in gene_descriptions.items():
        if keyword.lower() in description.lower():
            matching_genes.append(gene)
    return matching_genes

def extract_sequences(genes, gene_sequences):
    sequences = []
    for gene in genes:
        if gene in gene_sequences:
            sequences.append(gene_sequences[gene])
    return sequences

def main(panaroo_output, keyword, output_file):
    # Step 1: Create an index of genes from pan_genome_reference.fa
    gene_sequences = SeqIO.to_dict(SeqIO.parse(f"{panaroo_output}/pan_genome_reference.fa", "fasta"))

    # Step 2: Map gene names to descriptions from gene_data.csv
    gene_descriptions = {}
    with open(f"{panaroo_output}/gene_data.csv", "r") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            gene_name = row["gene_name"] if row["gene_name"] else row["annotation_id"]
            gene_descriptions[gene_name] = row["description"]

    # Extract matching sequences
    matching_genes = search_genes(keyword, gene_descriptions)
    matching_sequences = extract_sequences(matching_genes, gene_sequences)

    # Save matching sequences to a file
    with open(output_file, "w") as out_handle:
        SeqIO.write(matching_sequences, out_handle, "fasta")

    print(f"Found {len(matching_sequences)} genes matching '{keyword}'.")
    print(f"Sequences have been saved to {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: extract_keyword_genes.py <panaroo_output_dir> <keyword> <output_file>")
        sys.exit(1)
    
    panaroo_output = sys.argv[1]
    keyword = sys.argv[2]
    output_file = sys.argv[3]
    main(panaroo_output, keyword, output_file)