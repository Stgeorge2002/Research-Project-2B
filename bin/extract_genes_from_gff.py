#!/usr/bin/env python3
import sys
from Bio import SeqIO
import re
import gffutils

def extract_acetyltransferase_genes(gff_file, fasta_file, keyword, output_file):
    # Create a database from the GFF file
    db = gffutils.create_db(gff_file, dbfn=':memory:', force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)

    acetyltransferase_genes = []

    # Iterate through the features in the database
    for feature in db.all_features(featuretype='CDS'):
        product = feature.attributes.get('product', [''])[0]
        if re.search(keyword, product, re.IGNORECASE):
            gene_id = feature.attributes.get('locus_tag', ['Unknown'])[0]
            start = feature.start
            end = feature.end
            strand = feature.strand
            
            # Collect additional information for the header
            gene_name = feature.attributes.get('gene', [''])[0]
            protein_id = feature.attributes.get('protein_id', [''])[0]
            
            # Extract the sequence from the corresponding FASTA file
            with open(fasta_file, "r") as fasta_handle:
                for record in SeqIO.parse(fasta_handle, "fasta"):
                    if record.id == feature.seqid:
                        gene_seq = record.seq[start-1:end]
                        if strand == '-':
                            gene_seq = gene_seq.reverse_complement()
                        
                        # Create a comprehensive header
                        header = f"{gene_id} | {gene_name} | {protein_id} | {product} | {feature.seqid}:{start}-{end} ({strand})"
                        acetyltransferase_genes.append((header, gene_seq))
                        break

    # Write the extracted genes to a FASTA file
    with open(output_file, "w") as out_handle:
        for header, gene_seq in acetyltransferase_genes:
            SeqIO.write(SeqIO.SeqRecord(gene_seq, id=header, description=""), out_handle, "fasta")

    print(f"Extracted {len(acetyltransferase_genes)} acetyltransferase genes")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python extract_acetyltransferase_genes.py <gff_file> <fasta_file> <keyword> <output_file>")
        sys.exit(1)
    
    gff_file, fasta_file, keyword, output_file = sys.argv[1:]
    extract_acetyltransferase_genes(gff_file, fasta_file, keyword, output_file)