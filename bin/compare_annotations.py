#!/usr/bin/env python3

import argparse
import gffutils
from Bio import SeqIO
from collections import defaultdict

def extract_features(gff_file, fna_file):
    db = gffutils.create_db(gff_file, dbfn=':memory:', force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True)
    features = defaultdict(dict)
    
    for feature in db.all_features(featuretype='CDS'):
        gene_id = feature.attributes.get('locus_tag', ['Unknown'])[0]
        gene_name = feature.attributes.get('gene', [''])[0]
        product = feature.attributes.get('product', [''])[0]
        protein_id = feature.attributes.get('protein_id', [''])[0]
        start = feature.start
        end = feature.end
        strand = feature.strand
        
        # Extract sequence
        with open(fna_file, "r") as fna_handle:
            for record in SeqIO.parse(fna_handle, "fasta"):
                if record.id == feature.seqid:
                    gene_seq = record.seq[start-1:end]
                    if strand == '-':
                        gene_seq = gene_seq.reverse_complement()
                    break
        
        features[gene_id] = {
            "gene": gene_name,
            "product": product,
            "protein_id": protein_id,
            "start": start,
            "end": end,
            "strand": strand,
            "sequence": str(gene_seq)
        }
    
    return features

def compare_annotations(query_features, subject_features, query_fasta, subject_fasta):
    comparisons = []
    
    query_seqs = SeqIO.to_dict(SeqIO.parse(query_fasta, "fasta"))
    subject_seqs = SeqIO.to_dict(SeqIO.parse(subject_fasta, "fasta"))
    
    for q_id, q_feature in query_features.items():
        comparison = {
            "query_id": q_id,
            "query_gene": q_feature['gene'],
            "query_product": q_feature['product'],
            "query_protein_id": q_feature['protein_id'],
            "query_location": f"{q_feature['start']}-{q_feature['end']} ({q_feature['strand']})",
            "query_length": len(q_feature['sequence']),
        }
        
        # Find corresponding subject feature
        s_feature = None
        for s_id, feature in subject_features.items():
            if feature['gene'] == q_feature['gene'] or feature['product'] == q_feature['product']:
                s_feature = feature
                comparison["subject_id"] = s_id
                break
        
        if s_feature:
            comparison.update({
                "subject_gene": s_feature['gene'],
                "subject_product": s_feature['product'],
                "subject_protein_id": s_feature['protein_id'],
                "subject_location": f"{s_feature['start']}-{s_feature['end']} ({s_feature['strand']})",
                "subject_length": len(s_feature['sequence']),
            })
            
            # Compare sequences from BLAST results
            query_blast_seq = str(query_seqs[q_id].seq)
            subject_blast_seq = str(subject_seqs[s_id].seq)
            seq_identity = sum(q == s for q, s in zip(query_blast_seq, subject_blast_seq)) / max(len(query_blast_seq), len(subject_blast_seq)) * 100
            comparison["sequence_identity"] = f"{seq_identity:.2f}%"
        else:
            comparison["subject_match"] = "No corresponding feature found"
        
        comparisons.append(comparison)
    
    return comparisons

def write_comparison(comparisons, output_file):
    with open(output_file, "w") as outfile:
        for comp in comparisons:
            outfile.write(f"Query: {comp['query_id']}\n")
            outfile.write(f"  Gene: {comp['query_gene']}\n")
            outfile.write(f"  Product: {comp['query_product']}\n")
            outfile.write(f"  Protein ID: {comp['query_protein_id']}\n")
            outfile.write(f"  Location: {comp['query_location']}\n")
            outfile.write(f"  Sequence length: {comp['query_length']}\n")
            
            if "subject_id" in comp:
                outfile.write(f"Subject: {comp['subject_id']}\n")
                outfile.write(f"  Gene: {comp['subject_gene']}\n")
                outfile.write(f"  Product: {comp['subject_product']}\n")
                outfile.write(f"  Protein ID: {comp['subject_protein_id']}\n")
                outfile.write(f"  Location: {comp['subject_location']}\n")
                outfile.write(f"  Sequence length: {comp['subject_length']}\n")
                outfile.write(f"  Sequence identity: {comp['sequence_identity']}\n")
            else:
                outfile.write(f"{comp['subject_match']}\n")
            
            outfile.write("\n")

def main():
    parser = argparse.ArgumentParser(description="Compare annotations between query and subject")
    parser.add_argument("--query_gff", required=True, help="Query GFF file")
    parser.add_argument("--query_fna", required=True, help="Query FNA file")
    parser.add_argument("--subject_gff", required=True, help="Subject GFF file")
    parser.add_argument("--subject_fna", required=True, help="Subject FNA file")
    parser.add_argument("--query_fasta", required=True, help="Query FASTA file from BLAST hits")
    parser.add_argument("--subject_fasta", required=True, help="Subject FASTA file from BLAST hits")
    parser.add_argument("--output", required=True, help="Output comparison file")
    args = parser.parse_args()

    query_features = extract_features(args.query_gff, args.query_fna)
    subject_features = extract_features(args.subject_gff, args.subject_fna)
    
    comparisons = compare_annotations(query_features, subject_features, args.query_fasta, args.subject_fasta)
    write_comparison(comparisons, args.output)

if __name__ == "__main__":
    main()