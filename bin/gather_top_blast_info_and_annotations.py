#!/usr/bin/env python3
import csv
import re
from Bio import SeqIO, Entrez
from BCBio import GFF
import sys

# Add this line to define TRUNCATION_BUFFER
TRUNCATION_BUFFER = 30

# Get command-line arguments
sampleName = sys.argv[1]
blast_results_with_seq = sys.argv[2]
extracted_genes = sys.argv[3]
pangenome_alignment_results = sys.argv[4]
snippy_analysis = sys.argv[5]
query_gff = sys.argv[6]
email = sys.argv[7]
pangenome_reference = sys.argv[8]  # New argument for pangenome reference file

# Set your email for NCBI queries
Entrez.email = email

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

def parse_extracted_genes(fasta_file):
    extracted_genes = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        gene_id = record.id
        header = record.description
        extracted_genes[gene_id] = header
    return extracted_genes

def parse_query_gff(gff_file):
    query_annotations = {}
    for rec in GFF.parse(gff_file):
        for feature in rec.features:
            if feature.type == "CDS":
                gene_id = feature.id
                gene_name = feature.qualifiers.get("gene", [""])[0]
                product = feature.qualifiers.get("product", ["No product"])[0]
                query_annotations[gene_id] = (gene_name, product)
    return query_annotations

def parse_pangenome_alignment(alignment_file):
    alignments = {}
    with open(alignment_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 12:
                query_id = fields[0]
                alignments[query_id] = {
                    'pangenome_subject_id': fields[1],
                    'pangenome_identity': fields[2],
                    'pangenome_alignment_length': fields[3],
                    'pangenome_qlen': fields[3],  # Query length from the 4th column
                    'pangenome_mismatches': fields[4],
                    'pangenome_gaps': fields[5],
                    'pangenome_qstart': fields[6],
                    'pangenome_qend': fields[7],
                    'pangenome_sstart': fields[8],
                    'pangenome_send': fields[9],
                    'pangenome_slen': fields[9],  # Subject length from the 10th column
                    'pangenome_evalue': fields[10],
                    'pangenome_bitscore': fields[11]
                }
    return alignments

def parse_snippy_analysis(analysis_file):
    snippy_data = {}
    current_gene = None
    with open(analysis_file, 'r') as f:
        for line in f:
            if line.startswith("Gene ID:"):
                current_gene = line.split(":")[1].strip()
                snippy_data[current_gene] = [line.strip()]
            elif current_gene and line.strip():
                snippy_data[current_gene].append(line.strip())
    return snippy_data

# Add this new function to check for truncation
def is_truncated(query_length, subject_length):
    return abs(query_length - subject_length) > TRUNCATION_BUFFER

# New function to parse pangenome reference sequences
def parse_pangenome_reference(fasta_file):
    pangenome_sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        pangenome_sequences[record.id] = str(record.seq)
        # Also add entries for each individual gene name
        for gene in record.id.split('~~~'):
            pangenome_sequences[gene] = str(record.seq)
    return pangenome_sequences

print("Parsing extracted genes...")
extracted_genes = parse_extracted_genes(extracted_genes)

print("Parsing query GFF file...")
query_annotations = parse_query_gff(query_gff)

print("Parsing pangenome alignment results...")
pangenome_alignments = parse_pangenome_alignment(pangenome_alignment_results)

print("Parsing Snippy analysis results...")
snippy_data = parse_snippy_analysis(snippy_analysis)

print("Parsing pangenome reference sequences...")
pangenome_sequences = parse_pangenome_reference(pangenome_reference)

print("Processing BLAST results...")
top_hits = []
truncated_sequences = []
with open(blast_results_with_seq, 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    for row in reader:
        qseqid = row[0]
        if not any(hit[0] == qseqid for hit in top_hits):
            top_hits.append(row)
            query_length = int(row[3])  # 4th column for query length
            subject_length = int(row[5])  # 6th column for subject length
            if is_truncated(query_length, subject_length):
                truncated_sequences.append(qseqid)

print("Writing output file...")
with open(f"{sampleName}_top_blast_info_and_annotations.tsv", 'w') as out_f:
    # Write truncation summary
    out_f.write(f"Number of truncated sequences in custom BLAST DB: {len(truncated_sequences)}\n")
    out_f.write("Truncated sequence IDs in custom BLAST DB:\n")
    for seq_id in truncated_sequences:
        out_f.write(f"{seq_id}\n")
    out_f.write("\n" + "="*50 + "\n\n")

    for hit in top_hits:
        qseqid, sseqid, ident, align_len, qlen, slen, mismatches, gap_opens, qstart, qend, sstart, send, evalue, bitscore, qseq, sseq = hit[:16]
        
        query_gene_name, query_product = query_annotations.get(qseqid, ("Unknown", "No annotation found"))
        
        out_f.write("BLAST Custom DB Results:\n")
        out_f.write(f"Query ID: {qseqid}\n")
        out_f.write(f"Query Gene Name: {query_gene_name}\n")
        out_f.write(f"Query Product: {query_product}\n")
        out_f.write(f"Subject ID: {sseqid}\n")
        out_f.write(f"Identity: {ident}%\n")
        out_f.write(f"Alignment Length: {align_len}\n")
        out_f.write(f"Mismatches: {mismatches}\n")
        out_f.write(f"Gap Opens: {gap_opens}\n")
        out_f.write(f"Query Start-End: {qstart}-{qend}\n")
        out_f.write(f"Subject Start-End: {sstart}-{send}\n")
        out_f.write(f"E-value: {evalue}\n")
        out_f.write(f"Bit Score: {bitscore}\n")
        out_f.write(f"Query Length: {qlen}\n")
        out_f.write(f"Subject Length: {slen}\n")
        
        # Add truncation information
        if qseqid in truncated_sequences:
            out_f.write(f"Truncation: Yes (difference: {abs(int(qlen) - int(slen))})\n")
        else:
            out_f.write("Truncation: No\n")
        
        # Fetch GenBank data for the subject sequence
        subject_genbank_id = extract_genbank_id(sseqid)
        if subject_genbank_id:
            try:
                subject_record = fetch_genbank_data(subject_genbank_id)
                subject_gene_name = get_gene_name(subject_record)
                subject_description = subject_record.description
                out_f.write(f"Subject GenBank Gene Name: {subject_gene_name}\n")
                out_f.write(f"Subject GenBank Description: {subject_description}\n")
            except Exception as e:
                out_f.write(f"Error fetching Subject GenBank data: {str(e)}\n")
        
        out_f.write("\n")  # Add a gap
        
        # Add pangenome alignment information
        out_f.write("BLAST Pangenome Reference Results:\n")
        pangenome_subject_id = None
        if qseqid in pangenome_alignments:
            pan_align = pangenome_alignments[qseqid]
            pangenome_subject_id = pan_align['pangenome_subject_id']
            out_f.write(f"Pangenome Subject ID: {pan_align['pangenome_subject_id']}\n")
            out_f.write(f"Pangenome Identity: {pan_align['pangenome_identity']}%\n")
            out_f.write(f"Pangenome Alignment Length: {pan_align['pangenome_alignment_length']}\n")
            out_f.write(f"Pangenome Query Length: {pan_align['pangenome_qlen']}\n")
            
            # Use actual subject length from pangenome reference
            pangenome_subject_length = len(pangenome_sequences.get(pangenome_subject_id, ""))
            out_f.write(f"Pangenome Subject Length: {pangenome_subject_length}\n")
            
            out_f.write(f"Pangenome Mismatches: {pan_align['pangenome_mismatches']}\n")
            out_f.write(f"Pangenome Gaps: {pan_align['pangenome_gaps']}\n")
            out_f.write(f"Pangenome Query Start-End: {pan_align['pangenome_qstart']}-{pan_align['pangenome_qend']}\n")
            out_f.write(f"Pangenome Subject Start-End: {pan_align['pangenome_sstart']}-{pan_align['pangenome_send']}\n")
            out_f.write(f"Pangenome E-value: {pan_align['pangenome_evalue']}\n")
            out_f.write(f"Pangenome Bit Score: {pan_align['pangenome_bitscore']}\n")
            
            # Add truncation information for pangenome alignment
            pan_query_length = int(pan_align['pangenome_qlen'])
            if is_truncated(pan_query_length, pangenome_subject_length):
                out_f.write(f"Pangenome Truncation: Yes (difference: {abs(pan_query_length - pangenome_subject_length)})\n")
            else:
                out_f.write("Pangenome Truncation: No\n")
        else:
            out_f.write("No pangenome alignment found for this query\n")
        
        out_f.write("\n")  # Add a gap
        
        # Add Snippy analysis information
        out_f.write("Snippy Analysis Results:\n")
        if pangenome_subject_id and pangenome_subject_id in snippy_data:
            snippy_info = snippy_data[pangenome_subject_id]
        elif qseqid in snippy_data:
            snippy_info = snippy_data[qseqid]
        else:
            snippy_info = None

        if snippy_info:
            for line in snippy_info:
                if not line.startswith("Gene Name:") and not line.startswith("Product:"):
                    out_f.write(f"{line}\n")
            
            # Check for premature stop codons
            premature_stops = [pos for pos in snippy_info if "Premature stop codons:" in pos]
            if premature_stops:
                out_f.write(f"{premature_stops[0]}\n")
            else:
                out_f.write("Premature stop codons: 0\n")
        else:
            out_f.write("No Snippy analysis found for this query\n")
        
        out_f.write("\n" + "-" * 60 + "\n\n")

print(f"Done! Output written to {sampleName}_top_blast_info_and_annotations.tsv")

# New section to print unmatched alignments and analyses
print("\nUnmatched alignments and analyses:")

# Set of all query IDs
all_query_ids = set(hit[0] for hit in top_hits)

# Unmatched pangenome alignments
unmatched_pangenome = set(pangenome_alignments.keys()) - all_query_ids
if unmatched_pangenome:
    print("\nUnmatched pangenome alignments:")
    for query_id in unmatched_pangenome:
        print(f"  {query_id}")

# Unmatched Snippy analyses
unmatched_snippy = set(snippy_data.keys()) - all_query_ids - set(pangenome_alignments.keys())
if unmatched_snippy:
    print("\nUnmatched Snippy analyses:")
    for query_id in unmatched_snippy:
        print(f"  {query_id}")

# BLAST alignments without pangenome or Snippy matches
blast_without_matches = all_query_ids - set(pangenome_alignments.keys()) - set(snippy_data.keys())
if blast_without_matches:
    print("\nBLAST alignments without pangenome or Snippy matches:")
    for query_id in blast_without_matches:
        print(f"  {query_id}")

print("\nAnalysis complete.")