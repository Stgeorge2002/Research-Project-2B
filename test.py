
#!/usr/bin/env python3
import csv
from collections import OrderedDict, defaultdict
import re
import sys
from Bio import SeqIO

# Hardcoded file paths for testing
blast_results = "results/blast_results/non_mutated/all_non_mutated_extracted_genes_blast_results_with_seq.tsv"
gene_info = "results/Custom_BLAST_DB/gene_info.tsv"
reference_alignments = "results/aligned_acetyltransferases/acetyltransferases_alignments_with_seq.txt"
gene_data_csv = "results/panaroo_non_mutated/panaroo_output/gene_data.csv"
pangenome_reference = "results/panaroo_non_mutated/panaroo_output/pan_genome_reference.fa"
extracted_genes_file = "results/extracted_genes/all_non_mutated_extracted_genes.fasta"

# Define the truncation buffer
TRUNCATION_BUFFER = 30

# Read the extracted genes file
extracted_genes = {}
with open(extracted_genes_file, "r") as f:
    for record in SeqIO.parse(f, "fasta"):
        extracted_genes[record.id] = str(record.seq)

# Read gene_data.csv and create annotation_id to gff_file mapping
annotation_to_gff = {}
with open(gene_data_csv, 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        annotation_id = row['annotation_id']
        gff_file = row['gff_file']
        annotation_to_gff[annotation_id] = gff_file

# Read gene info
gene_info_dict = {}
with open(gene_info, 'r') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) == 2:
            gene_id = parts[0].split()[0]
            gene_info_dict[gene_id] = parts

# Process pangenome reference
pangenome_sequences = {}
for record in SeqIO.parse(pangenome_reference, "fasta"):
    pangenome_sequences[record.id] = str(record.seq)
    # Also add entries for each individual gene name
    for gene in record.id.split('~~~'):
        pangenome_sequences[gene] = str(record.seq)

# Initialize dictionaries for truncated sequences
truncated_sequences = {}
ref_truncated_sequences = {}

# Process BLAST results and compare with gene info
results = OrderedDict()
with open(blast_results, 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    for row in reader:
        query_id = row[0]
        if query_id not in results:
            results[query_id] = row
            query_length = int(row[4])  # Use query length from BLAST results
            subject_length = int(row[5])  # Use subject length from BLAST results
            length_difference = subject_length - query_length
            if abs(length_difference) > TRUNCATION_BUFFER:
                truncated_sequences[query_id] = length_difference
            
            # Apply the same logic to find GFF file for each query
            annotation_id = query_id.split('_', 2)[1] + '_' + query_id.split('_', 2)[2].split('_')[0]
            gff_file = annotation_to_gff.get(annotation_id, "GFF file not found")
            results[query_id].append(gff_file)  # Add GFF file to the results

# Process reference alignments
ref_alignments = {}
with open(reference_alignments, 'r') as f:
    reader = csv.reader(f, delimiter='\t')
    for row in reader:
        query_id = row[0]
        ref_alignments[query_id] = row
        
        # Check for truncation in reference alignments
        query_length = len(extracted_genes.get(query_id, ""))  # Use length of extracted sequence
        subject_id = row[1]
        subject_length = len(pangenome_sequences.get(subject_id, ""))  # Use pangenome reference length
        length_difference = subject_length - query_length
        if abs(length_difference) > TRUNCATION_BUFFER:
            ref_truncated_sequences[query_id] = length_difference

# Write comparison results
output_file = f"{blast_results.split('/')[-1].split('.')[0]}_with_gene_info_and_reference.txt"
with open(output_file, 'w') as f:
    # Write truncation summary for custom BLAST DB
    f.write(f"Number of truncated sequences in custom BLAST DB: {len(truncated_sequences)}\n")
    f.write(f"Truncated sequence IDs in custom BLAST DB (ID: truncation amount, GFF file):\n")
    for seq_id, truncation in truncated_sequences.items():
        annotation_id = seq_id.split('_', 2)[1] + '_' + seq_id.split('_', 2)[2].split('_')[0]
        gff_file = annotation_to_gff.get(annotation_id, "GFF file not found")
        f.write(f"{seq_id}: {truncation}, {gff_file}\n")
    f.write("\n" + "="*50 + "\n\n")

    # Write truncation summary for reference alignments
    f.write(f"Number of truncated sequences in reference alignments: {len(ref_truncated_sequences)}\n")
    f.write(f"Truncated sequence IDs in reference alignments (ID: truncation amount, GFF file):\n")
    for seq_id, truncation in ref_truncated_sequences.items():
        annotation_id = seq_id.split('_', 2)[1] + '_' + seq_id.split('_', 2)[2].split('_')[0]
        gff_file = annotation_to_gff.get(annotation_id, "GFF file not found")
        f.write(f"{seq_id}: {truncation}, {gff_file}\n")
    f.write("\n" + "="*50 + "\n\n")

    for query_id, entry in results.items():
        # Parse annotation ID from query_id
        annotation_id = query_id.split('_', 2)[1] + '_' + query_id.split('_', 2)[2].split('_')[0]
        gff_file = annotation_to_gff.get(annotation_id, "GFF file not found")
        f.write(f"Query ID: {query_id} (GFF file: {gff_file})\n")
        f.write(f"Subject ID: {entry[1]}\n")
        f.write(f"Percentage Identity: {entry[2]}\n")
        f.write(f"Alignment Length: {entry[3]}\n")
        f.write(f"Query Length: {entry[4]}\n")
        f.write(f"Subject Length: {entry[5]}\n")
        f.write(f"Mismatches: {entry[6]}\n")
        f.write(f"Gap Opens: {entry[7]}\n")
        f.write(f"Query Start: {entry[8]}\n")
        f.write(f"Query End: {entry[9]}\n")
        f.write(f"Subject Start: {entry[10]}\n")
        f.write(f"Subject End: {entry[11]}\n")
        f.write(f"E-value: {entry[12]}\n")
        f.write(f"Bit Score: {entry[13]}\n")

        # Add original sequence from extracted genes file
        if query_id in extracted_genes:
            f.write("Original sequence from extracted genes file:\n")
            f.write(f"{extracted_genes[query_id]}\n")
        else:
            f.write("Original sequence not found in extracted genes file\n")

        subject_id = entry[1].split()[0]
        if subject_id in gene_info_dict:
            f.write("Subject sequence header:\n")
            f.write(f"{gene_info_dict[subject_id][0]}\n")
            f.write("Subject sequence from gene info:\n")
            f.write(f"{gene_info_dict[subject_id][1]}\n")

        # Add reference alignment data
        if query_id in ref_alignments:
            ref_entry = ref_alignments[query_id]
            f.write("Reference alignment data:\n")
            f.write(f"Query ID: {ref_entry[0]}\n")
            f.write(f"Subject ID: {ref_entry[1]}\n")
            f.write(f"Percentage Identity: {ref_entry[2]}\n")
            f.write(f"Alignment Length: {ref_entry[3]}\n")
            f.write(f"Mismatches: {ref_entry[4]}\n")
            f.write(f"Gap Opens: {ref_entry[5]}\n")
            f.write(f"Query Start: {ref_entry[6]}\n")
            f.write(f"Query End: {ref_entry[7]}\n")
            f.write(f"Subject Start: {ref_entry[8]}\n")
            f.write(f"Subject End: {ref_entry[9]}\n")
            f.write(f"E-value: {ref_entry[10]}\n")
            f.write(f"Bit Score: {ref_entry[11]}\n")
            f.write(f"Query Sequence Length: {len(extracted_genes.get(query_id, ''))}\n")
            f.write(f"Subject Sequence Length: {len(pangenome_sequences.get(ref_entry[1], ''))}\n")
            
            # Add pangenome reference sequence if available
            subject_id = ref_entry[1]
            pangenome_seq = pangenome_sequences.get(subject_id)
            
            if pangenome_seq:
                f.write("Pangenome reference sequence:\n")
                f.write(f"{pangenome_seq}\n")
                f.write(f"Pangenome reference sequence length: {len(pangenome_seq)}\n")
                length_diff = len(pangenome_seq) - len(extracted_genes.get(query_id, ''))
                f.write(f"Length difference (Pangenome vs Reference Query): {length_diff}\n")
                
                if abs(length_diff) > TRUNCATION_BUFFER:
                    f.write(f"Warning: Possible truncation detected. Length difference: {length_diff}\n")
            else:
                f.write(f"Warning: No pangenome reference sequence found for Subject ID: {subject_id}\n")
        else:
            f.write(f"Warning: No reference alignment found for query ID {query_id}\n\n")

        f.write("-"*50 + "\n\n")

    # Write warnings for missing gene info
    for query_id, entry in results.items():
        subject_id = entry[1].split()[0]
        if subject_id not in gene_info_dict:
            f.write(f"Warning: No gene info found for subject ID {subject_id}\n")

# Create a set of unique truncated GFF file names
truncated_gff_files = set()
for seq_dict in [truncated_sequences, ref_truncated_sequences]:
    for seq_id, length_diff in seq_dict.items():
        annotation_id = seq_id.split('_', 2)[1] + '_' + seq_id.split('_', 2)[2].split('_')[0]
        gff_file = annotation_to_gff.get(annotation_id, "GFF file not found")
        truncated_gff_files.add(gff_file)

# Write truncated GFF file names to a file
with open("truncated_genes.txt", 'w') as f:
    for gff_file in truncated_gff_files:
        f.write(f"{gff_file}\n")

# Debug: Print the contents of the truncated_genes.txt file
print("Contents of truncated_genes.txt:")
with open("truncated_genes.txt", 'r') as f:
    print(f.read())

# New section to print mismatched alignments
print("\nMismatched alignments:")
for query_id, blast_entry in results.items():
    if query_id in ref_alignments:
        ref_entry = ref_alignments[query_id]
        blast_identity = float(blast_entry[2])
        ref_identity = float(ref_entry[2])
        
        if abs(blast_identity - ref_identity) > 0.1:  # 0.1% threshold for mismatch
            print(f"\nQuery ID: {query_id}")
            print(f"BLAST alignment: Identity {blast_identity}%, Length {blast_entry[3]}, Subject {blast_entry[1]}")
            print(f"Reference alignment: Identity {ref_identity}%, Length {ref_entry[3]}, Subject {ref_entry[1]}")
    else:
        print(f"\nQuery ID: {query_id}")
        print(f"BLAST alignment: Identity {blast_entry[2]}%, Length {blast_entry[3]}, Subject {blast_entry[1]}")
        print("No reference alignment found")

print("\nReference alignments without corresponding BLAST results:")
for query_id in ref_alignments:
    if query_id not in results:
        ref_entry = ref_alignments[query_id]
        print(f"\nQuery ID: {query_id}")
        print(f"Reference alignment: Identity {ref_entry[2]}%, Length {ref_entry[3]}, Subject {ref_entry[1]}")
        print("No BLAST alignment found")