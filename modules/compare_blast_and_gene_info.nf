// modules/compare_blast_and_gene_info.nf

process COMPARE_BLAST_AND_GENE_INFO {
    tag "Comparing BLAST results with gene info and reference alignments"
    publishDir "${params.outputDir}/blast_comparison", mode: 'copy'

    input:
    path blast_results
    path gene_info
    path reference_alignments
    path gene_data_csv

    output:
    path "${blast_results.baseName}_with_gene_info_and_reference.txt", emit: comparison_results

    script:
    """
    #!/usr/bin/env python3
    import csv
    from collections import OrderedDict
    import re

    # Define the truncation buffer
    TRUNCATION_BUFFER = 30

    # Read gene_data.csv and create annotation_id to gff_file mapping
    annotation_to_gff = {}
    with open("${gene_data_csv}", 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            annotation_id = row['annotation_id']
            numeric_id = re.search(r'\\d+', annotation_id).group()  # Extract numeric part
            gff_file = row['gff_file']
            annotation_to_gff[numeric_id] = gff_file

    # Read gene info
    gene_info_dict = {}
    with open("${gene_info}", 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) == 2:
                gene_id = parts[0].split()[0]
                gene_info_dict[gene_id] = parts

    # Process BLAST results and compare with gene info
    results = OrderedDict()
    truncated_sequences = {}
    with open("${blast_results}", 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            query_id = row[0]
            if query_id not in results:
                results[query_id] = row
                query_length = int(row[3])
                subject_length = int(row[4])
                length_difference = abs(query_length - subject_length)
                if length_difference > TRUNCATION_BUFFER:
                    truncated_sequences[query_id] = length_difference

    # Process reference alignments
    ref_alignments = {}
    ref_truncated_sequences = {}
    with open("${reference_alignments}", 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            query_id = row[0]
            ref_alignments[query_id] = row
            query_end = int(row[7])
            subject_end = int(row[9])
            length_difference = abs(query_end - subject_end)
            if length_difference > TRUNCATION_BUFFER:
                ref_truncated_sequences[query_id] = length_difference

    # Write comparison results
    with open("${blast_results.baseName}_with_gene_info_and_reference.txt", 'w') as f:
        # Write truncation summary for custom BLAST DB
        f.write(f"Number of truncated sequences in custom BLAST DB: {len(truncated_sequences)}\\n")
        f.write(f"Truncated sequence IDs in custom BLAST DB (ID: truncation amount, GFF file):\\n")
        for seq_id, truncation in truncated_sequences.items():
            numeric_id = re.search(r'\\d+', seq_id).group()  # Extract numeric part
            gff_file = annotation_to_gff.get(numeric_id, "GFF file not found")
            f.write(f"{seq_id}: {truncation}, {gff_file}\\n")
        f.write("\\n" + "="*50 + "\\n\\n")

        # Write truncation summary for reference alignments
        f.write(f"Number of truncated sequences in reference alignments: {len(ref_truncated_sequences)}\\n")
        f.write(f"Truncated sequence IDs in reference alignments (ID: truncation amount, GFF file):\\n")
        for seq_id, truncation in ref_truncated_sequences.items():
            numeric_id = re.search(r'\\d+', seq_id).group()  # Extract numeric part
            gff_file = annotation_to_gff.get(numeric_id, "GFF file not found")
            f.write(f"{seq_id}: {truncation}, {gff_file}\\n")
        f.write("\\n" + "="*50 + "\\n\\n")

        for query_id, entry in results.items():
            # Parse annotation ID from query_id
            match = re.search(r'(_[A-Z]+_\\d+)', query_id)
            if match:
                annotation_id = match.group(1)[1:]  # Remove the leading underscore
                if annotation_id in gene_info_dict:
                    gff_file = gene_info_dict[annotation_id][1]
                    f.write(f"Query ID: {query_id} (GFF file: {gff_file})\\n")
                else:
                    f.write(f"Query ID: {query_id} (GFF file not found)\\n")
            else:
                f.write(f"Query ID: {query_id}\\n")

            f.write(f"Subject ID: {entry[1]}\\n")
            f.write(f"Percentage Identity: {entry[2]}\\n")
            f.write(f"Alignment Length: {entry[3]}\\n")
            f.write(f"Query Length: {entry[3]}\\n")
            f.write(f"Subject Length: {entry[4]}\\n")
            f.write(f"Mismatches: {entry[5]}\\n")
            f.write(f"Gap Opens: {entry[6]}\\n")
            f.write(f"Query Start: {entry[7]}\\n")
            f.write(f"Query End: {entry[8]}\\n")
            f.write(f"Subject Start: {entry[9]}\\n")
            f.write(f"Subject End: {entry[10]}\\n")
            f.write(f"E-value: {entry[11]}\\n")
            f.write(f"Bit Score: {entry[12]}\\n")

            subject_id = entry[1].split()[0]
            if subject_id in gene_info_dict:
                f.write("Subject sequence header:\\n")
                f.write(f"{gene_info_dict[subject_id][0]}\\n")
                f.write("Subject sequence from gene info:\\n")
                f.write(f"{gene_info_dict[subject_id][1]}\\n\\n")
            else:
                f.write(f"Warning: No gene info found for subject ID {subject_id}\\n\\n")

            # Add reference alignment data
            if query_id in ref_alignments:
                ref_entry = ref_alignments[query_id]
                f.write("Reference alignment data:\\n")
                f.write(f"Query ID: {ref_entry[0]}\\n")
                f.write(f"Subject ID: {ref_entry[1]}\\n")
                f.write(f"Percentage Identity: {ref_entry[2]}\\n")
                f.write(f"Alignment Length: {ref_entry[3]}\\n")
                f.write(f"Mismatches: {ref_entry[4]}\\n")
                f.write(f"Gap Opens: {ref_entry[5]}\\n")
                f.write(f"Query Start: {ref_entry[6]}\\n")
                f.write(f"Query End: {ref_entry[7]}\\n")
                f.write(f"Subject Start: {ref_entry[8]}\\n")
                f.write(f"Subject End: {ref_entry[9]}\\n")
                f.write(f"E-value: {ref_entry[10]}\\n")
                f.write(f"Bit Score: {ref_entry[11]}\\n")
            else:
                f.write(f"Warning: No reference alignment found for query ID {query_id}\\n\\n")

            f.write("-"*50 + "\\n\\n")

        # Write warnings for missing gene info
        with open("${blast_results.baseName}_with_gene_info_and_reference.txt", 'a') as f:
            for query_id, entry in results.items():
                subject_id = entry[1].split()[0]
                if subject_id not in gene_info_dict:
                    f.write(f"Warning: No gene info found for subject ID {subject_id}\\n")
    """
}