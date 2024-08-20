// modules/extract_unique_genes.nf

process EXTRACT_UNIQUE_GENES {
    publishDir "${params.outputDir}/unique_genes", mode: 'copy'

    input:
    path gene_data_csv
    path pan_genome_reference

    output:
    path "unique_genes.fa", emit: unique_genes
    path "matched_genes.txt", emit: matched_genes
    path "debug_info.txt", emit: debug_info

    script:
    """
    #!/usr/bin/env python3
    import csv
    from Bio import SeqIO

    keyword = "${params.keyword}".lower()

    # Extract gene names for genes with the keyword in the description
    gene_names = set()
    matched_genes = []
    with open("${gene_data_csv}", 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if 'description' in row and keyword in row['description'].lower():
                if 'gene_name' in row:
                    gene_names.add(row['gene_name'])
                    matched_genes.append(f"{row['gene_name']}: {row['description']}")

    # Write matched genes to file
    with open("matched_genes.txt", "w") as f:
        for gene in matched_genes:
            f.write(f"{gene}\\n")

    # Extract sequences from pan_genome_reference.fa
    unique_sequences = []
    found_genes = set()
    for record in SeqIO.parse("${pan_genome_reference}", "fasta"):
        record_ids = record.id.split('~~~')
        for record_id in record_ids:
            base_id = record_id.split('_')[0]  # Remove any suffixes after underscore
            if base_id in gene_names or record_id in gene_names:
                unique_sequences.append(record)
                for id in record_ids:
                    found_genes.add(id)
                    found_genes.add(id.split('_')[0])
                break  # Stop checking other IDs if we found a match

    not_found_genes = gene_names - found_genes

    # Write unique sequences to file
    SeqIO.write(unique_sequences, "unique_genes.fa", "fasta")

    # Write debug info
    with open("debug_info.txt", "w") as f:
        f.write(f"Total matched genes: {len(gene_names)}\\n")
        f.write(f"Genes found in FASTA: {len(found_genes)}\\n")
        f.write(f"Genes not found in FASTA: {len(not_found_genes)}\\n")
        f.write("\\nGenes not found:\\n")
        for gene in sorted(not_found_genes):
            f.write(f"{gene}\\n")
        f.write("\\nAll found genes:\\n")
        for gene in sorted(found_genes):
            f.write(f"{gene}\\n")

    print(f"Found {len(gene_names)} genes matching the keyword '{keyword}'")
    print(f"Extracted {len(unique_sequences)} sequences")
    print(f"Check debug_info.txt for more details")
    """
}