import re
import argparse
from Bio import Entrez, SeqIO

# Set your email for NCBI queries
Entrez.email = "your_email@example.com"

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

def reorganize_results(combined_results):
    reorganized = {}
    for gene_id, data in combined_results.items():
        fasta_desc = data['FASTA']
        gene_name = fasta_desc.split('|')[1].strip() if '|' in fasta_desc else gene_id
        
        if gene_name not in reorganized:
            reorganized[gene_name] = {'FASTA': [], 'Snippy': [], 'BLAST': []}
        
        reorganized[gene_name]['FASTA'].append((gene_id, fasta_desc))
        if data['Snippy']:
            reorganized[gene_name]['Snippy'].append((gene_id, data['Snippy']))
        if data['BLAST']:
            reorganized[gene_name]['BLAST'].append((gene_id, data['BLAST']))
    
    return reorganized

def resolve_gene_names(reorganized):
    resolved = {}
    for gene_name, data in reorganized.items():
        base_name = re.sub(r'_\d+$', '', gene_name)
        if base_name not in resolved:
            resolved[base_name] = {'FASTA': [], 'Snippy': [], 'BLAST': []}
        
        resolved[base_name]['FASTA'].extend(data['FASTA'])
        resolved[base_name]['Snippy'].extend(data['Snippy'])
        resolved[base_name]['BLAST'].extend(data['BLAST'])
    
    return resolved

def main(snippy_file, blast_file, fasta_file, output_file):
    # Read FASTA file for gene descriptions
    fasta_data = {}
    with open(fasta_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                fields = line.strip().split()
                gene_id = fields[0][1:]  # Remove the '>' character
                description = " ".join(fields[1:])
                fasta_data[gene_id] = description

    # Read Snippy analysis
    snippy_data = {}
    with open(snippy_file, "r") as f:
        current_gene = None
        for line in f:
            line = line.strip()
            if line.startswith("Gene ID:"):
                current_gene = line.split(":")[1].strip()
                snippy_data[current_gene] = {}
            elif current_gene and ":" in line:
                key, value = line.split(":", 1)
                snippy_data[current_gene][key.strip()] = value.strip()

    # Read BLAST analysis
    blast_data = {}
    with open(blast_file, "r") as f:
        current_gene = None
        for line in f:
            line = line.strip()
            if line.startswith("Reference:"):
                current_gene = line.split(":")[1].strip()
                blast_data[current_gene] = {}
            elif current_gene and ":" in line:
                key, value = line.split(":", 1)
                blast_data[current_gene][key.strip()] = value.strip()

    # Combine results
    combined_results = {}
    for gene in set(list(snippy_data.keys()) + list(blast_data.keys()) + list(fasta_data.keys())):
        combined_results[gene] = {
            "FASTA": fasta_data.get(gene, "Unknown"),
            "Snippy": snippy_data.get(gene, {}),
            "BLAST": blast_data.get(gene, {})
        }

    # Reorganize results
    reorganized_results = reorganize_results(combined_results)

    # Resolve gene names
    resolved_results = resolve_gene_names(reorganized_results)

    # Output results to a file
    with open(output_file, "w") as f:
        f.write("Combined Results:\n\n")
        for gene_name, data in resolved_results.items():
            f.write(f"Gene ID: {gene_name}\n")
            
            for gene_id, fasta_desc in data['FASTA']:
                f.write(f"FASTA Description: {fasta_desc}\n")
            
            for gene_id, snippy in data['Snippy']:
                f.write("Snippy Analysis:\n")
                for key, value in snippy.items():
                    if key not in ["Gene Name", "Product"]:
                        f.write(f"  {key}: {value}\n")
            
            for gene_id, blast in data['BLAST']:
                f.write("BLAST Analysis:\n")
                for key, value in blast.items():
                    if key not in ["Reference Gene", "Reference Product"]:
                        f.write(f"  {key}: {value}\n")
                
                # Fetch GenBank data for the best database match
                best_match = blast.get('Best Database Match', '')
                genbank_id = extract_genbank_id(best_match)
                if genbank_id:
                    try:
                        record = fetch_genbank_data(genbank_id)
                        db_gene_name = get_gene_name(record)
                        db_description = record.description
                        f.write(f"  Database Gene Name: {db_gene_name}\n")
                        f.write(f"  Database Description: {db_description}\n")
                    except Exception as e:
                        f.write(f"  Error fetching GenBank data: {str(e)}\n")
            
            f.write("-" * 50 + "\n\n")

    print(f"Results have been written to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine Snippy, BLAST, and FASTA results")
    parser.add_argument("--snippy", required=True, help="Path to Snippy analysis file")
    parser.add_argument("--blast", required=True, help="Path to BLAST analysis file")
    parser.add_argument("--fasta", required=True, help="Path to extracted genes FASTA file")
    parser.add_argument("--output", required=True, help="Path to output file")
    args = parser.parse_args()

    main(args.snippy, args.blast, args.fasta, args.output)