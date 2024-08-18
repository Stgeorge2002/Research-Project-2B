process GATHER_TOP_BLAST_INFO_AND_ANNOTATIONS {
    tag "Gathering top BLAST info and annotations for ${sampleName}"
    publishDir "${params.outputDir}/top_blast_info_and_annotations", mode: 'copy'

    input:
    tuple val(sampleName), path(blast_results_with_seq)
    tuple val(sampleName), path(extracted_genes)
    path pangenome_alignment_results
    path snippy_analysis
    path query_gff

    output:
    tuple val(sampleName), path("${sampleName}_top_blast_info_and_annotations.tsv"), emit: top_blast_info_and_annotations

    script:
    """
    #!/usr/bin/env python3
    import csv
    import re
    from Bio import SeqIO, Entrez
    from BCBio import GFF

    # Set your email for NCBI queries
    Entrez.email = "${params.email}"

    def extract_genbank_id(blast_hit):
        match = re.search(r'cds_(\\w+\\.\\d+)', blast_hit)
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

    print("Parsing extracted genes...")
    extracted_genes = parse_extracted_genes("${extracted_genes}")

    print("Parsing query GFF file...")
    query_annotations = parse_query_gff("${query_gff}")

    print("Parsing pangenome alignment results...")
    pangenome_alignments = parse_pangenome_alignment("${pangenome_alignment_results}")

    print("Parsing Snippy analysis results...")
    snippy_data = parse_snippy_analysis("${snippy_analysis}")

    print("Processing BLAST results...")
    top_hits = []
    with open("${blast_results_with_seq}", 'r') as f:
        reader = csv.reader(f, delimiter='\\t')
        for row in reader:
            qseqid = row[0]
            if not any(hit[0] == qseqid for hit in top_hits):
                top_hits.append(row)

    print("Writing output file...")
    with open("${sampleName}_top_blast_info_and_annotations.tsv", 'w') as out_f:
        for hit in top_hits:
            qseqid, sseqid, ident, align_len, qlen, slen, mismatches, gap_opens, qstart, qend, sstart, send, evalue, bitscore, qseq, sseq = hit[:16]
            
            query_gene_name, query_product = query_annotations.get(qseqid, ("Unknown", "No annotation found"))
            
            out_f.write("BLAST Custom DB Results:\\n")
            out_f.write(f"Query ID: {qseqid}\\n")
            out_f.write(f"Query Gene Name: {query_gene_name}\\n")
            out_f.write(f"Query Product: {query_product}\\n")
            out_f.write(f"Subject ID: {sseqid}\\n")
            out_f.write(f"Identity: {ident}%\\n")
            out_f.write(f"Alignment Length: {align_len}\\n")
            out_f.write(f"Mismatches: {mismatches}\\n")
            out_f.write(f"Gap Opens: {gap_opens}\\n")
            out_f.write(f"Query Start-End: {qstart}-{qend}\\n")
            out_f.write(f"Subject Start-End: {sstart}-{send}\\n")
            out_f.write(f"E-value: {evalue}\\n")
            out_f.write(f"Bit Score: {bitscore}\\n")
            out_f.write(f"Query Length: {qlen}\\n")
            out_f.write(f"Subject Length: {slen}\\n")
            
            # Fetch GenBank data for the subject sequence
            subject_genbank_id = extract_genbank_id(sseqid)
            if subject_genbank_id:
                try:
                    subject_record = fetch_genbank_data(subject_genbank_id)
                    subject_gene_name = get_gene_name(subject_record)
                    subject_description = subject_record.description
                    out_f.write(f"Subject GenBank Gene Name: {subject_gene_name}\\n")
                    out_f.write(f"Subject GenBank Description: {subject_description}\\n")
                except Exception as e:
                    out_f.write(f"Error fetching Subject GenBank data: {str(e)}\\n")
            
            out_f.write("\\n")  # Add a gap
            
            # Add pangenome alignment information
            out_f.write("BLAST Pangenome Reference Results:\\n")
            pangenome_subject_id = None
            if qseqid in pangenome_alignments:
                pan_align = pangenome_alignments[qseqid]
                pangenome_subject_id = pan_align['pangenome_subject_id']
                out_f.write(f"Pangenome Subject ID: {pan_align['pangenome_subject_id']}\\n")
                out_f.write(f"Pangenome Identity: {pan_align['pangenome_identity']}%\\n")
                out_f.write(f"Pangenome Alignment Length: {pan_align['pangenome_alignment_length']}\\n")
                out_f.write(f"Pangenome Query Length: {pan_align['pangenome_qlen']}\\n")
                out_f.write(f"Pangenome Subject Length: {pan_align['pangenome_send']}\\n")  # Use 'pangenome_send' for subject length
                out_f.write(f"Pangenome Mismatches: {pan_align['pangenome_mismatches']}\\n")
                out_f.write(f"Pangenome Gaps: {pan_align['pangenome_gaps']}\\n")
                out_f.write(f"Pangenome Query Start-End: {pan_align['pangenome_qstart']}-{pan_align['pangenome_qend']}\\n")
                out_f.write(f"Pangenome Subject Start-End: {pan_align['pangenome_sstart']}-{pan_align['pangenome_send']}\\n")
                out_f.write(f"Pangenome E-value: {pan_align['pangenome_evalue']}\\n")
                out_f.write(f"Pangenome Bit Score: {pan_align['pangenome_bitscore']}\\n")
            else:
                out_f.write("No pangenome alignment found for this query\\n")
            
            out_f.write("\\n")  # Add a gap
            
            # Add Snippy analysis information
            out_f.write("Snippy Analysis Results:\\n")
            if pangenome_subject_id and pangenome_subject_id in snippy_data:
                for line in snippy_data[pangenome_subject_id]:
                    if not line.startswith("Gene Name:") and not line.startswith("Product:"):
                        out_f.write(f"{line}\\n")
            elif qseqid in snippy_data:
                for line in snippy_data[qseqid]:
                    if not line.startswith("Gene Name:") and not line.startswith("Product:"):
                        out_f.write(f"{line}\\n")
            else:
                out_f.write("No Snippy analysis found for this query\\n")
            
            out_f.write("\\n" + "-" * 60 + "\\n\\n")

    print(f"Done! Output written to ${sampleName}_top_blast_info_and_annotations.tsv")
    """
}