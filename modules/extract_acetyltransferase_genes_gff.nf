process EXTRACT_ACETYLTRANSFERASE_GENES_GFF {
    tag "Extracting acetyltransferase genes from ${sampleName}"
    publishDir "${params.outputDir}/extracted_acetyltransferase_genes", mode: 'copy'

    input:
    tuple val(sampleName), path(gff_file)
    val secondStageKeyword

    output:
    tuple val(sampleName), path("${sampleName}_acetyltransferase_genes.fasta"), emit: extracted_genes

    script:
    """
    #!/usr/bin/env python3
    from Bio import SeqIO
    from BCBio import GFF
    import re

    def extract_sequence(parent_seq, start, end, strand):
        seq = parent_seq[start-1:end]
        return seq if strand == 1 else seq.reverse_complement()

    acetyltransferase_genes = []
    with open("${gff_file}", "r") as gff_handle:
        for rec in GFF.parse(gff_handle):
            for feature in rec.features:
                if feature.type == "CDS":
                    product = feature.qualifiers.get("product", [""])[0]
                    if re.search("${secondStageKeyword}", product, re.IGNORECASE):
                        gene_seq = extract_sequence(rec.seq, feature.location.start.position + 1, feature.location.end.position, feature.strand)
                        gene_id = feature.qualifiers.get("locus_tag", ["Unknown"])[0]
                        acetyltransferase_genes.append((gene_id, gene_seq))

    with open("${sampleName}_acetyltransferase_genes.fasta", "w") as out_handle:
        for gene_id, gene_seq in acetyltransferase_genes:
            SeqIO.write(SeqIO.SeqRecord(gene_seq, id=gene_id, description=""), out_handle, "fasta")

    print(f"Extracted {len(acetyltransferase_genes)} acetyltransferase genes from ${sampleName}")
    """
}