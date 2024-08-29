process PANAROO_NON_MUTATED {
    tag "PANAROO on non-mutated genomes"
    publishDir "${params.outputDir}/panaroo_non_mutated", mode: 'copy'

    input:
    path gff_files

    output:
    path "panaroo_output", emit: panaroo_results
    path "panaroo_output/pan_genome_reference.fa", emit: pan_genome_reference

    script:
    """
    echo "Input GFF files:"
    ls -l ${gff_files}
    panaroo -i ${gff_files.join(' ')} \
            -o panaroo_output \
            --clean-mode ${params.panaroo_clean_mode} \
            --merge_paralogs \
            --threads ${task.cpus}
    
    # Check if the file exists
    if [ ! -f panaroo_output/pan_genome_reference.fa ]; then
        echo "Error: pan_genome_reference.fa not found!"
        exit 1
    fi
    
    # Post-processing to consolidate gene names
    python3 ${projectDir}/bin/consolidate_gene_names.py panaroo_output/pan_genome_reference.fa panaroo_output/consolidated_pan_genome_reference.fa
    
    # Replace the original file with the consolidated one
    mv panaroo_output/consolidated_pan_genome_reference.fa panaroo_output/pan_genome_reference.fa
    """
}