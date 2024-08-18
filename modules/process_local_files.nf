process PROCESS_LOCAL_FILES {
    publishDir "${params.outputDir}/processed_local_files", mode: 'copy'

    input:
    path local_files

    output:
    path "processed_files/*.fa", emit: processed_fasta

    script:
    """
    mkdir -p processed_files
    echo "Contents of input files:"
    ls -l ${local_files}
    
    for file in ${local_files}; do
        new_name=\$(basename "\$file" | sed 's/\\.contigs_velvet\\.fa\$/.fa/')
        cp "\$file" "processed_files/\$new_name"
    done
    
    echo "Processed files:"
    ls -l processed_files/
    
    if [ ! -f processed_files/*.fa ]; then
        echo "No .fa files found. Creating a dummy file."
        touch processed_files/dummy.fa
    fi
    """
}