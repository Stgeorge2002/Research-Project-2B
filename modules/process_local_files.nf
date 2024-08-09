process PROCESS_LOCAL_FILES {
    tag "Processing local files from ${localPath}"
    publishDir "${params.outputDir}/processed_local", mode: 'copy'

    input:
    path localPath

    output:
    path "flattened_files/*.fa", emit: processed_fasta

    script:
    """
    mkdir -p flattened_files
    find "${localPath}" -type f \\( -name "*.fa" -o -name "*.fasta" -o -name "*.fna" \\) -exec cp {} flattened_files/ \\;
    """
}