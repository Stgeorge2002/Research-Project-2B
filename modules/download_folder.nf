process DOWNLOAD_FOLDER {
    tag "Downloading ${url}"
    publishDir "${params.outputDir}/downloaded", mode: 'copy'

    input:
    val(url)

    output:
    path "flattened_files/*.fa", emit: downloaded_fasta

    script:
    """
    set -e
    mkdir -p unzipped_files
    mkdir -p flattened_files
    wget -O downloaded.zip "${url}"
    unzip -d unzipped_files downloaded.zip
    find unzipped_files -type f -name '*.fa' -exec cp {} flattened_files/ \\;
    """
}