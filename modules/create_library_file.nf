process CREATE_LIBRARY_FILE {
    tag "Creating library file for ${sampleName}"
    
    input:
    tuple val(sampleName), path(read1), path(read2)

    output:
    tuple val(sampleName), path('libraries.txt'), emit: library_file

    script:
    """
    echo "lib1 ${read1} ${read2} 500 0.25 FR" > libraries.txt
    """
}