process SSPACE_SCAFFOLDING {
    tag "SSPACE on ${sampleName}"
    container 'tabath123/sspace-gapfiller:latest'
    
    input:
    tuple val(sampleName), path(assembly), path(read1), path(read2)

    output:
    tuple val(sampleName), path("${sampleName}_scaffolds.fasta"), emit: scaffolds
    
    script:
    """
    echo "lib1 ${read1} ${read2} 500 0.75 FR" > library.txt
    
    SSPACE_Standard_v3.0.pl -l library.txt -s ${assembly} -x 0 -k 5 -b ${sampleName}_scaffolding
    
    mv ${sampleName}_scaffolding/standard_output.final.scaffolds.fasta ${sampleName}_scaffolds.fasta
    """
}