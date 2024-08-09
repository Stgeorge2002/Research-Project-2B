process GAPFILLER {
    tag "GapFiller on ${sampleName}"
    container 'tabath123/sspace-gapfiller:latest'
    
    input:
    tuple val(sampleName), path(scaffolds), path(read1), path(read2)

    output:
    tuple val(sampleName), path("${sampleName}_gap_filled_scaffolds.fasta"), emit: gap_filled_scaffolds
    
    script:
    """
    echo "lib1 ${read1} ${read2} 500 0.75 FR" > library.txt
    
    GapFiller.pl -l library.txt -s ${scaffolds} -b ${sampleName}_gapfiller
    
    mv ${sampleName}_gapfiller/standard_output.gapfilled.final.fa ${sampleName}_gap_filled_scaffolds.fasta
    """
}