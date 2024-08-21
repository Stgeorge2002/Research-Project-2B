# PanReadPipe

## Diagram of the whole pipeline:

![image](https://github.com/user-attachments/assets/1f246edd-45b2-410e-a6c2-4bad5ebb0c97)

# To run the pipeline:

PanReadPipe is a two-stage process, note all custom Docker Images are pushed to my Dockerhub account 'tabath123'

## First stage (PanG_DB):

note, the / indicates 'either'

Use 'nextflow run PanG_DB.nf -profile PanG_DB --profile PanGDB --useLocalFiles true/false --outputdir >your run name<' for the main process.

### Diagram of PanG_DB
![image](https://github.com/user-attachments/assets/e456019f-51cb-48dd-a387-53d23de4a6b1)

This run command has parameters onlyCore = true/false and exTrun = true/false, of which both are set default true, as this bypasses the positive control generation of mutated genomes and also allows for the automatic grabbing of the detected genomes with truncated genes of interest. Here the parameter --search_term >your search-term< to target genes of interest. 

After running this, it is here that the genenomes selected in the /">your run name<_further_analysis" and the output data in /">your run name<_PanG_DB_Output", (see /example_outputs/Mutated_Output_PanG_DB_Example.txt), and review the genomes for selection. /example/Mutated_Output_PanG_DB_Example). 

Here, I would also recommend removing the input files (if done locally) that have the same naming as the selected trucated detected, and then re run this first step -profile PanG_DB --profile PanGDB --useLocalFiles true --exTrun false --outputdir >different your run name<. (see /example_outputs/Non_Mutated_Output_PanG_DB_Example.txt) as this will give better results for the BLAST pangenome alignment in the second step, but will not affect the BLAST NCBI reference alignment or SNIPPY results/detection. 

## Second stage (XtraDetect):

Now for the second step Use 'nextflow run XtraDetect.nf -profile XtraDetect --reads/fa --outputDir >the same different run name<' 

### Diagram of XtraDetect
![image](https://github.com/user-attachments/assets/a7f84e63-495a-4b3f-bf76-a3ca9d35bb24)


--fa uses all of the genomes located in >your run name<_further_analysis<, ## PLEASE NOTE ##, when using --fa. Reads can also be used if nessecary, and they must be paired and located in the /real_reals directory unless parameterized otherwise. ## PLEASE NOTE ## when using --reads please set --skipUnicycler false (this will be set to default soon)

This will produce a file containing all alignments and trucations reports, aswell as the custom SNIPPY premature stop codon detection (see see /example_outputs/Mutated_Output_XtraDetect_Example.txt .

### For a view of the wider parameters used here please see the nextflow.config file, as these can all be tuned in the script. 
