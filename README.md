# PanReadPipe.

## Workflow process.

![image](https://github.com/user-attachments/assets/8fb04729-f62d-404e-be60-49cea307ff8d)

# Commands:
### To run the pipeline

All custom Docker Images are posted to my Dockerhub account 'tabath123'

Use 'nextflow run main.nf' for the main process.

Use 'nextflow run main.nf --runRealReadsAnalysis true --useRealReads True --readsDir >Your_directory<' to run additional new raw reads, or with '--useRealReads false' to use preloaded DWGSIM reads.  

Use 'nextflow run main.nf --runRealReadsAnalysis true --assembledGenomeDir >Your_directory<' to run an additional new 'real' assembled FASTA file.  
\
\
## Here is a list of other parameters that can be specified at the command line:

  ### Search term used for all processes, default to 'acetyltransferase:
  --search_term 'acetyltransferase' 

  ## Sets the main input .FA files and result output.  
  --inputFile "$projectDir/input_urls.txt" \
  --outputDir "results" 

  ### Mutation parameters:
  --num_genes_to_mutate 10 \
  --mutation_mode 'mutate' 

  ### DWGSIM parameters:
  --dwgsimErrorRate 0.0001 \
  --dwgsimOuterDistance 500 \
  --dwgsimStdDev 50 \
  --dwgsimNumReads 1000000 \ 
  --dwgsimReadLength 150 

  ### Prokka parameters:
  --prokka_genus 'Streptococcus' \
  --prokka_species 'pneumoniae' \
  --prokka_strain 'SPNEU' \

  ### Panaroo parameters:
  --panaroo_clean_mode 'strict' 

  ### Blast parameters:
  --blast_max_target_seqs 500 \
  --blast_word_size 7 \
  --blast_evalue "1e-10" \
  --blast_perc_identity 90 \
  --blast_dust "no" \
  --blast_soft_masking false \
  --blast_xdrop_ungap 0 \
  --blast_xdrop_gap 0 \
  --blast_gapopen 0 \
  --blast_gapextend 0 
  
  ### Local file processing:
  --useLocalFiles = false \
  --localFilePath = "$projectDir/local_genomes"

  ### New parameters for real reads analysis:
  --readsDir = "${launchDir}/real_reads" \
  --runRealReadsAnalysis = false \
  --useRealReads = false \
  --assembledGenomesDir = null \
  --genus = 'Streptococcus' \
  --species = 'pneumoniae' \
  --strain = 'SPNEU'

  ### New parameter for acetyltransferase search:
  --acetyltransferase_search_term = "Streptococcus pneumoniae[Organism] AND ${params.search_term}[Protein] AND gene[Feature Key]"

  ### Unicycler parameters:
  --unicycler_mode = 'bold' \
  --unicycler_min_fasta_length = 300 \
  --unicycler_kmer_count = 4 \
  --unicycler_min_component_size = 300 \
  --unicycler_min_dead_end_size = 300 \
  --unicycler_min_bridge_qual = 3 \
  --unicycler_keep_temp = false \
  --unicycler_no_correct = false \
  --unicycler_verbosity = 2 

  ### FASTA preprocessing parameters:
  --min_length = 200 \
  --remove_duplicates = true

  ### FASTP parameters:
  --qualified_quality_phred = 15 \
  --unqualified_percent_limit = 40 \
  --cut_mean_quality = 20 \
  --cut_window_size = 4 \
  --cut_front = true \
  --cut_tail = true \
  --cut_front_window_size = 1 \
  --cut_front_mean_quality = 20 \
  --cut_tail_window_size = 1 \
  --cut_tail_mean_quality = 20 \
  --length_required = 50
  

