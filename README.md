Please see below a diagram to explain the workflow process. tbr ~ to-be-removed 

![image](https://github.com/user-attachments/assets/9c8c4d6e-3fef-47cb-8e32-c6aa7b170006)

This diagram shows only the processes that we discussed previously, and they all work. Below is a outline of the parameters assessible via command line. 

IMPORTANT:Note ~ To run this, use the command "nextflow run main.nf --runRealReadsAnalysis false -with-docker" to run the pipeline as discussed last meeting. This is done as the new read reads part (with spades, sspace and gapfiller) isnt working yet. 


Commands:
# To run the pipeline
nextflow run main.nf \
# Parameters:
  # Search term used for all extraction, alignment and mutation processes, to focus on genes of interest, default to 'acetyltransferase.
  --search_term 'acetyltransferase' \
  --inputFile "$projectDir/input_urls.txt" \
  --outputDir "results" \

  # Mutation parameters
  --num_genes_to_mutate 10 \
  --mutation_mode 'mutate' \

  # DWGSIM parameters
  --dwgsimErrorRate 0.0001 \
  --dwgsimOuterDistance 500 \
  --dwgsimStdDev 50 \
  --dwgsimNumReads 1000000 \
  --dwgsimReadLength 150 \

  # Snippy parameters (tbr)
  --snippyCpus 4 \
 
  # Prokka parameters
  --prokka_genus 'Streptococcus' \
  --prokka_species 'pneumoniae' \
  --prokka_strain 'SPNEU' \

  # Panaroo parameters
  --panaroo_clean_mode 'strict' \

  # Blast parameters
  --blast_max_target_seqs 500 \
  --blast_word_size 7 \
  --blast_evalue "1e-10" \
  --blast_perc_identity 90 \
  --blast_dust "no" \
  --blast_soft_masking false \
  --blast_xdrop_ungap 0 \
  --blast_xdrop_gap 0 \
  --blast_gapopen 0 \
  --blast_gapextend 0 \

  # Enable true for local file use, defaults to file download from input_urls.txt
  --useLocalFiles false \

  # Creates a path for local files
  --localFilePath "$projectDir/local_genomes" \

  # Creates a path for reads. 
  --readsDir "${launchDir}/real_reads" \

  # Set to true for reads analysis enabling. 
  --runRealReadsAnalysis false \

  # Set to true for 
  --useRealReads false \
  -with-docker
