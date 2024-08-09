Please see below a diagram to explain the workflow process.


This diagram omits the extraction of keyword genes, and also the python script to scan for premature stop codons. 


![image](https://github.com/user-attachments/assets/9c8c4d6e-3fef-47cb-8e32-c6aa7b170006)




nextflow run main.nf \
  --search_term 'acetyltransferase' \
  --inputFile "$projectDir/input_urls.txt" \
  --outputDir "results" \
  --num_genes_to_mutate 10 \
  --mutation_mode 'mutate' \
  --dwgsimErrorRate 0.0001 \
  --dwgsimOuterDistance 500 \
  --dwgsimStdDev 50 \
  --dwgsimNumReads 1000000 \
  --dwgsimReadLength 150 \
  --snippyCpus 4 \
  --prokka_genus 'Streptococcus' \
  --prokka_species 'pneumoniae' \
  --prokka_strain 'SPNEU' \
  --panaroo_clean_mode 'strict' \
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
  --useLocalFiles false \
  --localFilePath "$projectDir/local_genomes" \
  --readsDir "${launchDir}/real_reads" \
  --runRealReadsAnalysis false \
  --useRealReads false \
  -with-docker
