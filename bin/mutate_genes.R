#!/usr/bin/env Rscript

# Load required libraries
library(GenomicFeatures)
library(Biostrings)
library(rtracklayer)
library(parallel)

args <- commandArgs(trailingOnly = TRUE)
gff_file <- args[1]
fasta_file <- args[2]
output_file <- args[3]
combined_genes_file <- args[4]
mutation_info_file <- args[5]
num_genes_to_mutate <- as.numeric(args[6])
keyword <- args[7]
mode <- args[8]  # Add mode parameter


# Parse the GFF file to get acetyltransferase genes
parse_gff_by_keyword <- function(gff_file, keyword) {
  gff <- import(gff_file, format = "gff3")
  fields <- c("Name", "gene", "product")
  matching_records <- NULL
  for (field in fields) {
    if (field %in% colnames(mcols(gff))) {
      matching_records <- gff[grepl(keyword, mcols(gff)[[field]], ignore.case = TRUE)]
      if (length(matching_records) > 0) {
        cat("Found matching records in field:", field, "\n")
        break
      }
    }
  }
  if (is.null(matching_records) || length(matching_records) == 0) {
    cat("No matching records found in GFF file for keyword:", keyword, "\n")
  }
  return(matching_records)
}

# Extract gene sequences from FASTA based on GFF records
extract_gene_sequences <- function(gff_records, fasta_file) {
  fasta <- readDNAStringSet(fasta_file)
  gene_sequences <- list()
  for (i in seq_along(gff_records)) {
    record <- gff_records[i]
    seqid <- as.character(seqnames(record))
    gene_name <- as.character(mcols(record)$Name)
    strand <- as.character(strand(record))
    if (seqid %in% names(fasta)) {
      start <- start(record)
      end <- end(record)
      seq <- subseq(fasta[[seqid]], start=start, end=end)
      gene_sequences[[paste0(gene_name, ":", seqid, ":", start, "-", end)]] <- list(seq=seq, start=start, end=end, strand=strand, gene_name=gene_name, seqid=seqid)
    } else {
      cat("Sequence ID", seqid, "not found in FASTA file.\n")
    }
  }
  return(gene_sequences)
}

# Function to mutate the sequence
mutate_sequence <- function(seq_info, mode = "mutate") {
  seq <- seq_info$seq
  strand <- seq_info$strand
  
  mutable_seq <- as.character(seq)
  seq_length <- nchar(mutable_seq)
  
  # Split the sequence into codons
  codons <- strsplit(mutable_seq, "(?<=...)", perl=TRUE)[[1]]
  
  if (mode == "mutate") {
    # Choose a random codon to mutate (excluding the first and last codons)
    codon_index <- sample(2:(length(codons)-1), 1)
    
    # Determine the appropriate stop codon based on the strand
    stop_codon <- if(strand == "+") "TAA" else "TTA"
    
    # Replace the chosen codon with the stop codon
    original_codon <- codons[codon_index]
    codons[codon_index] <- stop_codon
    
    # Reconstruct the mutated sequence
    mutated_seq <- paste(codons, collapse="")
    
    # Calculate the position of the mutation in the original sequence
    position <- (codon_index - 1) * 3 + 1
    
    mutation_type <- "mutation"
  } else {
    stop("Invalid mode. Use 'mutate'.")
  }
  
  cat("Mutating codon at index:", codon_index, "in sequence of length:", seq_length, "on strand:", strand, "\n")
  cat("Mutation type:", mutation_type, "\n")
  cat("Introduced stop codon:", stop_codon, "\n")
  
  mutation_info <- list(
    seq = DNAString(mutated_seq), 
    position = position,
    original = original_codon, 
    strand = strand,
    type = mutation_type,
    introduced_codon = stop_codon
  )
  
  return(mutation_info)
}

# Print lengths of the contigs
print_contig_lengths <- function(fasta_sequences) {
  cat("Contig lengths:\n")
  for (name in names(fasta_sequences)) {
    cat(name, ":", nchar(fasta_sequences[[name]]), "bp\n")
  }
}

# Verify the mutation
verify_mutation <- function(original_seq, mutated_seq, mutation_info, seq_info) {
  strand <- seq_info$strand
  expected <- mutation_info$introduced_codon
  
  original_seq <- as.character(original_seq)
  mutated_seq <- as.character(mutated_seq)
  position <- mutation_info$position
  
  found <- substr(mutated_seq, position, position + 2)
  
  original_codon <- mutation_info$original
  
  if (found != expected) {
    cat("Mutation verification failed at position", position, "\n")
    cat("Gene:", seq_info$gene_name, "Contig:", seq_info$seqid, "Strand:", strand, "\n")
    cat("Expected", expected, "found", found, "\n")
    cat("Original codon:", original_codon, "\n")
    cat("Original sequence context:", substr(original_seq, max(1, position - 6), min(nchar(original_seq), position + 8)), "\n")
    cat("Mutated sequence context:", substr(mutated_seq, max(1, position - 6), min(nchar(mutated_seq), position + 8)), "\n")
    return(FALSE)
  } else {
    cat("Mutation verification succeeded at position", position, "\n")
    cat("Gene:", seq_info$gene_name, "Contig:", seq_info$seqid, "Strand:", strand, "\n")
    cat("Expected", expected, "found", found, "\n")
    cat("Original codon:", original_codon, "\n")
    cat("Original sequence context:", substr(original_seq, max(1, position - 6), min(nchar(original_seq), position + 8)), "\n")
    cat("Mutated sequence context:", substr(mutated_seq, max(1, position - 6), min(nchar(mutated_seq), position + 8)), "\n")
  }
  return(TRUE)
}

# Function to wrap sequences to a specified width
wrap_sequence <- function(sequence, width = 60) {
  sapply(seq(1, nchar(sequence), by = width), function(x) substr(sequence, x, min(nchar(sequence), x + width - 1)))
}

# Write sequences to a FASTA file with specified width
write_fasta_with_wrap <- function(fasta, filepath, width = 60) {
  con <- file(filepath, "w")
  for (i in seq_along(fasta)) {
    # Clean the name to remove any leading/trailing whitespace and any random spaces
    name <- gsub("^\\s+|\\s+$", "", names(fasta)[i])
    cat(">", name, "\n", file = con)
    
    # Clean the sequence to remove any leading/trailing whitespace
    sequence <- gsub("^\\s+|\\s+$", "", as.character(fasta[[i]]))
    wrapped_seq <- wrap_sequence(sequence, width)
    writeLines(wrapped_seq, con)
  }
  close(con)
}

# Write combined original and mutated gene sequences to a file
write_combined_genes <- function(original_genes, mutated_genes, combined_file) {
  con <- file(combined_file, "w")
  for (name in names(original_genes)) {
    # Write original gene
    cat(">", name, "_original\n", file = con)
    wrapped_seq <- wrap_sequence(as.character(original_genes[[name]]), 60)
    writeLines(wrapped_seq, con)
    
    # Write mutated gene
    cat(">", name, "_mutated\n", file = con)
    wrapped_seq <- wrap_sequence(as.character(mutated_genes[[name]]), 60)
    writeLines(wrapped_seq, con)
  }
  close(con)
}

# Write mutation information to a file
write_mutation_info <- function(mutation_info, filepath) {
  con <- file(filepath, "w")
  for (info in mutation_info) {
    writeLines(paste("Gene:", info$gene_name), con)
    writeLines(paste("Contig:", info$seqid), con)
    writeLines(paste("Position:", info$position), con)
    writeLines(paste("Strand:", info$strand), con)
    writeLines(paste("Original Codon:", info$original_codon), con)
    writeLines(paste("Mutated Codon: TAA"), con)
    writeLines("\n", con)
  }
  close(con)
}

# Function for global verification of all mutations
verify_all_mutations <- function(original_fasta, mutated_sequences, gene_sequences) {
  for (name in names(mutated_sequences)) {
    seq_info <- gene_sequences[[name]]
    seqid <- seq_info$seqid
    start <- seq_info$start
    end <- seq_info$end
    strand <- seq_info$strand
    mutated_seq <- as.character(mutated_sequences[[name]]$seq)
    current_seq <- as.character(subseq(original_fasta[[seqid]], start=start, end=end))
    
    if (strand == "-") {
      current_seq <- as.character(reverseComplement(DNAString(current_seq)))
    }
    
    if (current_seq != mutated_seq) {
      cat("Global verification failed for", name, "\n")
      cat("Strand:", strand, "\n")
      cat("Expected:", mutated_seq, "\n")
      cat("Found   :", current_seq, "\n")
      cat("Mismatch positions:\n")
      for (i in 1:nchar(mutated_seq)) {
        if (substr(mutated_seq, i, i) != substr(current_seq, i, i)) {
          cat("Position", i, ": Expected", substr(mutated_seq, i, i), "Found", substr(current_seq, i, i), "\n")
        }
      }
      stop(paste("Global verification failed for", name))
    }
  }
  cat("Global verification passed for all mutated sequences\n")
}

# Main function to execute the workflow
main <- function(gff_file, fasta_file, output_file, combined_genes_file, mutation_info_file, num_genes_to_mutate, keyword, mode = "mutate") {
  gff_records <- parse_gff_by_keyword(gff_file, keyword)
  
  if (is.null(gff_records) || length(gff_records) == 0) {
    cat("No matching records found in GFF file for keyword:", keyword, "\n")
    return()
  }
  
  cat("Number of matching GFF records:", length(gff_records), "\n")
  
  gene_sequences <- extract_gene_sequences(gff_records, fasta_file)
  
  if (length(gene_sequences) == 0) {
    cat("No gene sequences extracted from FASTA file.\n")
    return()
  }
  
  cat("Number of gene sequences extracted:", length(gene_sequences), "\n")
  
  # Limit the number of genes to mutate
  genes_to_mutate <- sample(names(gene_sequences), min(num_genes_to_mutate, length(gene_sequences)))
  
  # Read the original FASTA file
  original_fasta <- readDNAStringSet(fasta_file)
  
  cat("Original contig lengths:\n")
  print_contig_lengths(original_fasta)
  
  # Introduce one mutation per selected gene sequence
  mutated_sequences <- lapply(gene_sequences[genes_to_mutate], function(seq) mutate_sequence(seq, mode = mode))
  
  # Containers for original and mutated full gene sequences and mutation information
  original_genes <- list()
  mutated_genes <- list()
  mutation_info_list <- list()
  
  # Replace the original sequences with the mutated sequences
  for (name in names(mutated_sequences)) {
    seq_info <- gene_sequences[[name]]
    seqid <- seq_info$seqid
    start <- seq_info$start
    end <- seq_info$end
    strand <- seq_info$strand
    mutated_seq_info <- mutated_sequences[[name]]
    mutated_seq <- mutated_seq_info$seq
    mutation_position <- mutated_seq_info$position
    original_seq <- as.character(subseq(original_fasta[[seqid]], start=start, end=end))
    
    # Length check before replacement
    if (nchar(as.character(mutated_seq)) != (end - start + 1)) {
      stop(paste("Length mismatch for", name, ": original length =", (end - start + 1), "mutated length =", nchar(as.character(mutated_seq))))
    }
    
    # Replace original sequence with the mutated sequence
    cat("Replacing sequence in contig:", seqid, "from position", start, "to", end, "\n")
    if (strand == "+") {
      original_fasta[[seqid]] <- replaceAt(original_fasta[[seqid]], IRanges(start=start, end=end), mutated_seq)
    } else {
      rev_comp_mutated <- reverseComplement(mutated_seq)
      original_fasta[[seqid]] <- replaceAt(original_fasta[[seqid]], IRanges(start=start, end=end), rev_comp_mutated)
    }
    
    # Post-replacement verification
    post_replacement_seq <- as.character(subseq(original_fasta[[seqid]], start=start, end=end))
    if (strand == "+") {
      if (post_replacement_seq != as.character(mutated_seq)) {
        stop(paste("Post-replacement verification failed for", name))
      }
    } else {
      if (post_replacement_seq != as.character(rev_comp_mutated)) {
        stop(paste("Post-replacement verification failed for", name))
      }
    }
    
    # Verification step
    cat("Verifying mutations for", name, "...\n")
    if (!verify_mutation(original_seq, as.character(mutated_seq_info$seq), mutated_seq_info, seq_info)) {
      cat("Verification failed for", name, "\n")
    }
    
    # Store original and mutated gene sequences
    original_genes[[name]] <- DNAString(original_seq)
    mutated_genes[[name]] <- DNAString(as.character(mutated_seq))
    
    # Store mutation information
    mutation_info_list[[name]] <- list(
      gene_name = seq_info$gene_name,
      seqid = seq_info$seqid,
      position = mutation_position,
      original_codon = mutated_seq_info$original,
      strand = seq_info$strand,
      type = mutated_seq_info$type
    )
  }
  
  # Write to a temporary file
  temp_output_file <- tempfile(pattern = "temp_mutated_", fileext = ".fna")
  write_fasta_with_wrap(original_fasta, filepath = temp_output_file)
  
  # Wait for the temporary file to be fully written
  while(file.size(temp_output_file) == 0) {
    Sys.sleep(0.1)  # Sleep for 100ms
  }
  
  # Ensure the file content is stable (not still being written)
  last_size <- -1
  current_size <- file.size(temp_output_file)
  while(last_size != current_size) {
    Sys.sleep(0.5)  # Sleep for 500ms
    last_size <- current_size
    current_size <- file.size(temp_output_file)
  }
  
  # Copy the temporary file to the final output file
  file.copy(temp_output_file, output_file, overwrite = TRUE)
  
  # Remove the temporary file
  file.remove(temp_output_file)
  
  cat("Saved mutated sequences to", output_file, "\n")
  
  # Global verification step
  verify_all_mutations(original_fasta, mutated_sequences, gene_sequences)
  
  # Read the mutated FASTA file and print lengths
  mutated_fasta <- readDNAStringSet(output_file)
  cat("Mutated contig lengths:\n")
  print_contig_lengths(mutated_fasta)
  
  # Write combined original and mutated gene sequences to a single file
  write_combined_genes(original_genes, mutated_genes, combined_genes_file)
  
  # Write mutation information to a file
  write_mutation_info(mutation_info_list, mutation_info_file)
  
  # Create the sentinel file
  sentinel_file <- paste0(output_file, ".done")
  file.create(sentinel_file)
  cat("Created sentinel file:", sentinel_file, "\n")
  
  # Close all open connections
  closeAllConnections()
}

# Execute the main function
main(gff_file, fasta_file, output_file, combined_genes_file, mutation_info_file, num_genes_to_mutate, keyword, mode)