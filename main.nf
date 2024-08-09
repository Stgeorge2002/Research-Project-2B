#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import modules
include { DOWNLOAD_FOLDER } from './modules/download_folder'
include { PROKKA_INITIAL } from './modules/prokka_initial'
include { PROKKA_MUTATED } from './modules/prokka_mutated'
include { MUTATE_GENES } from './modules/mutate_genes'
include { PANAROO_NON_MUTATED } from './modules/panaroo_non_mutated'
include { PANAROO_MUTATED } from './modules/panaroo_mutated'
include { DWGSIM_MUTATED } from './modules/dwgsim_mutated'
include { DWGSIM_NON_MUTATED } from './modules/dwgsim_non_mutated'
include { EXTRACT_NON_MUTATED_GENOME_KEYWORD_GENES } from './modules/extract_non_mutated_genome_keyword_genes'
include { EXTRACT_MUTATED_GENOME_KEYWORD_GENES } from './modules/extract_mutated_genome_keyword_genes'
include { SNIPPY_MUTATED } from './modules/snippy_mutated'
include { SNIPPY_NON_MUTATED } from './modules/snippy_non_mutated'
include { ANALYZE_NON_MUTATED_GENES } from './modules/analyze_non_mutated_genes'
include { ANALYZE_MUTATED_GENES } from './modules/analyze_mutated_genes'
include { BLAST_ACETYLTRANSFERASES_MUTATED } from './modules/blast_acetyltransferases_mutated'
include { BLAST_ACETYLTRANSFERASES_NON_MUTATED } from './modules/blast_acetyltransferases_non_mutated'
include { CREATE_CUSTOM_ACETYLTRANSFERASE_DB } from './modules/create_custom_acetyltransferase_db'
include { PROCESS_LOCAL_FILES } from './modules/process_local_files'
include { EXTRACT_ALL_KEYWORD_GENES_NON_MUTATED } from './modules/extract_all_keyword_genes_non_mutated'
include { ALIGN_KEYWORD_GENES_TO_REFERENCE } from './modules/align_keyword_genes_to_reference'
include { EXTRACT_ALL_KEYWORD_GENES_MUTATED } from './modules/extract_all_keyword_genes_mutated'
include { ALIGN_KEYWORD_GENES_TO_REFERENCE_MUTATED } from './modules/align_keyword_genes_to_reference_mutated'
include { SNIPPY_REAL_READS } from './modules/snippy_real_reads'
include { EXTRACT_ACETYLTRANSFERASE_GENES_GFF as EXTRACT_ACETYLTRANSFERASE_GENES_GFF_REAL } from './modules/extract_acetyltransferase_genes_gff'
include { EXTRACT_ACETYLTRANSFERASE_GENES_GFF as EXTRACT_ACETYLTRANSFERASE_GENES_GFF_ASSEMBLED } from './modules/extract_acetyltransferase_genes_gff'
include { PROKKA_ASSEMBLED as PROKKA_ASSEMBLED_REAL } from './modules/prokka_assembled'
include { PROKKA_ASSEMBLED as PROKKA_ASSEMBLED_PREASSEMBLED } from './modules/prokka_assembled'
include { ANALYZE_PREMATURE_STOPS_ASSEMBLED as ANALYZE_PREMATURE_STOPS_REAL } from './modules/analyze_premature_stops_assembled'
include { ANALYZE_PREMATURE_STOPS_ASSEMBLED as ANALYZE_PREMATURE_STOPS_PREASSEMBLED } from './modules/analyze_premature_stops_assembled'
include { BLAST_ACETYLTRANSFERASES as BLAST_ACETYLTRANSFERASES_REAL } from './modules/blast_acetyltransferases'
include { BLAST_ACETYLTRANSFERASES as BLAST_ACETYLTRANSFERASES_PREASSEMBLED } from './modules/blast_acetyltransferases'
include { SPADES_ASSEMBLY } from './modules/spades_assembly'
include { SSPACE_SCAFFOLDING } from './modules/sspace_scaffolding'
include { GAPFILLER } from './modules/gapfiller'

// Log info
log.info """
         P A N G E N O M E   A N A L Y S I S   P I P E L I N E    
         ===================================================
         Input file         : ${params.inputFile}
         Output directory   : ${params.outputDir}
         Non-mutated dir    : ${params.nonMutatedDir}
         Mutated dir        : ${params.mutatedDir}
         
         Mutation parameters:
         Num genes to mutate: ${params.num_genes_to_mutate}
         Mutation keyword   : ${params.mutation_keyword}
         
         Gene extraction:
         First stage keyword : ${params.keyword}
         Second stage keyword: ${params.secondStageKeyword}
         
         DWGSIM parameters:
         Error rate         : ${params.dwgsimErrorRate}
         Outer distance     : ${params.dwgsimOuterDistance}
         Standard deviation : ${params.dwgsimStdDev}
         Number of reads    : ${params.dwgsimNumReads}
         Read length        : ${params.dwgsimReadLength}
         
         Resource parameters:
         Max CPUs           : ${params.maxCpus}
         Max memory         : ${params.maxMemory}
         Max time           : ${params.maxTime}
         
         Prokka parameters:
         Genus              : ${params.prokka_genus}
         Species            : ${params.prokka_species}
         Strain             : ${params.prokka_strain}
         
         Panaroo parameters:
         Clean mode         : ${params.panaroo_clean_mode}
         
         BLAST parameters:
         BLAST database     : ${params.blast_db}
         Species taxid      : ${params.blast_species_taxid}
         Word size          : ${params.blast_word_size}
         E-value            : ${params.blast_evalue}
         Percent identity   : ${params.blast_perc_identity}
         Max target seqs    : ${params.blast_max_target_seqs}
         Dust               : ${params.blast_dust}
         Soft masking       : ${params.blast_soft_masking}
         X-drop ungapped    : ${params.blast_xdrop_ungap}
         X-drop gap         : ${params.blast_xdrop_gap}
         Gap open           : ${params.blast_gapopen}
         Gap extend         : ${params.blast_gapextend}
"""

// Parameter validation
if (!params.inputFile || !params.outputDir) {
    error "Input file and output directory parameters are required: --inputFile and --outputDir"
}

// Main workflow
workflow {

    if (params.useLocalFiles) {
        // Process local files
        PROCESS_LOCAL_FILES(params.localFilePath)
        fasta_files = PROCESS_LOCAL_FILES.out.processed_fasta.flatten().map { file ->
            def sampleName = file.baseName
            return tuple(sampleName, file)
        }
    } else {
        // Download and process files from URLs
        input_channel = Channel.fromPath(params.inputFile).splitText().map{it.trim()}
        DOWNLOAD_FOLDER(input_channel)
        fasta_files = DOWNLOAD_FOLDER.out.downloaded_fasta.flatten().map { file -> 
            def sampleName = file.baseName
            return tuple(sampleName, file)
        }
    }
    
    PROKKA_INITIAL(fasta_files)
    PROKKA_INITIAL.out.prokka_results.view { "PROKKA output: $it" }

    all_genomes = PROKKA_INITIAL.out.prokka_results
    selected_genome = all_genomes.randomSample(1)

    // Debug output
    all_genomes.view { "Genome for selection: $it" }
    selected_genome.view { "Selected genome: $it" }
    
    // Mutate the selected genome
    MUTATE_GENES(selected_genome)

    // Debug output for MUTATE_GENES
    MUTATE_GENES.out.mutated_genome.view { "Mutated genome: $it" }
    MUTATE_GENES.out.combined_genes.view { "Combined genes: $it" }
    MUTATE_GENES.out.mutation_info.view { "Mutation info: $it" }

    // Run Prokka on the "mutated" genome
    PROKKA_MUTATED(MUTATE_GENES.out.mutated_genome.map { it -> [it[0], it[2]] })

    // Debug output for DWGSIM input
    PROKKA_MUTATED.out.fna_files.view { "Debug - FNA file for DWGSIM (mutated): $it" }

    // Run DWGSIM on the mutated genome
    DWGSIM_MUTATED(PROKKA_MUTATED.out.fna_files.map { file -> 
        def sampleName = file.getName().toString().tokenize('.')[0]
        return tuple(sampleName, file)
    })

    // Run DWGSIM on the non-mutated genome
    DWGSIM_NON_MUTATED(
        selected_genome
            .map { sampleName, gff, fna -> 
                return tuple(sampleName, file(fna))
            }
    )

    // Prepare inputs for Panaroo
    non_mutated_gffs = PROKKA_INITIAL.out.gff_files
    mutated_gff = PROKKA_MUTATED.out.prokka_mutated_results.map { it[1] }  // Extract only the GFF file

    // Run Panaroo for non-mutated genomes
    PANAROO_NON_MUTATED(non_mutated_gffs.collect())

    // Prepare input for mutated Panaroo
    mutated_base_names = mutated_gff.map { file -> 
        def base_name = file.name.replaceFirst(/_mutated\.gff$/, '')
        return tuple(base_name, file)
    }

    all_gffs_for_mutated = non_mutated_gffs
        .map { file -> 
            def base_name = file.name.replaceFirst(/\.gff$/, '')
            return tuple(base_name, file)
        }
        .join(mutated_base_names, remainder: true)
        .map { base_name, non_mutated, mutated -> 
            mutated ?: non_mutated
        }

    PANAROO_MUTATED(all_gffs_for_mutated.collect())

    EXTRACT_NON_MUTATED_GENOME_KEYWORD_GENES(
        PANAROO_NON_MUTATED.out.panaroo_results
            .combine(selected_genome)
            .map { panaroo_result, selectedSampleName, selectedGff, selectedFna -> 
                tuple(selectedSampleName, file("${panaroo_result}/gene_data.csv"))
            },
        params.keyword
    )

    EXTRACT_MUTATED_GENOME_KEYWORD_GENES(
        PANAROO_MUTATED.out.panaroo_results
            .combine(selected_genome)
            .map { panaroo_result, selectedSampleName, selectedGff, selectedFna -> 
                tuple(selectedSampleName, file("${panaroo_result}/gene_data.csv"))
            },
        params.keyword
    )
    
    CREATE_CUSTOM_ACETYLTRANSFERASE_DB()

    // BLAST non-mutated acetyltransferases
    BLAST_ACETYLTRANSFERASES_NON_MUTATED(
        EXTRACT_NON_MUTATED_GENOME_KEYWORD_GENES.out.extracted_genes,
        CREATE_CUSTOM_ACETYLTRANSFERASE_DB.out.db.collect()
    )

    // BLAST mutated acetyltransferases
    BLAST_ACETYLTRANSFERASES_MUTATED(
        EXTRACT_MUTATED_GENOME_KEYWORD_GENES.out.extracted_genes,
        CREATE_CUSTOM_ACETYLTRANSFERASE_DB.out.db.collect()
    )

    // Extract all genes with keyword from non-mutated Panaroo results
    EXTRACT_ALL_KEYWORD_GENES_NON_MUTATED(
        PANAROO_NON_MUTATED.out.panaroo_results.map { it -> file("${it}/gene_data.csv") },
        params.keyword
    )

    // Align extracted acetyltransferases to pangenome reference
    ALIGN_KEYWORD_GENES_TO_REFERENCE(
        EXTRACT_ALL_KEYWORD_GENES_NON_MUTATED.out.extracted_genes,
        PANAROO_NON_MUTATED.out.panaroo_results.map { it -> file("${it}/pan_genome_reference.fa") }
    )

    // Extract all genes with keyword from mutated Panaroo results
    EXTRACT_ALL_KEYWORD_GENES_MUTATED(
        PANAROO_MUTATED.out.panaroo_results.map { it -> file("${it}/gene_data.csv") },
        params.keyword
    )

    // Align extracted mutated acetyltransferases to pangenome reference
    ALIGN_KEYWORD_GENES_TO_REFERENCE_MUTATED(
        EXTRACT_ALL_KEYWORD_GENES_MUTATED.out.extracted_genes_mutated,
        PANAROO_MUTATED.out.panaroo_results.map { it -> file("${it}/pan_genome_reference.fa") }
    )

    // Analyze non-mutated genes
    ANALYZE_NON_MUTATED_GENES(EXTRACT_NON_MUTATED_GENOME_KEYWORD_GENES.out.extracted_genes)

    // Analyze mutated genes
    ANALYZE_MUTATED_GENES(EXTRACT_MUTATED_GENOME_KEYWORD_GENES.out.extracted_genes)

    // Process analysis results
    ANALYZE_NON_MUTATED_GENES.out.analysis_results
        .mix(ANALYZE_MUTATED_GENES.out.analysis_results)
        .collectFile(name: "${params.outputDir}/combined_analysis_results.json", newLine: true, sort: true)
        .view { "Combined analysis results: $it" }

    // Update the DWGSIM_MUTATED and DWGSIM_NON_MUTATED output channel definitions
    DWGSIM_MUTATED.out.simulated_reads
        .map { sampleName, read1, read2 -> tuple(sampleName, file(read1), file(read2)) }
        .set { dwgsim_reads_mutated }

    DWGSIM_NON_MUTATED.out.simulated_reads
        .map { sampleName, read1, read2 -> tuple(sampleName, file(read1), file(read2)) }
        .set { dwgsim_non_mutated_reads }

    // Note: Unzipping of DWGSIM output is not necessary as downstream processes can handle gzipped files
    // If unzipping is required, consider adding a separate process for this step

    // Run SNIPPY for mutated comparison
    SNIPPY_MUTATED(
        dwgsim_reads_mutated,
        EXTRACT_NON_MUTATED_GENOME_KEYWORD_GENES.out.extracted_genes
    )

    // Run SNIPPY for non-mutated comparison
    SNIPPY_NON_MUTATED(
        dwgsim_non_mutated_reads,
        EXTRACT_MUTATED_GENOME_KEYWORD_GENES.out.extracted_genes
    )

    // Process or publish Snippy results
    SNIPPY_MUTATED.out.snippy_results.view { "Snippy mutated output: $it" }
    SNIPPY_NON_MUTATED.out.snippy_results.view { "Snippy non-mutated output: $it" }

    // Conditional execution of real reads analysis
    if (params.runRealReadsAnalysis) {
        // Prepare input for real reads analysis
        if (params.useRealReads && params.readsDir) {
            real_reads_ch = Channel.fromFilePairs("${params.readsDir}/*_{1,2}.fastq.gz")
            reference_genome = PANAROO_NON_MUTATED.out.pangenome_fasta
        } else {
            real_reads_ch = DWGSIM_MUTATED.out.simulated_reads
            reference_genome = EXTRACT_NON_MUTATED_GENOME_KEYWORD_GENES.out.extracted_genes
            log.info "Using DWGSIM mutated reads as no real reads were provided."
        }

        // Run Snippy on real or simulated reads
        SNIPPY_REAL_READS(real_reads_ch, reference_genome)

        // Run SPAdes on Snippy reads
        SPADES_ASSEMBLY(SNIPPY_REAL_READS.out.snippy_results)

        // Prepare input for SSPACE_SCAFFOLDING
        sspace_input = SPADES_ASSEMBLY.out.assembled_genome
            .join(SNIPPY_REAL_READS.out.snippy_results)
            .map { sampleName, assembly, snippy_dir ->
                def reads = file("${snippy_dir}/*.fastq.gz")
                tuple(sampleName, file(assembly), reads[0], reads[1])
            }

        // Add this before calling SSPACE_SCAFFOLDING
        sspace_input.view { "SSPACE input: $it" }

        // Run SSPACE_SCAFFOLDING
        SSPACE_SCAFFOLDING(sspace_input)

        // Add this after calling SSPACE_SCAFFOLDING
        SSPACE_SCAFFOLDING.out.scaffolds.view { "SSPACE output: $it" }

        // Prepare input for GAPFILLER
        gapfiller_input = SSPACE_SCAFFOLDING.out.scaffolds
            .join(SNIPPY_REAL_READS.out.snippy_results)
            .map { sampleName, scaffolds, snippy_dir ->
                def reads = file("${snippy_dir}/*.fastq.gz")
                tuple(sampleName, file(scaffolds), reads[0], reads[1])
            }

        // Run GAPFILLER
        GAPFILLER(gapfiller_input)

        // Run Prokka on gap-filled scaffolds
        PROKKA_ASSEMBLED_REAL(GAPFILLER.out.gap_filled_scaffolds)

        // Extract acetyltransferase genes from Prokka GFF files using the second stage keyword
        EXTRACT_ACETYLTRANSFERASE_GENES_GFF_REAL(PROKKA_ASSEMBLED_REAL.out.gff_file, params.secondStageKeyword)

        // Analyze premature stops in extracted acetyltransferase genes
        ANALYZE_PREMATURE_STOPS_REAL(EXTRACT_ACETYLTRANSFERASE_GENES_GFF_REAL.out.extracted_genes)

        // BLAST extracted acetyltransferase genes against the custom database
        BLAST_ACETYLTRANSFERASES_REAL(EXTRACT_ACETYLTRANSFERASE_GENES_GFF_REAL.out.extracted_genes, CREATE_CUSTOM_ACETYLTRANSFERASE_DB.out.db.collect())

        // You can add more analysis steps here, comparing results from both stages
    } else if (params.assembledGenomesDir) {
        // Use pre-assembled genomes if provided
        assembled_genomes_ch = Channel.fromPath("${params.assembledGenomesDir}/*.fa")
            .map { file -> tuple(file.baseName, file) }
        
        // Run Prokka on pre-assembled genomes
        PROKKA_ASSEMBLED_PREASSEMBLED(assembled_genomes_ch)

        // Extract acetyltransferase genes from Prokka GFF files using the second stage keyword
        EXTRACT_ACETYLTRANSFERASE_GENES_GFF_ASSEMBLED(PROKKA_ASSEMBLED_PREASSEMBLED.out.gff_file, params.secondStageKeyword)

        // Analyze premature stops in extracted acetyltransferase genes
        ANALYZE_PREMATURE_STOPS_PREASSEMBLED(EXTRACT_ACETYLTRANSFERASE_GENES_GFF_ASSEMBLED.out.extracted_genes)

        // BLAST extracted acetyltransferase genes against the custom database
        BLAST_ACETYLTRANSFERASES_PREASSEMBLED(EXTRACT_ACETYLTRANSFERASE_GENES_GFF_ASSEMBLED.out.extracted_genes, CREATE_CUSTOM_ACETYLTRANSFERASE_DB.out.db.collect())
    }
}

// Workflow completion handler
workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Execution duration: $workflow.duration"
}

// Workflow error handler
workflow.onError {
    log.error "Pipeline execution stopped with the following message: ${workflow.errorMessage}"
}