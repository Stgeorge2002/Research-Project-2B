#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import modules
include { DOWNLOAD_FOLDER } from './modules/download_folder'
include { FASTA_PREPROCESS as FASTA_PREPROCESS_INITIAL } from './modules/fasta_preprocess'
include { PROKKA as PROKKA_INITIAL } from './modules/prokka'
include { PANAROO_NON_MUTATED } from './modules/panaroo_non_mutated'
include { ALIGN_KEYWORD_GENES_TO_REFERENCE } from './modules/align_keyword_genes_to_reference'
include { CREATE_CUSTOM_ACETYLTRANSFERASE_DB } from './modules/create_custom_acetyltransferase_db'
include { BLAST_ACETYLTRANSFERASES_NON_MUTATED as BLAST_ALL_KEYWORD_GENES_NON_MUTATED } from './modules/blast_acetyltransferases_non_mutated'
include { ANALYZE_NON_MUTATED_GENES } from './modules/analyze_non_mutated_genes'
include { EXTRACT_ALL_KEYWORD_GENES_NON_MUTATED } from './modules/extract_all_keyword_genes_non_mutated'
include { MUTATE_GENES } from './modules/mutate_genes'
include { PROKKA as PROKKA_MUTATED } from './modules/prokka'
include { DWGSIM_MUTATED } from './modules/dwgsim_mutated'
include { DWGSIM_NON_MUTATED } from './modules/dwgsim_non_mutated'
include { PROCESS_LOCAL_FILES } from './modules/process_local_files'
include { COMPARE_BLAST_AND_GENE_INFO as COMPARE_BLAST_AND_GENE_INFO_NON_MUTATED } from './modules/compare_blast_and_gene_info'
include { SEARCH_TRUNCATED_GFF } from './modules/search_truncated_gff'


// Main workflow
workflow {
    // Core processes
    if (params.useLocalFiles) {
        Channel.fromPath("${params.localFilePath}/${params.localFilePattern}")
            .set { local_files }
        
        PROCESS_LOCAL_FILES(local_files.collect())
        
        fasta_files = PROCESS_LOCAL_FILES.out.processed_fasta
            .flatten()
            .filter { it.name != 'dummy.fa' }
            .map { file ->
                def sampleName = file.simpleName
                return tuple(sampleName, file)
            }
        fasta_files.view { "Processed file: $it" }
    } else {
        input_channel = Channel.fromPath(params.inputFile).splitText().map{it.trim()}
        DOWNLOAD_FOLDER(input_channel)
        fasta_files = DOWNLOAD_FOLDER.out.downloaded_fasta
            .flatten()
            .map { file -> 
                def sampleName = file.baseName
                return tuple(sampleName, file)
            }
    }
    
    FASTA_PREPROCESS_INITIAL(fasta_files)
    PROKKA_INITIAL(FASTA_PREPROCESS_INITIAL.out.cleaned_fasta)
    all_genomes = PROKKA_INITIAL.out.prokka_results
    selected_genome = all_genomes.randomSample(1).map { sampleName, gff, fna -> tuple(sampleName, gff, fna) }
    PANAROO_NON_MUTATED(PROKKA_INITIAL.out.gff_files.collect())
     
    EXTRACT_ALL_KEYWORD_GENES_NON_MUTATED(
        PANAROO_NON_MUTATED.out.panaroo_results.map { it -> file("${it}/gene_data.csv") },
        params.keyword
    )

    ALIGN_KEYWORD_GENES_TO_REFERENCE(
        EXTRACT_ALL_KEYWORD_GENES_NON_MUTATED.out.extracted_genes,
        PANAROO_NON_MUTATED.out.panaroo_results.map { it -> file("${it}/pan_genome_reference.fa") }
    )
    
    CREATE_CUSTOM_ACETYLTRANSFERASE_DB()

    BLAST_ALL_KEYWORD_GENES_NON_MUTATED(
        EXTRACT_ALL_KEYWORD_GENES_NON_MUTATED.out.extracted_genes,
        CREATE_CUSTOM_ACETYLTRANSFERASE_DB.out.db.collect()
    )
    
    COMPARE_BLAST_AND_GENE_INFO_NON_MUTATED(
        BLAST_ALL_KEYWORD_GENES_NON_MUTATED.out.blast_results,
        CREATE_CUSTOM_ACETYLTRANSFERASE_DB.out.gene_info,
        ALIGN_KEYWORD_GENES_TO_REFERENCE.out.blast_results,
        PANAROO_NON_MUTATED.out.panaroo_results.map { it -> file("${it}/gene_data.csv") },
        PANAROO_NON_MUTATED.out.panaroo_results.map { it -> file("${it}/pan_genome_reference.fa") }
    )
    
    ANALYZE_NON_MUTATED_GENES(EXTRACT_ALL_KEYWORD_GENES_NON_MUTATED.out.extracted_genes)

    // Collect all cleaned FASTA files into a single channel
    cleaned_fasta_files = FASTA_PREPROCESS_INITIAL.out.cleaned_fasta
        .map { it -> it[1] }
        .collect()

    // Additional processes for truncated gene analysis (only run if exTrun is true)
    if (params.exTrun) {
        SEARCH_TRUNCATED_GFF(
            COMPARE_BLAST_AND_GENE_INFO_NON_MUTATED.out.truncated_genes,
            cleaned_fasta_files
        )
    }
  
    // Additional processes (only run if coreOnly is false)
    if (!params.coreOnly) {
        // Always run DWGSIM_NON_MUTATED
        DWGSIM_NON_MUTATED(
            selected_genome
                .map { sampleName, gff, fna -> 
                    return tuple(sampleName, file(fna))
                }
        )
        
        DWGSIM_NON_MUTATED.out.simulated_reads
            .map { sampleName, read1, read2 -> tuple(sampleName, file(read1), file(read2)) }
            .set { dwgsim_non_mutated_reads }
        
        MUTATE_GENES(selected_genome)
        PROKKA_MUTATED(MUTATE_GENES.out.mutated_genome.map { it -> [it[0], it[2]] })
        
        DWGSIM_MUTATED(PROKKA_MUTATED.out.prokka_results.map { sampleName, gff, fna -> 
            tuple(sampleName, file(fna))
        })

        DWGSIM_MUTATED.out.simulated_reads
            .map { sampleName, read1, read2 -> tuple(sampleName, file(read1), file(read2)) }
            .set { dwgsim_reads_mutated }
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