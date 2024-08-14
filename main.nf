#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import modules
include { DOWNLOAD_FOLDER } from './modules/download_folder'
include { FASTA_PREPROCESS as FASTA_PREPROCESS_INITIAL } from './modules/fasta_preprocess'
include { FASTA_PREPROCESS as FASTA_PREPROCESS_ASSEMBLED } from './modules/fasta_preprocess'
include { FASTP } from './modules/fastp'
include { PROKKA as PROKKA_INITIAL } from './modules/prokka'
include { PROKKA as PROKKA_MUTATED } from './modules/prokka'
include { PROKKA as PROKKA_ASSEMBLED } from './modules/prokka'
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
include { ANALYZE_PREMATURE_STOPS_ASSEMBLED as ANALYZE_PREMATURE_STOPS_REAL } from './modules/analyze_premature_stops_assembled'
include { ANALYZE_PREMATURE_STOPS_ASSEMBLED as ANALYZE_PREMATURE_STOPS_PREASSEMBLED } from './modules/analyze_premature_stops_assembled'
include { BLAST_ACETYLTRANSFERASES as BLAST_ACETYLTRANSFERASES_REAL } from './modules/blast_acetyltransferases'
include { BLAST_ACETYLTRANSFERASES as BLAST_ACETYLTRANSFERASES_PREASSEMBLED } from './modules/blast_acetyltransferases'
include { UNICYCLER_ASSEMBLY } from './modules/unicycler_assembly'
include { APPLY_SNIPPY_CHANGES as APPLY_SNIPPY_CHANGES_MUTATED } from './modules/apply_snippy_changes'
include { APPLY_SNIPPY_CHANGES as APPLY_SNIPPY_CHANGES_NON_MUTATED } from './modules/apply_snippy_changes'
include { ANALYZE_BLAST_RESULTS as ANALYZE_BLAST_RESULTS_NON_MUTATED } from './modules/analyze_blast_results'
include { ANALYZE_BLAST_RESULTS as ANALYZE_BLAST_RESULTS_MUTATED } from './modules/analyze_blast_results'
include { ANALYZE_BLAST_RESULTS as ANALYZE_BLAST_RESULTS_REAL } from './modules/analyze_blast_results'
include { ANALYZE_BLAST_RESULTS as ANALYZE_BLAST_RESULTS_PREASSEMBLED } from './modules/analyze_blast_results'
include { EXTRACT_KEYWORD_GENES } from './modules/extract_keyword_genes'
include { APPLY_SNIPPY_CHANGES_REAL } from './modules/apply_snippy_changes_real'
include { COMBINE_REAL_RESULTS } from './modules/combine_real_results'

// Log info
log.info """
         P A N G E N O M E   A N A L Y S I S   P I P E L I N E    
         ===================================================
         Input file         : ${params.inputFile}
         Output directory   : ${params.outputDir}
         Non-mutated dir    : ${params.nonMutatedDir}
         Mutated dir        : ${params.mutatedDir}
         Use local files    : ${params.useLocalFiles}
         Local file path    : ${params.localFilePath}
         Core processes only: ${params.coreOnly}
         Run mutation analysis: ${params.run_mutation_analysis}
         
         Search term        : ${params.search_term}
         
         Mutation parameters:
         Num genes to mutate: ${params.num_genes_to_mutate}
         Mutation keyword   : ${params.mutation_keyword}
         Mutation mode      : ${params.mutation_mode}
         
         Gene extraction:
         First stage keyword : ${params.keyword}
         Second stage keyword: ${params.secondStageKeyword}
         
         DWGSIM parameters:
         Error rate         : ${params.dwgsimErrorRate}
         Outer distance     : ${params.dwgsimOuterDistance}
         Standard deviation : ${params.dwgsimStdDev}
         Number of reads    : ${params.dwgsimNumReads}
         Read length        : ${params.dwgsimReadLength}
         
         Snippy parameters:
         CPUs               : ${params.snippyCpus}
         
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
         Max target seqs    : ${params.blast_max_target_seqs}
         Word size          : ${params.blast_word_size}
         E-value            : ${params.blast_evalue}
         Percent identity   : ${params.blast_perc_identity}
         Dust               : ${params.blast_dust}
         Soft masking       : ${params.blast_soft_masking}
         X-drop ungapped    : ${params.blast_xdrop_ungap}
         X-drop gap         : ${params.blast_xdrop_gap}
         Gap open           : ${params.blast_gapopen}
         Gap extend         : ${params.blast_gapextend}
         
         Real reads analysis:
         Run real reads     : ${params.runRealReadsAnalysis}
         Use real reads     : ${params.useRealReads}
         Reads directory    : ${params.readsDir}
         Assembled genomes  : ${params.assembledGenomesDir}
         
         Acetyltransferase search:
         Search term        : ${params.acetyltransferase_search_term}
         Search retmax      : ${params.acetyltransferase_search_retmax}
         
         Unicycler parameters:
         Mode               : ${params.unicycler_mode}
         Min fasta length   : ${params.unicycler_min_fasta_length}
         Kmer count         : ${params.unicycler_kmer_count}
         Min component size : ${params.unicycler_min_component_size}
         Min dead end size  : ${params.unicycler_min_dead_end_size}
         Min bridge qual    : ${params.unicycler_min_bridge_qual}
         Keep temp          : ${params.unicycler_keep_temp}
         No correct         : ${params.unicycler_no_correct}
         Verbosity          : ${params.unicycler_verbosity}
         
         FASTA preprocessing:
         Min length         : ${params.fasta_preprocess.min_length}
         Remove duplicates  : ${params.fasta_preprocess.remove_duplicates}
         
         FASTP parameters:
         Qualified quality  : ${params.fastp.qualified_quality_phred}
         Unqualified limit  : ${params.fastp.unqualified_percent_limit}
         Cut mean quality   : ${params.fastp.cut_mean_quality}
         Cut window size    : ${params.fastp.cut_window_size}
         Cut front          : ${params.fastp.cut_front}
         Cut tail           : ${params.fastp.cut_tail}
         Cut front window   : ${params.fastp.cut_front_window_size}
         Cut front quality  : ${params.fastp.cut_front_mean_quality}
         Cut tail window    : ${params.fastp.cut_tail_window_size}
         Cut tail quality   : ${params.fastp.cut_tail_mean_quality}
         Length required    : ${params.fastp.length_required}
         
"""

// Parameter validation
if (!params.inputFile || !params.outputDir) {
    error "Input file and output directory parameters are required: --inputFile and --outputDir"
}

// Main workflow
workflow {
    // Core processes
    if (params.useLocalFiles) {
        PROCESS_LOCAL_FILES(params.localFilePath)
        fasta_files = PROCESS_LOCAL_FILES.out.processed_fasta.flatten().map { file ->
            def sampleName = file.baseName
            return tuple(sampleName, file)
        }
    } else {
        input_channel = Channel.fromPath(params.inputFile).splitText().map{it.trim()}
        DOWNLOAD_FOLDER(input_channel)
        fasta_files = DOWNLOAD_FOLDER.out.downloaded_fasta.flatten().map { file -> 
            def sampleName = file.baseName
            return tuple(sampleName, file)
        }
    }
    
    FASTA_PREPROCESS_INITIAL(fasta_files)
    PROKKA_INITIAL(FASTA_PREPROCESS_INITIAL.out.cleaned_fasta)
    all_genomes = PROKKA_INITIAL.out.prokka_results
    selected_genome = all_genomes.randomSample(1).map { sampleName, gff, fna -> tuple(sampleName, gff, fna) }
    PANAROO_NON_MUTATED(PROKKA_INITIAL.out.gff_files.collect())
    
    EXTRACT_NON_MUTATED_GENOME_KEYWORD_GENES(
        PANAROO_NON_MUTATED.out.panaroo_results
            .combine(selected_genome)
            .map { panaroo_result, sampleName, gff, fna -> 
                tuple(sampleName, file("${panaroo_result}/gene_data.csv"))
            },
        params.keyword
    )
    
    ALIGN_KEYWORD_GENES_TO_REFERENCE(
        EXTRACT_NON_MUTATED_GENOME_KEYWORD_GENES.out.extracted_genes,
        PANAROO_NON_MUTATED.out.panaroo_results.map { it -> file("${it}/pan_genome_reference.fa") }
    )
    
    CREATE_CUSTOM_ACETYLTRANSFERASE_DB()
    
    BLAST_ACETYLTRANSFERASES_NON_MUTATED(
        EXTRACT_NON_MUTATED_GENOME_KEYWORD_GENES.out.extracted_genes,
        CREATE_CUSTOM_ACETYLTRANSFERASE_DB.out.db.collect()
    )
    
    ANALYZE_NON_MUTATED_GENES(EXTRACT_NON_MUTATED_GENOME_KEYWORD_GENES.out.extracted_genes)
    
    ANALYZE_BLAST_RESULTS_NON_MUTATED(
        BLAST_ACETYLTRANSFERASES_NON_MUTATED.out.blast_results
            .combine(EXTRACT_NON_MUTATED_GENOME_KEYWORD_GENES.out.extracted_genes)
            .map { blast_result, extracted_genes -> 
                def sampleName = blast_result.baseName.replaceFirst(/_blast_results$/, '')
                return tuple(sampleName, blast_result, extracted_genes)
            },
        CREATE_CUSTOM_ACETYLTRANSFERASE_DB.out.db.collect()
    )

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
    
    EXTRACT_ALL_KEYWORD_GENES_NON_MUTATED(
        PANAROO_NON_MUTATED.out.panaroo_results.map { it -> file("${it}/gene_data.csv") },
        params.keyword
    )

    MUTATE_GENES(selected_genome)
    PROKKA_MUTATED(MUTATE_GENES.out.mutated_genome.map { it -> [it[0], it[2]] })
    
    DWGSIM_MUTATED(PROKKA_MUTATED.out.prokka_results.map { sampleName, gff, fna -> 
        tuple(sampleName, file(fna))
    })

    DWGSIM_MUTATED.out.simulated_reads
        .map { sampleName, read1, read2 -> tuple(sampleName, file(read1), file(read2)) }
        .set { dwgsim_reads_mutated }
    
    // Additional processes (only run if coreOnly is false)
    if (!params.coreOnly) {
        mutated_gff = PROKKA_MUTATED.out.prokka_results.map { it[1] }
        
        mutated_base_names = mutated_gff.map { file -> 
            def base_name = file.name.replaceFirst(/_mutated\.gff$/, '')
            return tuple(base_name, file)
        }
        
        all_gffs_for_mutated = PROKKA_INITIAL.out.gff_files
            .map { file -> 
                def base_name = file.name.replaceFirst(/\.gff$/, '')
                return tuple(base_name, file)
            }
            .join(mutated_base_names, remainder: true)
            .map { base_name, non_mutated, mutated -> 
                mutated ?: non_mutated
            }
        
        PANAROO_MUTATED(all_gffs_for_mutated.collect())
        
        EXTRACT_MUTATED_GENOME_KEYWORD_GENES(
            PANAROO_MUTATED.out.panaroo_results
                .combine(selected_genome)
                .map { panaroo_result, selectedSampleName, selectedGff, selectedFna -> 
                    tuple(selectedSampleName, file("${panaroo_result}/gene_data.csv"))
                },
            params.keyword
        )
        
        BLAST_ACETYLTRANSFERASES_MUTATED(
            EXTRACT_MUTATED_GENOME_KEYWORD_GENES.out.extracted_genes,
            CREATE_CUSTOM_ACETYLTRANSFERASE_DB.out.db.collect()
        )
        
        ANALYZE_BLAST_RESULTS_MUTATED(
            BLAST_ACETYLTRANSFERASES_MUTATED.out.blast_results
                .combine(EXTRACT_MUTATED_GENOME_KEYWORD_GENES.out.extracted_genes)
                .map { blast_result, extracted_genes -> 
                    def sampleName = blast_result.baseName.replaceFirst(/_blast_results$/, '')
                    return tuple(sampleName, blast_result, extracted_genes)
                },
            CREATE_CUSTOM_ACETYLTRANSFERASE_DB.out.db.collect()
        )
        
        EXTRACT_ALL_KEYWORD_GENES_MUTATED(
            PANAROO_MUTATED.out.panaroo_results.map { it -> file("${it}/gene_data.csv") },
            params.keyword
        )
        
        ALIGN_KEYWORD_GENES_TO_REFERENCE_MUTATED(
            EXTRACT_ALL_KEYWORD_GENES_MUTATED.out.extracted_genes_mutated,
            PANAROO_MUTATED.out.panaroo_results.map { it -> file("${it}/pan_genome_reference.fa") }
        )
        
        ANALYZE_MUTATED_GENES(EXTRACT_MUTATED_GENOME_KEYWORD_GENES.out.extracted_genes)
        
        ANALYZE_NON_MUTATED_GENES.out.analysis_results
            .mix(ANALYZE_MUTATED_GENES.out.analysis_results)
            .collectFile(name: "${params.outputDir}/combined_analysis_results.json", newLine: true, sort: true)
        
        SNIPPY_MUTATED(
            dwgsim_reads_mutated,
            EXTRACT_NON_MUTATED_GENOME_KEYWORD_GENES.out.extracted_genes
        )
        
        APPLY_SNIPPY_CHANGES_MUTATED(
            SNIPPY_MUTATED.out.snippy_results
        )
        
        SNIPPY_NON_MUTATED(
            dwgsim_non_mutated_reads,
            EXTRACT_MUTATED_GENOME_KEYWORD_GENES.out.extracted_genes
        )
        
        APPLY_SNIPPY_CHANGES_NON_MUTATED(
            SNIPPY_NON_MUTATED.out.snippy_results
        )
    }

    // Real reads analysis
    if (params.runRealReadsAnalysis) {
        if (params.useRealReads && params.readsDir) {
            real_reads_ch = Channel.fromFilePairs("${params.readsDir}/*_{1,2}.fastq.gz")
            reference_genome = PANAROO_NON_MUTATED.out.pangenome_fasta
        } else {
            real_reads_ch = DWGSIM_MUTATED.out.simulated_reads
            reference_genome = selected_genome.map { sampleName, gff, fna -> file(fna) }
        }
        
        FASTP(real_reads_ch)
        cleaned_reads = FASTP.out.cleaned_reads
        
        UNICYCLER_ASSEMBLY(cleaned_reads)
        PROKKA_ASSEMBLED(UNICYCLER_ASSEMBLY.out.assembly)
        
        EXTRACT_ACETYLTRANSFERASE_GENES_GFF_REAL(
            PROKKA_ASSEMBLED.out.prokka_results
                .map { sample_id, gff, fna -> tuple(sample_id, gff, fna) }, 
            params.secondStageKeyword
        )
        
        EXTRACT_KEYWORD_GENES(PANAROO_NON_MUTATED.out.panaroo_results, params.keyword)

        snippy_input = FASTP.out.cleaned_reads
            .combine(EXTRACT_KEYWORD_GENES.out.extracted_genes)
            .map { sampleName, read1, read2, extracted_genes -> 
                tuple(sampleName, read1, read2, extracted_genes)
            }

        SNIPPY_REAL_READS(snippy_input)

        APPLY_SNIPPY_CHANGES_REAL(SNIPPY_REAL_READS.out.snippy_results)
        
        ANALYZE_PREMATURE_STOPS_REAL(EXTRACT_ACETYLTRANSFERASE_GENES_GFF_REAL.out.extracted_genes)
        
        BLAST_ACETYLTRANSFERASES_REAL(EXTRACT_ACETYLTRANSFERASE_GENES_GFF_REAL.out.extracted_genes, CREATE_CUSTOM_ACETYLTRANSFERASE_DB.out.db.collect())
        
        ANALYZE_BLAST_RESULTS_REAL(
            BLAST_ACETYLTRANSFERASES_REAL.out.blast_results
                .join(EXTRACT_ACETYLTRANSFERASE_GENES_GFF_REAL.out.extracted_genes)
                .map { sampleName, blast_result, extracted_genes -> 
                    return tuple(sampleName, blast_result, extracted_genes)
                },
            CREATE_CUSTOM_ACETYLTRANSFERASE_DB.out.db.collect()
        )

        COMBINE_REAL_RESULTS(
            APPLY_SNIPPY_CHANGES_REAL.out.analysis_results,
            BLAST_ACETYLTRANSFERASES_REAL.out.blast_results,
            EXTRACT_ACETYLTRANSFERASE_GENES_GFF_REAL.out.extracted_genes
        )
    } else if (params.assembledGenomesDir) {
        assembled_genomes_ch = Channel.fromPath("${params.assembledGenomesDir}/*.fa")
            .map { file -> tuple(file.baseName, file) }
        
        FASTA_PREPROCESS_ASSEMBLED(assembled_genomes_ch)
        PROKKA_ASSEMBLED(FASTA_PREPROCESS_ASSEMBLED.out.cleaned_fasta)
        
        EXTRACT_ACETYLTRANSFERASE_GENES_GFF_ASSEMBLED(
            PROKKA_ASSEMBLED.out.prokka_results
                .map { sample_id, gff, fna -> tuple(sample_id, gff, fna) }, 
            params.secondStageKeyword
        )
        
        ANALYZE_PREMATURE_STOPS_PREASSEMBLED(EXTRACT_ACETYLTRANSFERASE_GENES_GFF_ASSEMBLED.out.extracted_genes)
        
        BLAST_ACETYLTRANSFERASES_PREASSEMBLED(EXTRACT_ACETYLTRANSFERASE_GENES_GFF_ASSEMBLED.out.extracted_genes, CREATE_CUSTOM_ACETYLTRANSFERASE_DB.out.db.collect())
        
        ANALYZE_BLAST_RESULTS_PREASSEMBLED(
            BLAST_ACETYLTRANSFERASES_PREASSEMBLED.out.blast_results
                .join(EXTRACT_ACETYLTRANSFERASE_GENES_GFF_ASSEMBLED.out.extracted_genes)
                .map { sampleName, blast_result, extracted_genes -> 
                    return tuple(sampleName, blast_result, extracted_genes)
                },
            CREATE_CUSTOM_ACETYLTRANSFERASE_DB.out.db.collect()
        )
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
