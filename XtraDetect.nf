#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

// Import modules
include { FASTP } from './modules/fastp'
include { UNICYCLER_ASSEMBLY } from './modules/unicycler_assembly'
include { PROKKA as PROKKA_ASSEMBLED } from './modules/prokka'
include { EXTRACT_ACETYLTRANSFERASE_GENES_GFF as EXTRACT_ACETYLTRANSFERASE_GENES_GFF_REAL } from './modules/extract_acetyltransferase_genes_gff'
include { EXTRACT_KEYWORD_GENES } from './modules/extract_keyword_genes'
include { SNIPPY_REAL_READS } from './modules/snippy_real_reads'
include { APPLY_SNIPPY_CHANGES_REAL } from './modules/apply_snippy_changes_real'
include { ANALYZE_PREMATURE_STOPS_ASSEMBLED as ANALYZE_PREMATURE_STOPS_REAL } from './modules/analyze_premature_stops_assembled'
include { BLAST_ACETYLTRANSFERASES as BLAST_ACETYLTRANSFERASES_REAL } from './modules/blast_acetyltransferases'
include { FASTA_PREPROCESS as FASTA_PREPROCESS_ASSEMBLED } from './modules/fasta_preprocess'
include { EXTRACT_ACETYLTRANSFERASE_GENES_GFF as EXTRACT_ACETYLTRANSFERASE_GENES_GFF_ASSEMBLED } from './modules/extract_acetyltransferase_genes_gff'
include { ANALYZE_PREMATURE_STOPS_ASSEMBLED as ANALYZE_PREMATURE_STOPS_PREASSEMBLED } from './modules/analyze_premature_stops_assembled'
include { BLAST_ACETYLTRANSFERASES as BLAST_ACETYLTRANSFERASES_PREASSEMBLED } from './modules/blast_acetyltransferases'
include { GATHER_TOP_BLAST_INFO_AND_ANNOTATIONS as GATHER_TOP_BLAST_INFO_AND_ANNOTATIONS_READS} from './modules/gather_top_blast_info_and_annotations'
include { GATHER_TOP_BLAST_INFO_AND_ANNOTATIONS as GATHER_TOP_BLAST_INFO_AND_ANNOTATIONS_ASSEMBLE} from './modules/gather_top_blast_info_and_annotations'
include { ALIGN_KEYWORD_GENES_TO_REFERENCE as ALIGN_ACETYLTRANSFERASES_TO_PANGENOME } from './modules/align_keyword_genes_to_reference'
include { DWGSIM_NON_MUTATED } from './modules/dwgsim_non_mutated'

// Main workflow
workflow {
    // Create channels for input files
    if (params.fa) {
        input_ch = Channel.fromPath(params.inputfa)
            .map { fa -> tuple(fa.simpleName, fa) }
    } else if (params.reads != null) {
        def patterns = [
            "*_{1,2}.{fastq,fq}.{gz,bz2}",
            "*_R{1,2}.{fastq,fq}.{gz,bz2}",
            "*_{1,2}.{fastq,fq}",
            "*_R{1,2}.{fastq,fq}",
            "*.{1,2}.{fastq,fq}.{gz,bz2}",
            "*.R{1,2}.{fastq,fq}.{gz,bz2}",
            "*.{1,2}.{fastq,fq}",
            "*.R{1,2}.{fastq,fq}",
            "*#*_{1,2}.fastq.gz"
        ]
        
        input_ch = Channel.empty()
        for (pattern in patterns) {
            def files = file("${params.inputreads}/${pattern}")
            if (files) {
                input_ch = Channel.fromFilePairs("${params.inputreads}/${pattern}", checkIfExists: true)
                    .map { sample_id, reads -> tuple(sample_id, reads[0], reads[1]) }
                break
            }
        }
        
        input_ch.ifEmpty { error "No files match any of the patterns in directory: ${params.inputreads}" }
    } else {
        error "Either '--fa' or '--reads' must be specified."
    }

    // Process each sample
    input_ch.each { sample ->
        PROCESS_SAMPLE(sample)
    }
}

// Sub-workflow to process a single sample
workflow PROCESS_SAMPLE {
    take:
    sample

    main:
    if (params.fa) {
        // Run DWGSIM on .fa files
        DWGSIM_NON_MUTATED(sample)
        real_reads_ch = DWGSIM_NON_MUTATED.out.simulated_reads
        input_fa = sample
    } else if (params.reads) {
        real_reads_ch = sample
        input_fa = Channel.empty()
    } else {
        error "Neither 'fa' nor 'reads' parameter is set. Please provide either --fa or --reads."
    }

    // Real reads analysis
    FASTP(real_reads_ch)
    cleaned_reads = FASTP.out.cleaned_reads
    
    if (!params.skipUnicycler) {
        UNICYCLER_ASSEMBLY(cleaned_reads)
        assembly_ch = UNICYCLER_ASSEMBLY.out.assembly
    } else {
        assembly_ch = input_fa
    }
    
    PROKKA_ASSEMBLED(assembly_ch)
    
    EXTRACT_ACETYLTRANSFERASE_GENES_GFF_REAL(
        PROKKA_ASSEMBLED.out.prokka_results
            .map { id, gff, fna -> tuple(id, gff, fna) }, 
        params.secondStageKeyword
    )
    
    panaroo_results = Channel.fromPath("${params.outputDir}/panaroo_non_mutated/panaroo_output")
    EXTRACT_KEYWORD_GENES(panaroo_results, params.keyword)

    // Define the pangenome reference file
    pangenome_reference = file("${params.outputDir}/panaroo_non_mutated/panaroo_output/pan_genome_reference.fa")

    ALIGN_ACETYLTRANSFERASES_TO_PANGENOME(
        EXTRACT_ACETYLTRANSFERASE_GENES_GFF_REAL.out.extracted_genes.map { id, fasta -> fasta },
        pangenome_reference
    )

    pangenome_reference_ch = Channel.fromPath(pangenome_reference)

    snippy_input = FASTP.out.cleaned_reads
        .combine(pangenome_reference_ch)
        .map { id, r1, r2, pangenome_reference -> 
            tuple(id, r1, r2, pangenome_reference)
        }

    SNIPPY_REAL_READS(snippy_input)

    APPLY_SNIPPY_CHANGES_REAL(SNIPPY_REAL_READS.out.snippy_results)

    ANALYZE_PREMATURE_STOPS_REAL(EXTRACT_ACETYLTRANSFERASE_GENES_GFF_REAL.out.extracted_genes)
    
    custom_db = Channel.fromPath("${params.outputDir}/Custom_BLAST_DB/*")
    BLAST_ACETYLTRANSFERASES_REAL(
        EXTRACT_ACETYLTRANSFERASE_GENES_GFF_REAL.out.extracted_genes,
        custom_db.collect()
    )
    
    gene_info = file("${params.outputDir}/Custom_BLAST_DB/gene_info.tsv")

    // Then, update the GATHER_TOP_BLAST_INFO_AND_ANNOTATIONS_READS process call
    GATHER_TOP_BLAST_INFO_AND_ANNOTATIONS_READS(
        BLAST_ACETYLTRANSFERASES_REAL.out.blast_results,
        EXTRACT_ACETYLTRANSFERASE_GENES_GFF_REAL.out.extracted_genes,
        ALIGN_ACETYLTRANSFERASES_TO_PANGENOME.out.blast_results,
        APPLY_SNIPPY_CHANGES_REAL.out.analysis_results.map { it -> it[1] },
        PROKKA_ASSEMBLED.out.prokka_results.map { it -> it[1] },
        pangenome_reference,
        gene_info
    )
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