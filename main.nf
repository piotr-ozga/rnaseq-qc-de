#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { INPUT_HANDLER } from './subworkflows/local/input_handler/main.nf' 
include { QC_AND_TRIMMING } from './subworkflows/local/qc_and_trimming/main.nf' 
include { SALMON_ALIGN } from './subworkflows/local/salmon_align/main.nf'
include { REPORTING } from './subworkflows/local/reporting/main.nf' 
include { RNASEQ_DEA } from './subworkflows/local/rnaseq_dea/main.nf'

workflow {
    // --- PARAMETER VALIDATION --- 
    def mandatory_files = [
        "Samplesheet": params.samplesheet,
        "Genome FASTA": params.genome_fasta,
        "GTF Annotation": params.gtf,
        "Transcriptome FASTA": params.transcriptome_fasta
    ]

    mandatory_files.each { name, path ->
        if (!path) {
            error "[ERROR]: ${name} is not defined. Please use --${name.toLowerCase().replace(' ', '_')} <path>"
        }
        if (!file(path).exists()) {
            error "[ERROR]: ${name} file not found: ${path}"
        }
    }

    // --- Pipeline summary ---
    log.info """
        ================================================
        R N A S E Q - Q C - D E  P I P E L I N E
        ================================================
        outdir          : ${params.outdir}
        samplesheet     : ${params.samplesheet}
        ================================================
        """.stripIndent()

    // --- INPUT CHANNEL ---
    ch_samples = channel
        .fromPath(params.samplesheet)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def meta = [ id: row.sample ]
            meta.single_end = row.read1 && !row.read2 ? true : false

            def r1 = row.read1 ? file(row.read1) : null
            def r2 = row.read2 ? file(row.read2) : null
            [ meta, row.srr, r1, r2 ]
        }

    INPUT_HANDLER(ch_samples)
    QC_AND_TRIMMING(INPUT_HANDLER.out.reads)
    SALMON_ALIGN(
        QC_AND_TRIMMING.out.trimmed_reads,
        file(params.genome_fasta),
        file(params.transcriptome_fasta),
        file(params.gtf),
        params.salmon_lib_type
    )

    RNASEQ_DEA(
        SALMON_ALIGN.out.results
            .map { meta, dir -> dir }
            .collect(),
        file(params.samplesheet),
        params.ref_level
    )

    REPORTING()
}