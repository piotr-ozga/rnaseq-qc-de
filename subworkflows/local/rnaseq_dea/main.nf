include { DESEQ2_ANALYSIS } from '../../../modules/local/deseq2/main.nf'
include { ANNOTATE_RESULTS } from '../../../modules/local/annotation/main.nf'

workflow RNASEQ_DEA {
    take:
    salmon_dir
    samplesheet
    reference_level
    gtf

    main:
    // Perform DEA
    DESEQ2_ANALYSIS(salmon_dir, samplesheet, reference_level)
    // Add annotations with GTF
    ANNOTATE_RESULTS(
        DESEQ2_ANALYSIS.out.results_tsv,
        gtf
    )

    emit:
    results_tsv        = DESEQ2_ANALYSIS.out.results_tsv
    dds_rds            = DESEQ2_ANALYSIS.out.dds_rds
    vst_rds            = DESEQ2_ANALYSIS.out.vst_rds
    annotated_tsv      = ANNOTATE_RESULTS.out.annotated_tsv
    annotation_summary = ANNOTATE_RESULTS.out.summary_txt 
}
