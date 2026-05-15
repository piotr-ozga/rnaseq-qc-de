include { DESEQ2_ANALYSIS } from '../../../modules/local/deseq2/main.nf'
include { ANNOTATE_RESULTS } from '../../../modules/local/annotation/main.nf'
include { VISUALIZE_DEA } from '../../../modules/local/visualization/main.nf'

// Runs DESeq2 differential expression analysis, annotates results with
// GTF gene information, and generates visualization plots
workflow RNASEQ_DEA {
    take:
    salmon_dir
    samplesheet
    reference_level
    gtf
    lfc_threshold
    padj_threshold
    heatmap_genes
    volcano_labels

    main:
    DESEQ2_ANALYSIS(
        salmon_dir,
        samplesheet,
        reference_level,
        lfc_threshold,
        padj_threshold
    )

    ANNOTATE_RESULTS(
        DESEQ2_ANALYSIS.out.results_tsv,
        gtf
    )

    VISUALIZE_DEA(
        ANNOTATE_RESULTS.out.annotated_tsv,
        DESEQ2_ANALYSIS.out.vst_rds,
        samplesheet,
        reference_level,
        lfc_threshold,
        padj_threshold,
        heatmap_genes,
        volcano_labels
    )

    emit:
    results_tsv        = DESEQ2_ANALYSIS.out.results_tsv
    dds_rds            = DESEQ2_ANALYSIS.out.dds_rds
    vst_rds            = DESEQ2_ANALYSIS.out.vst_rds
    annotated_tsv      = ANNOTATE_RESULTS.out.annotated_tsv
    annotation_summary = ANNOTATE_RESULTS.out.summary_txt
    plots              = VISUALIZE_DEA.out.plots
}
