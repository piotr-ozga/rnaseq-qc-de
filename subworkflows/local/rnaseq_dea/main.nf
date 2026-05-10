include { DESEQ2_ANALYSIS } from '../../../modules/local/deseq2/main.nf'

workflow RNASEQ_DEA {
    take:
    salmon_dir
    samplesheet
    reference_level

    main:
    DESEQ2_ANALYSIS(salmon_dir, samplesheet, reference_level)

    emit:
    results_tsv = DESEQ2_ANALYSIS.out.results_tsv
    dds_rds = DESEQ2_ANALYSIS.out.dds_rds
    vst_rds = DESEQ2_ANALYSIS.out.vst_rds
}
