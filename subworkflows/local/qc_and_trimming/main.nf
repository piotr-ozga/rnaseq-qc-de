include { FASTQC as FASTQC_RAW} from '../../../modules/local/fastqc/main.nf'
include { TRIMGALORE } from '../../../modules/local/trimgalore/main.nf'
include { FASTQC as FASTQC_TRIMMED} from '../../../modules/local/fastqc/main.nf'

workflow QC_AND_TRIMMING {
    take:
    reads

    main:
    FASTQC_RAW(reads)
    TRIMGALORE(reads)
    FASTQC_TRIMMED(TRIMGALORE.out.reads)

    emit:
    trimmed_reads         = TRIMGALORE.out.reads
    fastqc_zip_raw        = FASTQC_RAW.out.zip
    fastqc_zip_trimmed    = FASTQC_TRIMMED.out.zip
    trim_log              = TRIMGALORE.out.log
}