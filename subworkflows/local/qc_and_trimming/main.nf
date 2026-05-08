include { FASTQC as FASTQC_RAW} from '../../../modules/local/fastqc/main.nf'
workflow QC_AND_TRIMMING {
    take:
    reads

    main:
    FASTQC_RAW(reads)

    emit:
    fastqc_zip_raw        = FASTQC_RAW.out.zip
}