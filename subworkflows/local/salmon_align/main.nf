include { SALMON_INDEX } from '../../../modules/local/salmon/index/main.nf'
include { SALMON_QUANT } from '../../../modules/local/salmon/quant/main.nf'

workflow SALMON_ALIGN {
    take:
    reads
    genome_fasta
    transcript_fasta
    gtf
    lib_type

    main:
    SALMON_INDEX(genome_fasta, transcript_fasta)

    // collect() broadcast single index to all samples
    SALMON_QUANT(
        reads,
        SALMON_INDEX.out.index.collect(),
        gtf,
        lib_type
    )

    emit:
    results = SALMON_QUANT.out.results
}
