include { SRA_DOWNLOAD } from '../../../modules/local/sratools/main.nf'

workflow INPUT_HANDLER {
    take:
    ch_samples

    main:
    ch_branched = ch_samples.branch {
        local:  it[2] != null
        remote: it[2] == null
    }

    ch_local = ch_branched.local
        .map { meta, srr, r1, r2 ->
            def files = r2 ? [ r1, r2 ] : [ r1 ]
            [ meta, files ]
        }
    
    SRA_DOWNLOAD(ch_branched.remote.map { 
        meta, srr, r1, r2 -> [ meta, srr ] }
    )

    ch_reads = ch_local.mix(SRA_DOWNLOAD.out.reads)

    emit:
    reads    = ch_reads
    versions = SRA_DOWNLOAD.out.version
}