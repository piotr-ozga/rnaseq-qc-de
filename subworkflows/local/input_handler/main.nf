include { SRA_DOWNLOAD } from '../../../modules/local/sratools/main.nf'

workflow INPUT_HANDLER {
    take:
    samplesheet

    main:
    // Input channel
    ch_samples = channel
        .fromPath(samplesheet)
        .splitCsv(header: true, sep: ',')
        .map { row ->
            def meta = [ 
                id:        row.sample,
                condition: row.condition
            ]
            meta.single_end = row.read1 && !row.read2 ? true : false

            def r1 = row.read1 ? file("${workflow.projectDir}/${row.read1}", checkIfExists: true) : null
            def r2 = row.read2 ? file("${workflow.projectDir}/${row.read2}", checkIfExists: true) : null
            [ meta, row.srr, r1, r2 ]
        }

    // Validate samplesheet integrity
    ch_validated = ch_samples.map { meta, srr, r1, r2 ->
        if (!meta.id) {
            error "[SAMPLESHEET ERROR]: Missing 'sample' ID for row: [${srr}, ${r1}, ${r2}]" 
        }
        if (!meta.condition) {
            error "[SAMPLESHEET ERROR]: Sample '${meta.id}' must provide the 'condition' column value." 
        }
        if (!srr && !r1) {
            error "[SAMPLESHEET ERROR]: Sample '${meta.id}' must provide either an SRR ID or local read paths." 
        }
        
        return [ meta, srr, r1, r2 ]
    }

    // Split samples into local files and remote SRA accessions
    ch_branched = ch_validated.branch {
        local:  it[2] != null
        remote: it[2] == null
    }

    // Format local paths into a list of files
    ch_local = ch_branched.local
        .map { meta, srr, r1, r2 ->
            def files = r2 ? [ r1, r2 ] : [ r1 ]
            [ meta, files ]
        }
    
    // Download raw FASTQ files from NCBI SRA for remote accession
    SRA_DOWNLOAD(ch_branched.remote.map { 
        meta, srr, r1, r2 -> [ meta, srr ] }
    )

    // Merge local and downloaded reads into a single stream for processing
    ch_reads = ch_local.mix(SRA_DOWNLOAD.out.reads)

    emit:
    reads    = ch_reads
    versions = SRA_DOWNLOAD.out.version
}