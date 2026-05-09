include { MULTIQC } from '../../../modules/local/multiqc/main.nf'

workflow REPORTING {
    main:
    
    // Extract and format software_versions from the 'versions' topic
    ch_versions = channel.topic('versions')
        .map { proc, tool, ver ->
            def clean_ver = ver.toString().readLines()[0].trim().replaceAll(/^[vV]/,"")
            [ tool, clean_ver ]
        }
        .unique()
        .collectFile(
            name: 'software_versions_mqc.yml',
            sort: true,
            newLine: true,
            seed: "id: 'software_versions'\nsection_name: 'Software Versions'\nplot_type: 'table'\ndata:\n"
        ) { tool, ver -> " ${tool}:\n        version: '${ver}'\n"}

    // Collect all files from 'multiqc_files' topic
    // Use .collect() to ensure MultiQC qaits for all upstream processes to finish
    ch_multiqc_files = channel.topic('multiqc_files').collect()

    MULTIQC (
        ch_multiqc_files,
        ch_versions
    )

    emit:
    multiqc_report = MULTIQC.out.report
}