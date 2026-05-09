process MULTIQC {
    tag "report"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/multiqc:1.34--4fc8657c816047c0' :
        'community.wave.seqera.io/library/multiqc:1.34--db7c73dae76bc9e6' }"

    input:
    path 'inputs/*'
    path 'software_versions.yml'

    output:
    path "multiqc_report.html",     emit: report
    path "multiqc_data",            emit: data

    script:
    def args = task.ext.args ?: ""
    """
    multiqc . \\
        --comment "rnaseq-qc-de Pipeline Report" \\
        ${args} \\
    """

    stub:
    """
    touch multiqc_report.html
    mkdir multiqc_data
    """
}