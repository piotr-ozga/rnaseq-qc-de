process FASTQC {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/fastqc:0.12.1--104d26ddd9519960' :
        'community.wave.seqera.io/library/fastqc:0.12.1--af7a5314d5015c29' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"),     emit: html
    tuple val(meta), path("*.zip"),      emit: zip
    tuple val("${task.process}"), val('fastqc'), eval("fastqc --version | sed 's/FastQC v//'"), emit: version, topic: versions

    script:
    def args = task.ext.args ?: ""
    """
    fastqc \\
        ${args} \\
        --threads ${task.cpus} \\
        ${reads}
    """

    stub:
    """
    touch ${ [reads].flatten().collect { it.simpleName + "_fastqc.{html,zip}" }.join(' ') }
    """
}