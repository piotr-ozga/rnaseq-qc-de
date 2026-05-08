process SRA_DOWNLOAD {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/sra-tools_pigz:8fcdefce490d03a1'
        : 'community.wave.seqera.io/library/sra-tools_pigz:a2ab3e08e0821853'}"

    secret 'NCBI_API_KEY'

    input:
    tuple val(meta), val(srr_id)

    output:
    tuple val(meta), path("*.fastq.gz"),     emit: reads
    tuple val("${task.process}"), val('fasterq-dump'), eval('fasterq-dump --version | grep -Eo "[0-9]+(\\.[0-9]+)+"'), emit: version, topic: versions

    script:
    def args_pf = task.ext.args_pf ?: ""
    def args_fd = task.ext.args_fd ?: ""
    """
    prefetch \\
        ${args_pf} \\
        --output-directory . \\
        ${srr_id}

    fasterq-dump \\
        ${args_fd} \\
        --threads ${task.cpus} \\
        --temp . \\
        --split-files \\
        ${srr_id}/${srr_id}.sra

    pigz -p ${task.cpus} *.fastq
    rm -rf ${srr_id}
    """

    stub:
    """
    touch ${srr_id}_1.fastq.gz
    ${meta.single_end ? "" : "touch ${srr_id}_2.fastq.gz"}
    """
}