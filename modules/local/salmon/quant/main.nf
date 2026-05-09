process SALMON_QUANT {
    tag "${meta.id}"
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/salmon:1.10.3--726401738a281398'
        : 'community.wave.seqera.io/library/salmon:1.10.3--fcd0755dd8abb423'}"

    input:
    tuple val(meta), path(reads)
    path index
    path gtf
    val lib_type

    output:
    tuple val(meta), path("${meta.id}"),     emit: results
    path("${meta.id}"),                      emit: mqc, topic: multiqc_files
    tuple val("${task.process}"), val('salmon'), eval('salmon --version | grep -Eo "[0-9]+(\\.[0-9]+)+"'), emit: version, topic: versions

    script:
    def input_reads = meta.single_end ? "-r ${reads}" : "-1 ${reads[0]} -2 ${reads[1]}"
    def args = task.ext.args ?: ""
    """
    salmon quant \\
        ${args} \\
        --threads ${task.cpus} \\
        --libType ${lib_type} \\
        --index ${index} \\
        --geneMap ${gtf} \\
        ${input_reads} \\
        --validateMappings \\
        -o ${meta.id}
     """

    stub:
    """
    mkdir ${meta.id}
    touch ${meta.id}/quant.sf
    touch ${meta.id}/quant.genes.sf
    """
}