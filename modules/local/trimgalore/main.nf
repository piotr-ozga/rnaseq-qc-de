process TRIMGALORE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/trim-galore:2.1.0--9c3d977448fceb85'
        : 'community.wave.seqera.io/library/trim-galore:2.1.0--27e6376b8f6c1872'}"


    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.fq.gz"),         emit: reads
    tuple val(meta), path("*report.txt"),     emit: log
    tuple val("${task.process}"), val('trim_galore'), eval('trim_galore --version | grep -Eo "[0-9]+(\\.[0-9]+)+"'), emit: version, topic: versions

    script:
    def args = task.ext.args ?: ''
    
    if (args.contains('--fastqc')) {
        log.warn "[MODULE: TRIMGALORE] Sample: ${meta.id} > '--fastqc' from ext.args ignored! This pipeline uses a dedicated FASTQC module for reporting."
    }

    def filtered_args = args.tokenize().findAll { it != "--fastqc" }.join(' ')
    def paired = meta.single_end ? "" : "--paired"
    """
    trim_galore \\
        ${filtered_args} \\
        ${paired} \\
        --gzip \\
        --cores ${task.cpus} \\
        ${reads}
    """

    stub:
    if (meta.single_end) {
        """
        touch ${meta.id}_trimmed.fq.gz
        touch ${meta.id}_trimming_report.txt
        """
    } else {
        """
        touch ${meta.id}_val_1.fq.gz
        touch ${meta.id}_val_2.fq.gz
        touch ${meta.id}_trimming_report_1.txt
        touch ${meta.id}_trimming_report_2.txt
        """
    }
}