process SALMON_INDEX {
    tag "salmon_index"
    label "process_medium"

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/salmon:1.10.3--726401738a281398'
        : 'community.wave.seqera.io/library/salmon:1.10.3--fcd0755dd8abb423'}"

    input:
    path genome_fasta
    path transcript_fasta

    output:
    path "salmon_index",     emit: index
    tuple val("${task.process}"), val('salmon'), eval('salmon --version | grep -Eo "[0-9]+(\\.[0-9]+)+"'), emit: version, topic: versions

    script:
    def args = task.ext.args ?: ""
    """
    zgrep "^>" ${genome_fasta} | cut -d " " -f 1 | sed 's/>//g' > decoys.txt
    zcat -f ${transcript_fasta} > gentrome.fa
    zcat -f ${genome_fasta}     >> gentrome.fa

    salmon index \\
        ${args} \\
        --threads ${task.cpus} \\
        -t gentrome.fa \\
        -d decoys.txt \\
        -i salmon_index
    """

    stub:
    """
    mkdir salmon_index
    """
}