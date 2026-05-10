process DESEQ2_ANALYSIS {
    tag "DESeq2 on ${samplesheet.baseName}"
    label "process_low"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-deseq2_bioconductor-tximport_r-base_r-dplyr_pruned:2e8147b24a0a2537' :
        'community.wave.seqera.io/library/bioconductor-deseq2_bioconductor-tximport_r-base_r-dplyr_pruned:3aca25f1b464c8f1' }"

    input:
    path "inputs/*"
    path samplesheet
    val reference_level

    output:
    path "results/results.tsv", emit: results_tsv
    path "results/dds.rds",     emit: dds_rds
    path "results/vst.rds",     emit: vst_rds
    tuple val("${task.process}"), val('r-base'), eval('R --version | grep -Eo "[0-9]+(\\.[0-9]+)+"'),                                                     emit: v_r,        topic: versions
    tuple val("${task.process}"), val('tximport'), eval("Rscript -e 'cat(as.character(packageVersion(\"tximport\")))' | grep -Eo \"[0-9]+(\\.[0-9]+)+\""), emit: v_tximport, topic: versions
    tuple val("${task.process}"), val('DESeq2'), eval("Rscript -e 'cat(as.character(packageVersion(\"DESeq2\")))' | grep -Eo \"[0-9]+(\\.[0-9]+)+\""),     emit: v_deseq2,   topic: versions
    tuple val("${task.process}"), val('readr'), eval("Rscript -e 'cat(as.character(packageVersion(\"readr\")))' | grep -Eo \"[0-9]+(\\.[0-9]+)+\""),       emit: v_readr,    topic: versions
    tuple val("${task.process}"), val('dplyr'), eval("Rscript -e 'cat(as.character(packageVersion(\"dplyr\")))' | grep -Eo \"[0-9]+(\\.[0-9]+)+\""),       emit: v_dplyr,    topic: versions
    tuple val("${task.process}"), val('tibble'), eval("Rscript -e 'cat(as.character(packageVersion(\"tibble\")))' | grep -Eo \"[0-9]+(\\.[0-9]+)+\""),     emit: v_tibble,   topic: versions
    
    script:
    """
    deseq2_analysis.R \\
        inputs \\
        ${samplesheet} \\
        results \\
        ${reference_level}
    """
}