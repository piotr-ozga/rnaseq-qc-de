process VISUALIZE_DEA {
    tag "Visualizing ${annotated_tsv.baseName}"
    label "process_low"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-deseq2_r-base_r-dplyr_r-ggplot2_pruned:25d0ba752efaeac7' :
        'community.wave.seqera.io/library/bioconductor-deseq2_r-base_r-dplyr_r-ggplot2_pruned:fee1118c5f73317e' }"

    input:
    path annotated_tsv
    path vst_rds
    path samplesheet
    val ref_level

    output:
    path "plots/*.pdf",     emit: plots
    tuple val("${task.process}"), val('r-base'), eval('R --version | grep -Eo "[0-9]+(\\.[0-9]+)+"'),                                                          emit: v_r,        topic: versions
    tuple val("${task.process}"), val('ggplot2'), eval("Rscript -e 'cat(as.character(packageVersion(\"ggplot2\")))' | grep -Eo \"[0-9]+(\\.[0-9]+)+\""),       emit: v_ggplot2,  topic: versions
    tuple val("${task.process}"), val('dplyr'), eval("Rscript -e 'cat(as.character(packageVersion(\"dplyr\")))' | grep -Eo \"[0-9]+(\\.[0-9]+)+\""),           emit: v_dplyr,    topic: versions
    tuple val("${task.process}"), val('ggrepel'), eval("Rscript -e 'cat(as.character(packageVersion(\"ggrepel\")))' | grep -Eo \"[0-9]+(\\.[0-9]+)+\""),       emit: v_ggrepel,  topic: versions
    tuple val("${task.process}"), val('readr'), eval("Rscript -e 'cat(as.character(packageVersion(\"readr\")))' | grep -Eo \"[0-9]+(\\.[0-9]+)+\""),           emit: v_readr,    topic: versions
    tuple val("${task.process}"), val('pheatmap'), eval("Rscript -e 'cat(as.character(packageVersion(\"pheatmap\")))' | grep -Eo \"[0-9]+(\\.[0-9]+)+\""),     emit: v_pheatmap, topic: versions
    
    script:
    """
    visualize_dea.R \\
        ${annotated_tsv} \\
        ${vst_rds} \\
        ${samplesheet} \\
        ${ref_level} \\
        plots
    """
}