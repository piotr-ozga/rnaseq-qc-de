process ANNOTATE_RESULTS {
    tag "Annotating ${results_tsv.baseName}"
    label "process_low"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bioconductor-rtracklayer_r-base_r-dplyr_r-readr:d171b9305a0d1976' :
        'community.wave.seqera.io/library/bioconductor-rtracklayer_r-base_r-dplyr_r-readr:0aa893da74c8680d' }"

    input:
    path results_tsv
    path gtf

    output:
    path "results/results_annotated.tsv", emit: annotated_tsv
    path "results/annotation_summary.txt",     emit: summary_txt
    tuple val("${task.process}"), val('r-base'), eval('R --version | grep -Eo "[0-9]+(\\.[0-9]+)+"'),                                                                emit: v_r,           topic: versions
    tuple val("${task.process}"), val('rtracklayer'), eval("Rscript -e 'cat(as.character(packageVersion(\"rtracklayer\")))' | grep -Eo \"[0-9]+(\\.[0-9]+)+\""),     emit: v_rtracklayer, topic: versions
    tuple val("${task.process}"), val('readr'), eval("Rscript -e 'cat(as.character(packageVersion(\"readr\")))' | grep -Eo \"[0-9]+(\\.[0-9]+)+\""),                 emit: v_readr,       topic: versions
    tuple val("${task.process}"), val('dplyr'), eval("Rscript -e 'cat(as.character(packageVersion(\"dplyr\")))' | grep -Eo \"[0-9]+(\\.[0-9]+)+\""),                 emit: v_dplyr,       topic: versions
    
    script:
    """
    annotate_results.R \\
        ${results_tsv} \\
        ${gtf} \\
        results
    """
}