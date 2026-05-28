process MULTIQC_PYTHON {
    label 'process_single'
    tag "${meta.id}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.31--pyhdfd78af_0' :
        'biocontainers/multiqc:1.30--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(multiqc_files)
    path multiqc_config

    output:
    tuple val(meta), path("*multiqc_report.html"), emit: report
    tuple val(meta), path("*_data")              , emit: data
    tuple val(meta), path("*_plots")             , optional:true, emit: plots
    tuple val("${task.process}"), val('multiqc'), eval('multiqc --version | sed "s/.* //g"'), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def config = multiqc_config ? "--config $multiqc_config" : ''
    def file_list = multiqc_files.join(' ')
    """
    multiqc_report.py \\
        ${file_list} \\
        ${config} \\
        --filename ${prefix}_multiqc_report.html \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}_multiqc_data
    mkdir ${prefix}_multiqc_plots
    touch ${prefix}_multiqc_report.html
    """
}
