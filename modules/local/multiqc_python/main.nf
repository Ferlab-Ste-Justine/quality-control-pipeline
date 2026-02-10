process MULTIQC_PYTHON {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.31--pyhdfd78af_0' :
        'biocontainers/multiqc:1.30--pyhdfd78af_0' }"

    input:
    path multiqc_files
    path multiqc_config

    output:
    path "*multiqc_report.html", emit: report
    path "*_data"              , emit: data
    path "*_plots"             , optional:true, emit: plots
    tuple val("${task.process}"), val('multiqc'), eval('multiqc --version | sed "s/.* //g"'), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // def prefix = task.ext.prefix ?: "${meta.id}"
    def config = multiqc_config ? "--config $multiqc_config" : ''
    def file_list = multiqc_files.join(' ')
    """
    multiqc_report.py \\
        ${file_list} \\
        ${config} \\
        $args
    """

    stub:
    """
    mkdir multiqc_data
    mkdir multiqc_plots
    touch multiqc_report.html
    """
}
