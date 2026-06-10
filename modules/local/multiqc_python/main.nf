process MULTIQC_PYTHON {
    label 'process_single'
    tag "${meta.id}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc:1.31--pyhdfd78af_0' :
        'biocontainers/multiqc:1.30--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(multiqc_files), path(ped_file)
    path multiqc_config
    path thresholds_file

    output:
    tuple val(meta), path("*multiqc_report.html"), emit: report
    tuple val(meta), path("*_data*")             , emit: data
    tuple val(meta), path("*_plots")             , optional:true, emit: plots
    tuple val(meta), path("qc_json/*.metrics.json"), optional:true, emit: qc_json
    tuple val("${task.process}"), val('multiqc'), eval('multiqc --version | sed "s/.* //g"'), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args         = task.ext.args ?: ''
    def prefix       = task.ext.prefix ?: "${meta.id}"
    def mode         = task.ext.mode ?: (meta.id == 'Cohort' ? 'cohort' : 'family')
    def config_arg   = multiqc_config ? "--config ${multiqc_config}" : ''
    def ped_arg      = ped_file       ? "--ped ${ped_file}"          : ''
    def thresh_arg   = thresholds_file? "--thresholds ${thresholds_file}" : ''
    def file_list    = multiqc_files.join(' ')
    """
    multiqc_report.py \\
        ${file_list} \\
        ${config_arg} \\
        ${ped_arg} \\
        ${thresh_arg} \\
        --mode ${mode} \\
        --filename ${prefix}_multiqc_report.html \\
        --json-out qc_json \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir ${prefix}_multiqc_data qc_json
    touch ${prefix}_multiqc_report.html qc_json/${prefix}.metrics.json
    """
}
