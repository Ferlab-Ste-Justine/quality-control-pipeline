process RAPPORT_HTML_FASTQ {
    tag 'qc_analysis_report'
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "docker.io/ferlabcrsj/qc-pipeline-python:1.1.0"

    input:
    path analyse_fastq
    path errors

    output:
    path "qc/qc_analysis_report.html", emit: report
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    html_generator_fastq.py \\
    $analyse_fastq \\
    $errors \\
    qc \\   
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qc_analysis_report: \$(python --version |& sed '1!d ; s/python //')
    END_VERSIONS
    """


    stub:
    def args = task.ext.args ?: ''
    """
    mkdir -p ${params.outdir}/qc
    touch ${params.outdir}/qc/qc_analysis_report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        qc_analysis_report: \$(python --version |& sed '1!d ; s/python //')
    END_VERSIONS
    """

}

