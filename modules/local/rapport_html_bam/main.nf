process RAPPORT_HTML_BAM {
    tag 'qc_analysis_report'
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "docker.io/ferlabcrsj/qc-pipeline-python:1.1.0"

    input:
    path mosdepth_txt
    path samtools_txt
    path plot_html
    path median_coverage_csv
    path region_coverage_csv
    path other_error_csv

    output:
    path "qc/qc_analysis_report.html", emit: report
    path "samtools_reads_info.png", emit: plot
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def median = median_coverage_csv ? median_coverage_csv : ""
    def region = region_coverage_csv ? region_coverage_csv : ""
    def other = other_error_csv ? other_error_csv : ""
    """
    html_generator_bam.py \\
        $mosdepth_txt \\
        $samtools_txt \\
        $plot_html \\
        qc \\
        ${params.outdir}/qc/qc_pass \\
        $median \\
        $region \\
        $other
    
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

