process MERGE_REGIONS_COV {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1':
        'biocontainers/python:3.9--1' }"

    input:
    val(meta), path(mean), path(median), path(pcov), path(count)
    val(region_bed)

    output:
    val(meta), path("*.txt") emit: report
    // tuple val("${task.process}"), val('regions_cov'), eval('regions_cov.py --version | sed "s/.* //g"'), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // def prefix = task.ext.prefix ?: "${meta.id}"
    def config = multiqc_config ? "--config $multiqc_config" : ''
    """
    echo -e "#Chr\tStart\tEnd\tmean\tmedian\tcount" > merged_regions_cov.tmp
    paste $mean <(cut -f4 $median) <(cut -f4 $pcov) <(cut -f4 $count) >> merged_regions_cov.tmp
    paste merged_regions_cov.tmp <(cut -f1-3 --complement  $pcov) > merged_regions_cov.txt
    """

    stub:
    """
    touch qc_regions.txt
    """
}
