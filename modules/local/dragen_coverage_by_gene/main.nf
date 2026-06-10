process DRAGEN_COVERAGE_BY_GENE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1':
        'biocontainers/pandas:2.2.1' }"

    input:
    tuple val(meta), path(cov_report), path(read_cov_report)

    output:
    tuple val(meta), path("*.coverage_by_gene.tsv"), emit: report, optional: true
    path "versions.yml"                            , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    dragen_coverage_by_gene.py \\
        --cov-report $cov_report \\
        --read-cov-report $read_cov_report \\
        --output ${prefix}.coverage_by_gene.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen_coverage_by_gene.py: \$(dragen_coverage_by_gene.py --version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.coverage_by_gene.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dragen_coverage_by_gene.py: \$(dragen_coverage_by_gene.py --version)
    END_VERSIONS
    """
}
