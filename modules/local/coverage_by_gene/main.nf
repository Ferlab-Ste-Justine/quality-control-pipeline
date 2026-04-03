process COVERAGE_BY_GENE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pandas:2.2.1':
        'biocontainers/pandas:2.2.1' }"

    input:
    tuple val(meta), path(mean), path(thresholds)

    output:
    tuple val(meta), path("*.tsv"), emit: report
    path "versions.yml"           , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    aggregate_mosdepth_by_gene.py \\
        --mean $mean \\
        --thresholds $thresholds \\
        --output ${prefix}.coverage_by_gene.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aggregate_mosdepth_by_gene.py: \$(aggregate_mosdepth_by_gene.py --version)
        python: \$(python3 --version | sed 's/Python //g')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch qc_regions.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aggregate_mosdepth_by_gene.py: \$(aggregate_mosdepth_by_gene.py --version)
        python: \$(python3 --version | sed 's/Python //g')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
