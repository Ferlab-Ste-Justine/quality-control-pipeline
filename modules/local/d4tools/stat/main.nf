process D4TOOLS_STAT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/d4tools:0.3.11--h3ab6199_2':
        'biocontainers/d4tools:0.3.11--h3ab6199_2' }"

    input:
    tuple val(meta), path(d4)
    tuple val(meta2), path(bed_regions)

    output:
    tuple val(meta), path("*.txt"), emit: stat
    path "versions.yml"           , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def regions_arg = bed_regions ? "--region ${bed_regions} " : ''
    """
    d4tools stat -H \\
        $args \\
        $regions_arg \\
        -t $task.cpus \\
        $d4 > ${prefix}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        d4tools: \$(d4tools --version)
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def regions_arg = bed_regions ? "--region ${bed_regions} " : ''
    """
    echo $args $regions_arg

    touch ${prefix}.d4stat.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        d4tools: \$(d4tools --version)
    END_VERSIONS
    """
}
