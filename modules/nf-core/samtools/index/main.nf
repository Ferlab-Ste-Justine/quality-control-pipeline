process SAMTOOLS_INDEX {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/samtools:1.24--h9dcdb79_1'
        : 'quay.io/biocontainers/samtools:1.24--h9dcdb79_1'}"

    input:
    tuple val(meta), path(input)

    output:
    tuple val(meta), path("*.{bai,csi,crai}"), emit: index
    tuple val("${task.process}"), val('samtools'), eval("samtools version | sed '1!d;s/.* //'"), emit: versions_samtools, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    samtools \\
        index \\
        -@ ${task.cpus-1} \\
        ${args} \\
        ${input}
    """

    stub:
    def args = task.ext.args ?: ''
    def extension = file(input).getExtension() == 'cram'
        ? "crai"
        : args.contains("-c") ? "csi" : "bai"
    """
    touch ${input}.${extension}
    """
}
