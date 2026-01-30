process SAMTOOLS_REHEADER {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.21--h96c455f_1':
        'biocontainers/samtools:1.21--h96c455f_1' }"

    input:
    tuple val(meta), path(input) // channel: [ val(meta), path(bam/cram) ]

    output:
    tuple val(meta), path("*.{bam,cram}"), emit: bam
    path "versions.yml"           , topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def file_type = input.getExtension()
    if ("$input" == "${prefix}.${file_type}") error "Input and output names are the same, set prefix in module configuration to disambiguate!"
    """
    samtools \\
        reheader \\
        $args \\
        $input > ${prefix}.${file_type}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def file_type = input.getExtension()
    """
    touch ${prefix}.${file_type}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version |& sed '1!d ; s/samtools //')
    END_VERSIONS
    """
}
