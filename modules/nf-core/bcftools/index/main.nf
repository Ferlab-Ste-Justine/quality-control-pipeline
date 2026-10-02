process BCFTOOLS_INDEX {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/bcftools:1.23.1--hb2cee57_0'
        : 'quay.io/biocontainers/bcftools:1.23.1--hb2cee57_0'}"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.{tbi,csi}"), emit: index, optional: true
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | sed '1!d; s/^.*bcftools //'"), topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    bcftools \\
        index \\
        ${args} \\
        --threads ${task.cpus} \\
        ${vcf}
    """

    stub:
    def args = task.ext.args ?: ''
    def extension = args.contains("--tbi") || args.contains("-t")
        ? "tbi"
        : "csi"
    """
    touch ${vcf}.${extension}
    """
}
