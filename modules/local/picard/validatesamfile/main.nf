process PICARD_VALIDATESAMFILE {
    tag "$meta.id"
    label 'process_low'
    errorStrategy 'ignore'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/picard:3.4.0--hdfd78af_0':
        'biocontainers/picard:3.4.0--hdfd78af_0' }" //picard:3.3.0--hdfd78af_0

    input:
    tuple val(meta), path(bam), path(bai)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)
    tuple val(meta4), path(dict)

    output:
    tuple val(meta), path("*.txt"), emit: txt
    // tuple val(meta), path(".exitcode"), optional:true,emit: exit
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference_file = fasta ? "--REFERENCE_SEQUENCE ${fasta}" : ""

    def avail_mem = 3072
    if (!task.memory) {
        log.info '[Picard ValidateSamFile] Available memory not known - defaulting to 3GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.mega*0.8).intValue()
    }

    """
    picard \\
        -Xmx${avail_mem}M \\
        ValidateSamFile \\
        ${args} \\
        ${reference_file} \\
        --INPUT ${bam} \\
        --OUTPUT ${prefix}.txt \\
        || true

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard ValidateSamFile --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.validateSamFile.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        picard: \$(picard ValidateSamFile --version 2>&1 | grep -o 'Version:.*' | cut -f2- -d:)
    END_VERSIONS
    """
}
