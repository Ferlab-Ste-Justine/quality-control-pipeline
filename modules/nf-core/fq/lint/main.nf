process FQ_LINT {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fq:0.12.0--h9ee0642_0':
        'biocontainers/fq:0.12.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("*.fq_lint.txt"), emit: lint
    tuple val(meta), path("*.fq_lint.status.yml"), emit: status
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    set +e
    FQLINT_EXIT_CODE=0

    echo "Running fq lint on ${fastq}"
    fq lint \\
        $args \\
        $fastq > ${prefix}.fq_lint.txt
    
    FQ_LINT_EXIT_CODE=\$?
    if [ \$FQ_LINT_EXIT_CODE -ne 0 ]; then
        echo -e "fq lint returned a non-zero exit status. \$FQ_LINT_EXIT_CODE"
        if [ \$FQ_LINT_EXIT_CODE -eq 1 ]; then
            FQ_STATUS=FAIL
            echo "ERROR: fq lint failed for ${prefix}"
            echo \$FQ_STATUS > ${prefix}.fq_lint.status.yml
        else
            exit \$FQ_LINT_EXIT_CODE
        fi
    else
        FQ_STATUS=PASS
        echo \$FQ_STATUS > ${prefix}.fq_lint.status.yml 
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fq: \$(echo \$(fq lint --version | sed 's/fq-lint //g'))
    END_VERSIONS

    exit 0
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fq_lint.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fq: \$(echo \$(fq lint --version | sed 's/fq-lint //g'))
    END_VERSIONS
    """
}
