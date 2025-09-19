process CREATE_FAMILY_PED {
    tag "$meta.familyId"
    label 'process_single'

    conda "conda-forge::python=3.13"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.13' :
        'biocontainers/python:3.13' }"

    input:
    tuple val(meta), val(ped_values)

    output:
    tuple val(meta), path("*.ped")  , emit: ped_files
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def familyId = meta.familyId
    def header = ['#family_id', 'name', 'paternal_id', 'maternal_id', 'sex', 'phenotype'].join('\t')
    def lines = ped_values.collect { sample ->
        ['family_id','name', 'paternal_id', 'maternal_id', 'sex', 'phenotype'].collect { sample[it] }.join('\t')
    } 
    def outfile_text = ( [header] + lines ).join('\n')
    """
    echo -e "$outfile_text" > ${familyId}.ped
  
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        create_family_ped: v1.0
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """

    stub:
    def familyId = meta.familyId
    """
    touch ${familyId}.ped

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        create_family_ped: v1.0
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
