// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/modules/nf-core/
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// TODO nf-core: Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

process QC_ANALYSIS_DEPTH_OF_SAMPLES {
    tag 'depth_analysis'
    label 'process_single'

    // TODO nf-core: List required Conda package(s).
    //               Software MUST be pinned to channel (i.e. "bioconda"), version (i.e. "1.10").
    //               For Conda, the build (i.e. "h9402c20_2") must be EXCLUDED to support installation on different operating systems.
    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "docker.io/ferlabcrsj/qc-pipeline-python:1.1.0"

    input:
    path multiqc_data

    output:
    path "errors/erreurs_de_couvertures_par_region.txt", emit: errors_by_region, optional: true
    path "errors/erreurs_de_couvertures_mediane.txt", emit: errors_by_median, optional: true
    path "errors/erreurs_de_couvertures_anormales.txt", emit: errors_by_outliers, optional: true
    path "qc_pass", emit: qc_pass, optional: true
    path "qc_fail", emit: qc_fail, optional: true
    path "command_output_analyse_mosdepth/mosdepth_analysis_output.txt", emit: mosdepth_analysis
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def coverage_per_contig = "$multiqc_data/mosdepth_perchrom.txt"
    def multiqc_general_stats = "$multiqc_data/multiqc_general_stats.txt"
    
    // TODO nf-core: Where possible, a command MUST be provided to obtain the version number of the software e.g. 1.10
    //               If the software is unable to output a version number on the command-line then it can be manually specified
    //               e.g. https://github.com/nf-core/modules/blob/master/modules/nf-core/homer/annotatepeaks/main.nf
    //               Each software used MUST provide the software name and version number in the YAML version file (versions.yml)
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "task.ext.args" directive
    // TODO nf-core: If the tool supports multi-threading then you MUST provide the appropriate parameter
    //               using the Nextflow "task" variable e.g. "--threads $task.cpus"
    // TODO nf-core: Please replace the example samtools command below with your module's command
    // TODO nf-core: Please indent the command appropriately (4 spaces!!) to help with readability ;)
    """
    mkdir -p command_output_analyse_mosdepth
    analyse_mosdepth.py ${args} $coverage_per_contig $multiqc_general_stats > command_output_analyse_mosdepth/mosdepth_analysis_output.txt 2>&1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        depth_analysis_of_samples: \$(python --version |& sed '1!d ; s/python //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    
    // TODO nf-core: A stub section should mimic the execution of the original module as best as possible
    //               Have a look at the following examples:
    //               Simple example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bcftools/annotate/main.nf#L47-L63
    //               Complex example: https://github.com/nf-core/modules/blob/818474a292b4860ae8ff88e149fbcda68814114d/modules/nf-core/bedtools/split/main.nf#L38-L54
    """
    touch erreurs_de_couvertures_par_region.txt
    touch erreurs_de_couvertures_mediane.txt
    touch erreurs_de_couvertures_anormales.txt
    touch mosdepth_analysis_output.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        depth_analysis_of_samples: \$(python --version |& sed '1!d ; s/python //')
    END_VERSIONS
    """
}
