/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'

include { CRAM_QC_MOSDEPTH_SAMTOOLS as CRAM_QC_NO_MD  } from '../subworkflows/local/cram_qc_mosdepth_samtools/main'
include { QC_ANALYSIS_FASTQC_OF_SAMPLES } from '../modules/local/qc_analysis_fastqc_of_samples'
include { QC_ANALYSIS_DEPTH_OF_SAMPLES } from '../modules/local/qc_analysis_depth_of_samples'
include { QC_ANALYSIS_SAMTOOLS_OF_SAMPLES } from '../modules/local/qc_analysis_samtools_of_samples'
include { RAPPORT_HTML_BAM } from '../modules/local/rapport_html_bam'
include { RAPPORT_HTML_FASTQ } from '../modules/local/rapport_html_fastq'

// QC for fastq files
include { FASTQC } from '../modules/nf-core/fastqc/main'

// Create multiqc report
include { MULTIQC } from '../modules/nf-core/multiqc/main'                                                                                                                           

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow QUALITY_CONTROL {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    fasta
    intervals_for_preprocessing

    main:
    // To gather all QC reports for MultiQC
    ch_multiqc_files = Channel.empty()
    multiqc_report   = Channel.empty()
    reports          = Channel.empty()
    versions         = Channel.empty()

    // Check if input contains FASTQ files
    def inputPath = params.input
    def file = new File(inputPath)

    // Read the first line of the file
    def firstLine = file.withReader { reader ->
        reader.readLine()
    }

    // Check if the first line contains "fastq"
    contains_fastq = firstLine?.contains('fastq')


    if (contains_fastq) {

        FASTQC(ch_samplesheet)

        reports = reports.mix(FASTQC.out.zip.collect{ meta, logs -> logs })
        versions = versions.mix(FASTQC.out.versions.first())
    
    } else {

        CRAM_QC_NO_MD(ch_samplesheet, fasta, intervals_for_preprocessing)

        // Gather QC reports
        reports = reports.mix(CRAM_QC_NO_MD.out.reports.collect{ meta, report -> [ report ] })

        // Gather used softwares versions
        versions = versions.mix(CRAM_QC_NO_MD.out.versions)

    }

    //
    // Collate and save software versions
    //
    version_yaml = Channel.empty()
    version_yaml = softwareVersionsToYAML(versions)
            .collectFile(storeDir: "${params.outdir}/pipeline_info", name: 'nf_core_quality_control_software_mqc_versions.yml', sort: true, newLine: true)

    println "Project : $workflow.projectDir"

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config                     = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
    ch_multiqc_logo                       = params.multiqc_logo ? Channel.fromPath(params.multiqc_logo, checkIfExists: true) : Channel.empty()
    summary_params                        = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary                   = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files                      = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files                      = ch_multiqc_files.mix(version_yaml)
    ch_multiqc_files                      = ch_multiqc_files.mix(reports)

    // ch_multiqc_files.collect().view()

    MULTIQC (
            ch_multiqc_files.collect(),
            ch_multiqc_config.toList(),
            ch_multiqc_custom_config.toList(),
            ch_multiqc_logo.toList()
        )

    ch_multiqc_data = MULTIQC.out.data

    if (contains_fastq){

        QC_ANALYSIS_FASTQC_OF_SAMPLES(ch_multiqc_data)

        RAPPORT_HTML_FASTQ(
            QC_ANALYSIS_FASTQC_OF_SAMPLES.out.fastqc_analysis,
            QC_ANALYSIS_FASTQC_OF_SAMPLES.out.unequal_total_sequences,
            )

    } else {

        QC_ANALYSIS_DEPTH_OF_SAMPLES(ch_multiqc_data)

        QC_ANALYSIS_SAMTOOLS_OF_SAMPLES(ch_multiqc_data)

        ch_mosdepth_analysis = QC_ANALYSIS_DEPTH_OF_SAMPLES.out.mosdepth_analysis
        ch_samtools_analysis = QC_ANALYSIS_SAMTOOLS_OF_SAMPLES.out.samtools_analysis
        ch_plot = QC_ANALYSIS_SAMTOOLS_OF_SAMPLES.out.reads_info
        ch_median_coverage = QC_ANALYSIS_DEPTH_OF_SAMPLES.out.errors_by_median
        ch_errors_by_region = QC_ANALYSIS_DEPTH_OF_SAMPLES.out.errors_by_region
        ch_errors_by_outliers = QC_ANALYSIS_DEPTH_OF_SAMPLES.out.errors_by_outliers

        RAPPORT_HTML_BAM(
            ch_mosdepth_analysis,
            ch_samtools_analysis,
            ch_plot,
            ch_median_coverage,
            ch_errors_by_region,
            ch_errors_by_outliers
            )
    
    }

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
