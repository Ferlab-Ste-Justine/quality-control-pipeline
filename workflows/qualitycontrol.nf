/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_quality-control-pipeline_pipeline'

include { FASTQC                 } from '../modules/nf-core/fastqc'
include { MULTIQC                } from '../modules/nf-core/multiqc'
include { CREATE_FAMILY_PED      } from '../modules/local/create_family_ped'
include { PICARD_VALIDATESAMFILE } from '../modules/local/picard/validatesamfile'
include { SAMTOOLS_SAMPLES       } from '../modules/local/samtools/samples'
include { QUALIMAP_BAMQCCRAM     } from '../modules/nf-core/qualimap/bamqccram'
include { VERIFYBAMID_VERIFYBAMID2 } from '../modules/nf-core/verifybamid/verifybamid2'

include { VCF_ID_REPAIR          } from '../subworkflows/local/vcf_id_repair'
include { BAM_ID_REPAIR          } from '../subworkflows/local/bam_id_repair'
include { CRAM_SOMALIER          } from '../subworkflows/local/cram_somalier'
include { FASTQ_QC                 } from '../subworkflows/local/fastq_qc/main'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow QUALITYCONTROL {

    take:
    ch_fastq
    ch_aln
    ch_gvcf // channel: samplesheet read in from --input

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_final_reports = Channel.empty()

    // inputs
    ch_fasta = Channel.value(file(params.fasta, checkIfExists:true)) 
    ch_fai   = Channel.value(file(params.fai, checkIfExists:true))
    ch_dict   = Channel.value(file(params.fasta_dict, checkIfExists:true))
    // somalier sites VCF
    ch_somalier_sites = params.somalier_sites ? Channel.value(file(params.somalier_sites, checkIfExists:true)) : Channel.of([])
    ch_ped = params.ped_file ? Channel.value(file(params.ped_file, checkIfExists:true)) : Channel.of([])
    // verifybamid SVD files
    ch_svd_ud  = params.verifybamid_svd_prefix ? Channel.value(file(params.verifybamid_svd_prefix + '.UD', checkIfExists:true)) : Channel.of([])
    ch_svd_mu  = params.verifybamid_svd_prefix ? Channel.value(file(params.verifybamid_svd_prefix + '.mu', checkIfExists:true)) : Channel.of([])
    ch_svd_bed  = params.verifybamid_svd_prefix ? Channel.value(file(params.verifybamid_svd_prefix + '.bed', checkIfExists:true)) : Channel.of([])


    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        FASTQ QUALITY CONTROL
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    FASTQ_QC ( ch_fastq )
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_QC.out.reports.collect{it})
        ch_versions = ch_versions.mix(FASTQ_QC.out.versions)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        BAM/CRAM QC
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    //
    // ----- Picard validation ------ 
    //
    PICARD_VALIDATESAMFILE(
        ch_aln,
        ch_fasta.map { it -> [ [id:"fasta"], it] }, 
        ch_fai.map { it -> [ [id:"fai"], it] },
        ch_dict.map { it -> [ [id:"dict"], it] }
    )

    ch_versions = ch_versions.mix(PICARD_VALIDATESAMFILE.out.versions)

    // Joining with validate output to create dependency
    ch_cram_crai = PICARD_VALIDATESAMFILE.out.txt
        .join(ch_aln)
        .map { meta, _txt, cram, crai ->
        [meta, cram, crai]
    }

    //
    // ----- Validate ID
    //
    BAM_ID_REPAIR(ch_cram_crai, [])
    ch_versions = ch_versions.mix(BAM_ID_REPAIR.out.versions)

    ch_cram_crai = BAM_ID_REPAIR.out.bam_bai

    //
    // ----- CRAM_SOMALIER -----
    //
    ch_cram_crai_somalier = ch_cram_crai
        .map { meta, cram, crai ->
        def subset_meta = meta.subMap(['id','familyId','sample'])
        [subset_meta, cram, crai, meta.sampleSize]
    }

    // Participant group - group aln files of the same participant together. Will not create file if resulting filter is empty
    ch_aln
        .map{ meta, _cram, _crai ->
            [ meta.participant, meta.sample]
        }
        .groupTuple()
        .filter{ participant, ch_samples -> ch_samples.size() > 1 } 
        .collectFile(name: 'sample_groups.txt') { 
            participant, samples ->
            "${samples.join(',')}\n"
        }
        .ifEmpty{[]}
        .set{ch_sample_groups}

    ch_somalier_input_ped =  ch_ped.map { it -> [ [id:"ped"], it ] }

    // If we want to run the analysis in a per-family basis or with the entire cohort
    if (params.somalier_perfamily) {
        ch_ped_grouped = ch_ped
            .splitCsv(sep: '\t', header: ["family_id","name","paternal_id", "maternal_id", "sex", "phenotype"], skip: 1)
            .map { row -> 
                [['familyId':row.family_id], row] }
            .groupTuple()

        CREATE_FAMILY_PED(ch_ped_grouped)

        ch_versions = ch_versions.mix(CREATE_FAMILY_PED.out.versions)

        ch_somalier_input_ped = CREATE_FAMILY_PED.out.ped_files
            .map { meta, ped_file -> 
                [ [id: meta.familyId] + meta, ped_file ]
            }
    }
    
    CRAM_SOMALIER(
            ch_cram_crai_somalier,
            ch_fasta.map { it -> [ [id:"fasta"], it] },
            ch_fai.map { it -> [ [id:"fai"], it] },
            ch_somalier_sites.map { it -> [ [id:"sites"], it] },
            ch_somalier_input_ped ?: ch_ped.map { it -> [ [id:"ped"], it ] },
            ch_sample_groups,
            'familyId'  // Common identifier to relate samples by (family ID in this case)
        )

    ch_versions = ch_versions.mix(CRAM_SOMALIER.out.versions)
    ch_final_reports = ch_final_reports.mix(CRAM_SOMALIER.out.html)

    // Validate header and merge multiple runs 
    // TODO: Add subworkflow for reheader and merging

    //
    // ----- CRAM QUALIMAP -----
    //
    QUALIMAP_BAMQCCRAM (
        ch_cram_crai, 
        params.regions_bed ? Channel.value(file( params.regions_bed, checkIfExists: true )).toList() : [],
        ch_fasta, ch_fai)

    ch_qualimap = QUALIMAP_BAMQCCRAM.out.results
    ch_versions = ch_versions.mix(QUALIMAP_BAMQCCRAM.out.versions)

    //
    // ----- VERIFYBAMID2 - Contamination -----
    //
    ch_svd_in = ch_svd_ud.combine(ch_svd_mu).combine(ch_svd_bed).collect()
    VERIFYBAMID_VERIFYBAMID2(
        ch_cram_crai, ch_svd_in, [], ch_fasta)

    ch_self_sm = VERIFYBAMID_VERIFYBAMID2.out.self_sm
    
    // //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'quality-control-pipeline_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    ch_multiqc_files = ch_multiqc_files.mix(PICARD_VALIDATESAMFILE.out.txt.map{it[1]}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(CRAM_SOMALIER.out.pairs_tsv.map { _meta, report -> report })
    ch_multiqc_files = ch_multiqc_files.mix(CRAM_SOMALIER.out.samples_tsv.map { _meta, report -> report })

    ch_multiqc_files = ch_multiqc_files.mix(QUALIMAP_BAMQCCRAM.out.results.map{it[1]}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(VERIFYBAMID_VERIFYBAMID2.out.self_sm.map{it[1]}.collect().ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    final_reports  = ch_final_reports
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
