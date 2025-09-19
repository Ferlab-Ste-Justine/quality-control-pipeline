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
include { SAMTOOLS_CONVERT       } from '../modules/nf-core/samtools/convert/main'
include { PICARD_VALIDATESAMFILE } from '../modules/local/picard/validatesamfile'
include { SAMTOOLS_SAMPLES       } from '../modules/local/samtools/samples'

include { VCF_ID_REPAIR          } from '../subworkflows/local/vcf_id_repair'
include { BAM_ID_REPAIR          } from '../subworkflows/local/bam_id_repair'
include { BAM_MERGE              } from '../subworkflows/local/bam_merge_reheader'
include { CRAM_SOMALIER          } from '../subworkflows/local/cram_somalier'
include { BAM_QC as BAM_QC_WGS   } from '../subworkflows/local/bam_qc'
include { BAM_QC as BAM_QC_TARGET   } from '../subworkflows/local/bam_qc'
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
    ch_dict   = params.fasta_dict ? Channel.value(file(params.fasta_dict, checkIfExists:true)) : Channel.value([])
    ch_intervals = params.regions_bed ? Channel.value(file(params.regions_bed, checkIfExists: true)) : Channel.value([])
    // somalier sites VCF
    ch_somalier_sites = params.somalier_sites ? Channel.value(file(params.somalier_sites, checkIfExists:true)) : Channel.value([])
    ch_ped = params.ped_file ? Channel.value(file(params.ped_file, checkIfExists:true)) : Channel.of([])
    // verifybamid SVD files
    ch_svd_ud  = params.verifybamid_svd_prefix ? Channel.value(file(params.verifybamid_svd_prefix + '.UD', checkIfExists:true)) : Channel.value([])
    ch_svd_mu  = params.verifybamid_svd_prefix ? Channel.value(file(params.verifybamid_svd_prefix + '.mu', checkIfExists:true)) : Channel.value([])
    ch_svd_bed = params.verifybamid_svd_prefix ? Channel.value(file(params.verifybamid_svd_prefix + '.bed', checkIfExists:true)) : Channel.value([])
    ch_svd_in = ch_svd_ud.mix(ch_svd_mu).mix(ch_svd_bed)

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
    ch_in_validate = ch_aln
        .map { meta, bam, bai ->
            def new_lane = meta.n_lanes > 1 ? meta.lane : ""
            def new_meta = meta + [id:meta.sample + (meta.n_lanes > 1 ? ".${meta.lane}" : "")] + [lane: new_lane]
            [new_meta, bam, bai]
        }

    PICARD_VALIDATESAMFILE( 
        ch_in_validate,
        ch_fasta.map { it -> [ [id:"fasta"], it] }, 
        ch_fai.map { it -> [ [id:"fai"], it] },
        ch_dict.map { it -> [ [id:"dict"], it] }
    )

    ch_versions = ch_versions.mix(PICARD_VALIDATESAMFILE.out.versions)

    // Joining with validate output to create dependency
    ch_cram_crai_validated = PICARD_VALIDATESAMFILE.out.txt
        .join(ch_in_validate)
        .map { meta, _txt, cram, crai ->
        [meta, cram, crai]
    }

    //
    // ----- Validate ID
    //
    BAM_ID_REPAIR(ch_cram_crai_validated, ch_fasta)
    ch_versions = ch_versions.mix(BAM_ID_REPAIR.out.versions)

    ch_cram_crai = BAM_ID_REPAIR.out.bam_bai 
        .map { meta, bam, bai ->
        def new_meta = meta + [id: meta.sample] // Reset ID to sample name
        [new_meta, bam, bai]
        }

    //
    // ----- CRAM_SOMALIER -----
    //
    ch_cram_crai_somalier = ch_cram_crai
        .map { meta, cram, crai ->
        [meta.sample, meta, cram, crai]
        }
        .groupTuple()
        .map { _sample, meta, cram, crai ->  [meta.withIndex(), cram, crai] }
        .transpose()
        .map { meta_idx, cram, crai -> [meta_idx[0] + [sample_idx:meta_idx[1]], cram, crai] }
        .map { meta, cram, crai ->
        def prefix = (meta.n_seqTypes > 1 ? "${meta.sequencingType}." : "") + (meta.n_lanes > 1 ? "${meta.lane}" : "")
        def samplename_somalier = meta.sample_idx > 0 ? "${meta.sample}.${prefix}" : meta.sample
        [meta + [prefix:prefix, samplename_somalier:samplename_somalier, id:samplename_somalier], cram, crai, meta.n_lanes]
    }

    // Participant group - group aln files of the same participant together. Will not create file if resulting filter is empty
    ch_cram_crai_somalier
        .map{ meta, _cram, _crai, _lanes ->
            [ meta.participant, meta.samplename_somalier]
        }
        .unique()
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

    // Merge multiple runs 
    //
    // ----- SAMTOOLS MERGE -----
    //
    bam_to_merge = ch_cram_crai
        .map { meta, cram, crai ->
        [groupKey(meta.subMap('id', 'participant', 'sample', 'sequencingType', 'status'), meta.n_lanes), cram, crai]
    }
    .groupTuple()

    BAM_MERGE(bam_to_merge, ch_fasta, ch_fai)
    ch_cram_crai_merged = BAM_MERGE.out.bam_bai
    ch_versions = ch_versions.mix(BAM_MERGE.out.versions)

    //
    // ----- ALIGNMENT QC -----
    //
    // separate qc for targeted seq vs wgs
    ch_bam_qc = ch_cram_crai_merged
        .branch { meta, bam, bai ->
            wgs: meta.sequencingType == 'WGS'
            targeted: meta.sequencingType != 'WGS'
        }

    BAM_QC_WGS(
        ch_bam_qc.wgs,
        ch_fasta, 
        ch_fai,
        [],
        ch_svd_in
    )

    BAM_QC_TARGET(
        ch_bam_qc.targeted,
        ch_fasta, 
        ch_fai,
        ch_intervals,
        ch_svd_in
    )

    ch_versions = ch_versions.mix(BAM_QC_WGS.out.versions)
    ch_versions = ch_versions.mix(BAM_QC_TARGET.out.versions)
    



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COLLECT SOFTWARE VERSIONS & MultiQC
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

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
    ch_multiqc_files = ch_multiqc_files.mix(BAM_QC_WGS.out.multiqc.map{it[1]}.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(BAM_QC_TARGET.out.multiqc.map{it[1]}.collect().ifEmpty([]))

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
