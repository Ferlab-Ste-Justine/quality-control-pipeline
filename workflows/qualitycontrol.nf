/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_quality-control-pipeline_pipeline'
include { FASTQ_QC               } from '../subworkflows/local/fastq_qc/main'
include { BAM_QC as BAM_QC_WGS   } from '../subworkflows/local/bam_qc/main'
include { BAM_QC as BAM_QC_TARGET   } from '../subworkflows/local/bam_qc/main'
include { VCF_QC                 } from '../subworkflows/local/vcf_qc/main'
include { BAM_MERGE              } from '../subworkflows/local/bam_merge'
include { CRAM_SOMALIER          } from '../subworkflows/local/cram_somalier'
include { GATK4_BEDTOINTERVALLIST } from '../modules/nf-core/gatk4/bedtointervallist/main'
include { CREATE_FAMILY_PED      } from '../modules/local/create_family_ped/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow QUALITYCONTROL {

    take:
    ch_samplesheet    // channel: [ val(meta), path(fastq/cram/bam/vcf files) ]

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // inputs
    ch_fasta = params.fasta ? channel.value(file(params.fasta, checkIfExists:true)) : channel.value([])
    ch_fai   = params.fai ? channel.value(file(params.fai, checkIfExists:true)) : channel.value([])
    ch_dict   = params.fasta_dict ? channel.value(file(params.fasta_dict, checkIfExists:true)) : channel.value([])
    ch_intervals = params.regions_bed ? channel.value(file(params.regions_bed, checkIfExists: true)) : channel.value([])
    qc_regions_1 = params.qc_coverage_region_1 ? channel.value(file(params.qc_coverage_region_1, checkIfExists:true)) : channel.value([])
    qc_regions_2 = params.qc_coverage_region_2 ? channel.value(file(params.qc_coverage_region_2, checkIfExists:true)) : channel.value([])
    ch_targets = params.targets_bed ? channel.value(file(params.targets_bed, checkIfExists:true)) : channel.value([])
    ch_exons = params.exons_bed ? channel.value(file(params.exons_bed, checkIfExists:true)) : channel.value([])

    // somalier sites VCF
    ch_somalier_sites = params.somalier_sites ? channel.value(file(params.somalier_sites, checkIfExists:true)) : channel.value([])
    ch_ped = params.ped_file ? channel.value(file(params.ped_file, checkIfExists:true)) : channel.of([])
    // verifybamid SVD files
    ch_svd_ud  = params.verifybamid_svd_prefix ? channel.value(file(params.verifybamid_svd_prefix + '.UD', checkIfExists:true)) : channel.value([])
    ch_svd_mu  = params.verifybamid_svd_prefix ? channel.value(file(params.verifybamid_svd_prefix + '.mu', checkIfExists:true)) : channel.value([])
    ch_svd_bed = params.verifybamid_svd_prefix ? channel.value(file(params.verifybamid_svd_prefix + '.bed', checkIfExists:true)) : channel.value([])
    ch_svd_in = ch_svd_ud.combine(ch_svd_mu).combine(ch_svd_bed).collect()
    ncm_snp_pt = params.ngscheckmate_snp_pt ? channel.value(file(params.ngscheckmate_snp_pt, checkIfExists:true)) : channel.value([])

    // Input files can be fastq, cram/bam, or vcf - separate different data types
    // Branch input based on file type
    ch_samplesheet_parsed = ch_samplesheet
        .branch { meta, files ->
        fastq: meta.fileType == "FASTQ"
        aln: meta.fileType in ["BAM", "CRAM"]
            [ meta - meta.subMap('lane','runId'), files[0], files[1] ]
        vcf: meta.fileType in ["VCF","GVCF"]
            [ meta - meta.subMap('lane','runId'), files[0], files[1] ]
        }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        FASTQ QUALITY CONTROL
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // If fastq files are provided, run FastQC, SeqFu, and ngsCheckMate
    FASTQ_QC ( ch_samplesheet_parsed.fastq, ncm_snp_pt.map { it -> [ [id:"snp_pt"], it] } )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_QC.out.reports.collect{it})
    ch_versions = ch_versions.mix(FASTQ_QC.out.versions.first())

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        BAM/CRAM QUALITY CONTROL
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    GATK4_BEDTOINTERVALLIST(ch_intervals.map { it -> [[id:"bed"], it] },
                            ch_dict.map { it -> [[id:"dict"], it] })

    ch_interval_list = GATK4_BEDTOINTERVALLIST.out.interval_list
                        .map { _meta, intervals -> intervals }
                        .collect()
                        .ifEmpty([])

    // Merge multiple runs
    //
    // ----- SAMTOOLS MERGE -----
    //
    bam_to_merge = ch_samplesheet_parsed.aln
        .map { meta, cram, crai ->
        [ groupKey(meta.subMap('id', 'participant', 'sample', 'sequencingType', 'status'), meta.n_lanes), cram, crai ]
    }
    .groupTuple()

    BAM_MERGE(bam_to_merge, ch_fasta, ch_fai)

    ch_versions = ch_versions.mix(BAM_MERGE.out.versions)

    //
    // ----- ALIGNMENT QC -----
    //
    // separate qc for targeted seq vs wgs
    ch_bam_qc = BAM_MERGE.out.bam_bai
        .map { groupKey, bam, bai ->
        [groupKey.target, bam, bai]
        }
        .branch { meta, bam, bai ->
            wgs: meta.sequencingType == 'WGS'
            targeted: meta.sequencingType != 'WGS'
        }

    BAM_QC_WGS(
        ch_bam_qc.wgs,
        ch_fasta,
        ch_fai,
        [],
        qc_regions_1,
        qc_regions_2,
        [],
        ch_svd_in
    )

    BAM_QC_TARGET(
        ch_bam_qc.targeted,
        ch_fasta,
        ch_fai,
        ch_intervals,
        qc_regions_1,
        qc_regions_2,
        ch_interval_list,
        ch_svd_in
    )

    ch_multiqc_files = ch_multiqc_files.mix(BAM_QC_WGS.out.reports)
    ch_multiqc_files = ch_multiqc_files.mix(BAM_QC_TARGET.out.reports)

    ch_versions = ch_versions.mix(BAM_QC_WGS.out.versions)
    ch_versions = ch_versions.mix(BAM_QC_TARGET.out.versions)

    //
    // ----- CRAM_SOMALIER -----
    //
    ch_cram_crai_somalier = BAM_MERGE.out.bam_bai
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

    ch_multiqc_files = ch_multiqc_files.mix(CRAM_SOMALIER.out.pairs_tsv.map { _meta, report -> report })
    ch_multiqc_files = ch_multiqc_files.mix(CRAM_SOMALIER.out.samples_tsv.map { _meta, report -> report })

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        VCF QUALITY CONTROL
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    VCF_QC (
        ch_samplesheet_parsed.vcf,
        ch_fasta,
        ch_intervals,
        ch_targets,
        ch_exons
    )

    ch_multiqc_files = ch_multiqc_files.mix(VCF_QC.out.vcf_stats.collect{_meta, report -> report})

    //
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

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
