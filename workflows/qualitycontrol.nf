/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC_PYTHON                } from '../modules/local/multiqc_python/main'
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_quality-control-pipeline_pipeline'
include { buildPedRowsForFamily  } from '../subworkflows/local/utils_nfcore_quality-control-pipeline_pipeline'
include { FASTQ_QC               } from '../subworkflows/local/fastq_qc/main'
include { BAM_QC as BAM_QC_WGS   } from '../subworkflows/local/bam_qc/main'
include { BAM_QC as BAM_QC_TARGET   } from '../subworkflows/local/bam_qc/main'
include { VCF_QC                 } from '../subworkflows/local/vcf_qc/main'
include { BAM_MERGE              } from '../subworkflows/local/bam_merge'
include { CRAM_SOMALIER          } from '../subworkflows/local/cram_somalier'
include { GATK4_BEDTOINTERVALLIST } from '../modules/nf-core/gatk4/bedtointervallist/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
workflow QUALITYCONTROL {

    take:
    ch_samplesheet    // channel: [ val(meta), path(fastq/cram/bam/vcf files) ]

    main:

    ch_versions = channel.empty()
    // ch_multiqc_files emits [val(meta), path(file)] tuples. meta.familyId tags the
    // family the file belongs to; cohort-wide files (no familyId) are routed to every
    // family report in per-family mode, or merged into the single report in cohort mode.
    ch_multiqc_files = channel.empty()

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
    ch_ped = params.ped_file ? channel.value(file(params.ped_file, checkIfExists:true)) : channel.value([])
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

    if (!params.dragen_metrics_dir) {
        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            FASTQ QUALITY CONTROL
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */

        // If fastq files are provided, run FastQC, SeqFu, and ngsCheckMate
        FASTQ_QC ( ch_samplesheet_parsed.fastq, ncm_snp_pt.map { it -> [ [id:"snp_pt"], it] } )
        ch_multiqc_files = ch_multiqc_files.mix(FASTQ_QC.out.reports)
        ch_versions = ch_versions.mix(FASTQ_QC.out.versions.first())

        /*
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            BAM/CRAM QUALITY CONTROL
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        */
        GATK4_BEDTOINTERVALLIST( ch_intervals.flatten().map { it -> [[id:"bed"], it] },
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
            [ groupKey(meta.subMap('id', 'participant', 'sample', 'familyId', 'sex', 'sequencingType', 'status', 'relationship_to_proband', 'affected_status'), meta.n_lanes), cram, crai ]
        }
        .groupTuple()

        BAM_MERGE(bam_to_merge, ch_fasta, ch_fai)

        ch_versions = ch_versions.mix(BAM_MERGE.out.versions)
        // Post-merge per-sample CRAM/CRAI feeds somalier in normal mode.
        ch_cram_crai_source = BAM_MERGE.out.bam_bai

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

        // HsMetrics needs bait + target intervals. Resolution order:
        //   bait  : params.targets_bed > params.regions_bed > [] (no HsMetrics input)
        //   target: params.exons_bed   > bait
        def hs_bait_path   = params.targets_bed ?: params.regions_bed
        def hs_target_path = params.exons_bed   ?: hs_bait_path
        ch_hs_bait   = hs_bait_path   ? channel.value(file(hs_bait_path,   checkIfExists:true)) : channel.value([])
        ch_hs_target = hs_target_path ? channel.value(file(hs_target_path, checkIfExists:true)) : channel.value([])


        BAM_QC_WGS(
            ch_bam_qc.wgs,
            ch_fasta,
            ch_fai,
            ch_dict,
            ch_intervals,
            qc_regions_1,
            qc_regions_2,
            ch_interval_list,
            ch_hs_bait,
            ch_hs_target,
            ch_svd_in,
            'WGS'
        )

        BAM_QC_TARGET(
            ch_bam_qc.targeted,
            ch_fasta,
            ch_fai,
            ch_dict,
            ch_intervals,
            qc_regions_1,
            qc_regions_2,
            ch_interval_list,
            ch_hs_bait,
            ch_hs_target,
            ch_svd_in,
            'TARGETED'
        )

        ch_multiqc_files = ch_multiqc_files.mix(BAM_QC_WGS.out.reports)
        ch_multiqc_files = ch_multiqc_files.mix(BAM_QC_TARGET.out.reports)

        ch_versions = ch_versions.mix(BAM_QC_WGS.out.versions)
        ch_versions = ch_versions.mix(BAM_QC_TARGET.out.versions)

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

        ch_multiqc_files = ch_multiqc_files.mix(VCF_QC.out.vcf_metrics)
        ch_multiqc_files = ch_multiqc_files.mix(VCF_QC.out.vcf_stats)

    } else {
        // DRAGEN mode skips BAM_MERGE; feed somalier the BAM/CRAM straight from
        // the samplesheet. Trim meta to the same fields BAM_MERGE keeps on its
        // GroupKey target so meta.n_lanes is absent — otherwise the per-family
        // groupTuple inside CRAM_SOMALIER emits each sample early (count=1) and
        // fails the per-family join.
        ch_cram_crai_source = ch_samplesheet_parsed.aln
            .map { meta, cram, crai ->
                [ meta.subMap('id', 'participant', 'sample', 'familyId', 'sex', 'sequencingType', 'status', 'relationship_to_proband', 'affected_status'), cram, crai ]
            }
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        PEDIGREE ANALYSIS - SOMALIER
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    ch_cram_crai_somalier = ch_cram_crai_source
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
        .first()  // value channel so the single sample_groups file is reused across every per-family relate
        .set{ch_sample_groups}

    //
    // Build pedigree input for somalier.
    //   1. If params.ped_file is set, use it (split per-family when !cohort_mode).
    //   2. Otherwise derive from samplesheet meta (relationship_to_proband, affected_status, sex).
    def pedHeader = ['#family_id','name','paternal_id','maternal_id','sex','phenotype'].join('\t')

    if (params.ped_file) {
        if (params.cohort_mode) {
            ch_somalier_input_ped = ch_ped.map { it -> [ [id:"ped"], it ] }
        } else {
            ch_somalier_input_ped = ch_ped
                .splitCsv(sep: '\t', header: ["family_id","name","paternal_id","maternal_id","sex","phenotype"], skip: 1)
                .map { row ->
                    def line = [row.family_id, row.name, row.paternal_id, row.maternal_id, row.sex, row.phenotype].join('\t')
                    [ "${row.family_id}.ped".toString(), line ]
                }
                .collectFile(
                    newLine: true,
                    sort: true,
                    seed: pedHeader
                )
                .map { ped_file -> [ [id: ped_file.baseName], ped_file ] }
        }
    } else {
        ch_somalier_input_ped = ch_samplesheet
            .map { meta, _files ->
                [ meta.sample, meta + [samplename_somalier: meta.sample, sample_idx: 0] ]
            }
            .unique { entry -> entry[0] }
            .map { _sample, meta -> [meta.familyId, meta] }
            .groupTuple()
            .flatMap { familyId, metas ->
                buildPedRowsForFamily(familyId, metas, !(params.cohort_mode as boolean))
            }
            .collectFile(
                storeDir: "${params.outdir}/reports/pedigree",
                newLine: true,
                sort: true,
                seed: pedHeader
            )
            .map { ped_file -> [ [id: ped_file.baseName], ped_file ] }
    }

    ch_somalier_input_ped_for_relate = ch_somalier_input_ped
    if (!params.cohort_mode) {
        ch_families_with_cram = ch_cram_crai_somalier
            .map { meta, _cram, _crai, _count -> meta.familyId }
            .unique()
        ch_somalier_input_ped_for_relate = ch_somalier_input_ped
            .map { meta, ped -> [ meta.id, meta, ped ] }
            .join(ch_families_with_cram.map { fid -> [ fid, true ] })
            .map { _id, meta, ped, _flag -> [ meta, ped ] }
    }

    CRAM_SOMALIER(
            ch_cram_crai_somalier,
            ch_fasta.map { it -> [ [id:"fasta"], it] },
            ch_fai.map { it -> [ [id:"fai"], it] },
            ch_somalier_sites.map { it -> [ [id:"sites"], it] },
            ch_somalier_input_ped_for_relate ?: ch_ped.map { it -> [ [id:"ped"], it ] },
            ch_sample_groups,
            'familyId'  // Common identifier to relate samples by (family ID in this case)
        )

    ch_versions = ch_versions.mix(CRAM_SOMALIER.out.versions)

    // Somalier outputs are already grouped: meta.id is the familyId in per-family mode,
    // or 'Cohort' in cohort mode. Map id -> familyId so the downstream grouping works.
    ch_multiqc_files = ch_multiqc_files.mix(
        CRAM_SOMALIER.out.pairs_tsv.map { meta, report ->
            [ params.cohort_mode ? meta : meta + [familyId: meta.id], report ]
        }
    )
    ch_multiqc_files = ch_multiqc_files.mix(
        CRAM_SOMALIER.out.samples_tsv.map { meta, report ->
            [ params.cohort_mode ? meta : meta + [familyId: meta.id], report ]
        }
    )

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        DRAGEN METRICS (alternative to BAM_QC / VCF_QC)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    if (params.dragen_metrics_dir) {
        // Build a sample -> familyId map from the samplesheet so DRAGEN files
        // pick up the right familyId for per-family report routing.
        ch_sample_family = ch_samplesheet
            .map { meta, _files -> [ meta.sample, meta.familyId ] }
            .unique()

        ch_dragen_files = channel.fromPath([
                "${params.dragen_metrics_dir}/*.csv",
                "${params.dragen_metrics_dir}/**/*.csv",
            ], checkIfExists: false)
            .map { f -> [ f.name.tokenize('.')[0], f ] }   // [sample, file]
            .combine(ch_sample_family, by: 0)              // [sample, file, familyId]
            .map { sample, f, familyId -> [ [id: sample, sample: sample, familyId: familyId], f ] }

        ch_multiqc_files = ch_multiqc_files.mix(ch_dragen_files)
    }

    // topic_versions = channel.topic('versions')
    topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }
    ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name:  'quality-control-pipeline_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

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
    ch_workflow_summary = channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    // Cohort-wide multiqc inputs (no familyId). These are wrapped as [meta, file]
    // with an empty meta so they merge cleanly with the per-sample/per-family files.
    ch_multiqc_cohort_files = channel.empty()
        .mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
        .mix(ch_collated_versions)
        .mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
        .map { f -> [ [:], f ] }

    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_cohort_files)

    if (params.cohort_mode) {
        // Single cohort-wide MULTIQC report.
        ch_multiqc_input = ch_multiqc_files
            .map { _meta, f -> f }
            .collect()
            .map { fs -> [ [id: 'Cohort'], fs.flatten() ] }
    }
    else {
        // One MULTIQC report per family. Cohort-wide files (no familyId) are
        // attached to every family report so each is self-contained.
        ch_multiqc_branched = ch_multiqc_files.branch { meta, mqc_file ->
            per_family: meta.familyId
                [ meta.familyId, mqc_file ]
            cohort: true
                mqc_file
        }

        // Wrap cohort_files in an extra list so .combine() treats it as a single
        // positional value (combine unpacks tuples/lists into positional args otherwise).
        ch_cohort_collected = ch_multiqc_branched.cohort.collect().map { fs -> [fs] }
        ch_multiqc_input = ch_multiqc_branched.per_family
            .groupTuple()
            .combine(ch_cohort_collected)
            .map { familyId, family_files, cohort_files ->
                // flatten: some emits (e.g. paired-end fastqc zips) are lists of paths.
                [ [id: familyId], (family_files + cohort_files).flatten() ]
            }
    }

    // Join PED to the multiqc input by id. In cohort mode the ped meta.id may be
    // 'ped' (from params.ped_file) or 'Cohort' (samplesheet-derived); normalize.
    ch_ped_by_id = ch_somalier_input_ped.map { ped_meta, ped ->
        def id = params.cohort_mode ? 'Cohort' : ped_meta.id
        [ id, ped ]
    }
    ch_multiqc_input_with_ped = ch_multiqc_input
        .map { meta, files -> [ meta.id, meta, files ] }
        .join(ch_ped_by_id, remainder: true)
        .map { _id, meta, files, ped -> [ meta, files, ped ?: [] ] }

    // Use the user-supplied thresholds when params.qc_thresholds is set; otherwise
    // fall back to the bundled defaults in assets/qc_thresholds.yml.
    ch_qc_thresholds = channel.value(file(
        params.qc_thresholds ?: "$projectDir/assets/qc_thresholds.yml",
        checkIfExists: true,
    ))

    MULTIQC_PYTHON (
        ch_multiqc_input_with_ped,
        ch_multiqc_config.toList(),
        ch_qc_thresholds
    )
    //,
    //     ch_multiqc_custom_config.toList(),
    //     ch_multiqc_logo.toList(),
    //     [],
    //     []
    // )

    emit:multiqc_report = MULTIQC_PYTHON.out.report.map { _meta, report -> report }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
