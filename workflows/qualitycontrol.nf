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
include { dragenFileMatchesSample } from '../subworkflows/local/utils_nfcore_quality-control-pipeline_pipeline'
include { FASTQ_QC               } from '../subworkflows/local/fastq_qc/main'
include { BAM_QC as BAM_QC_WGS   } from '../subworkflows/local/bam_qc/main'
include { BAM_QC as BAM_QC_TARGET  } from '../subworkflows/local/bam_qc/main'
include { VCF_QC as VCF_QC_WGS   } from '../subworkflows/local/vcf_qc/main'
include { VCF_QC as VCF_QC_TARGET  } from '../subworkflows/local/vcf_qc/main'
include { BAM_MERGE              } from '../subworkflows/local/bam_merge'
include { PED_FAMILY_CHECK       } from '../subworkflows/local/ped_family_check'
include { PED_MULTIQC_JOIN       } from '../subworkflows/local/ped_multiqc_join'
include { CRAM_SOMALIER          } from '../subworkflows/local/cram_somalier'
include { GATK4_BEDTOINTERVALLIST } from '../modules/nf-core/gatk4/bedtointervallist/main'
include { DRAGEN_COVERAGE_BY_GENE } from '../modules/local/dragen_coverage_by_gene/main'
include { UNTAR                  } from '../modules/nf-core/untar/main'

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
    ch_exons = params.exons_bed ? channel.value(file(params.exons_bed, checkIfExists:true)) : channel.value([])
    qc_regions_1 = params.qc_coverage_region_1 ? channel.value(file(params.qc_coverage_region_1, checkIfExists:true)) : channel.value([])
    qc_regions_2 = params.qc_coverage_region_2 ? channel.value(file(params.qc_coverage_region_2, checkIfExists:true)) : channel.value([])
    ch_targets = params.targets_bed ? channel.value(file(params.targets_bed, checkIfExists:true)) : channel.value([])
    ch_baits = params.baits_bed ? channel.value(file(params.baits_bed, checkIfExists:true)) : channel.value([])

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

    // meta.dragenMetrics (samplesheet column) selects, per sample and independently
    // per fileType row, whether that row is recomputed by the standard pipeline or
    // sourced from DRAGEN's own pre-computed metrics instead. A single run's
    // samplesheet can freely mix both kinds of samples/families.
    ch_fastq_branch = ch_samplesheet_parsed.fastq.branch { meta, files ->
        dragen: meta.dragenMetrics
        standard: true
    }
    ch_aln_branch = ch_samplesheet_parsed.aln.branch { meta, cram, crai ->
        dragen: meta.dragenMetrics
        standard: true
    }
    ch_vcf_branch = ch_samplesheet_parsed.vcf.branch { meta, vcf, tbi ->
        dragen: meta.dragenMetrics
        standard: true
    }

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        FASTQ QUALITY CONTROL
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */

    // If fastq files are provided, run FastQC, SeqFu, and ngsCheckMate
    FASTQ_QC ( ch_fastq_branch.standard, ncm_snp_pt.map { it -> [ [id:"fastq_checkmate"], it] } )
    ch_multiqc_files = ch_multiqc_files.mix(FASTQ_QC.out.reports)
    ch_versions = ch_versions.mix(FASTQ_QC.out.versions)

    /*
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        BAM/CRAM QUALITY CONTROL
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    // The interval_list is only read by PICARD_COLLECTWGSMETRICS (BAM_QC_WGS), so
    // only build it when at least one standard-branch WGS alignment row exists --
    // an all-DRAGEN (or all-targeted) run then schedules no GATK task at all.
    ch_has_standard_wgs_aln = ch_aln_branch.standard
        .filter { meta, _cram, _crai -> meta.sequencingType == 'WGS' }
        .first()
        .map { _row -> true }

    GATK4_BEDTOINTERVALLIST( ch_intervals.flatten()
                                .combine(ch_has_standard_wgs_aln)
                                .map { bed, _has_wgs -> [[id:"bed"], bed] },
                            ch_dict.map { it -> [[id:"dict"], it] })

    ch_interval_list = GATK4_BEDTOINTERVALLIST.out.interval_list
                        .map { _meta, intervals -> intervals }
                        .collect()
                        .ifEmpty([])

    // Merge multiple runs
    //
    // ----- SAMTOOLS MERGE -----
    //
    bam_to_merge = ch_aln_branch.standard
        .map { meta, cram, crai ->
        [ groupKey(meta.subMap('id', 'participant', 'sample', 'familyId', 'sex', 'sequencingType', 'status', 'relationship_to_proband', 'affected_status'), meta.n_lanes), cram, crai ]
    }
    .groupTuple()

    BAM_MERGE(bam_to_merge, ch_fasta, ch_fai)

    ch_versions = ch_versions.mix(BAM_MERGE.out.versions)

    // CRAM/CRAI source for somalier: standard-branch samples contribute their
    // merged file (from BAM_MERGE); dragen-branch samples were never merged, so
    // they feed somalier directly from the samplesheet, trimmed the same way
    // BAM_MERGE trims its GroupKey target so meta.n_lanes is absent — otherwise
    // the per-family groupTuple inside CRAM_SOMALIER emits each sample early
    // (count=1) and fails the per-family join.
    ch_cram_crai_source = BAM_MERGE.out.bam_bai
        .mix(
            ch_aln_branch.dragen.map { meta, cram, crai ->
                [ meta.subMap('id', 'participant', 'sample', 'familyId', 'sex', 'sequencingType', 'status', 'relationship_to_proband', 'affected_status'), cram, crai ]
            }
        )

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
        ch_dict,
        ch_intervals,
        qc_regions_1,
        qc_regions_2,
        ch_interval_list,
        [],
        [],
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
        ch_baits,
        ch_targets,
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

    // separate qc for targeted seq vs wgs
    ch_vcf_qc = ch_vcf_branch.standard
        .branch { meta, vcf, tbi ->
            wgs: meta.sequencingType == 'WGS'
            targeted: meta.sequencingType != 'WGS'
        }

    VCF_QC_WGS (
        ch_vcf_qc.wgs,
        ch_fasta,
        ch_intervals,
        channel.value([]),
        ch_exons
    )

    VCF_QC_TARGET (
        ch_vcf_qc.targeted,
        ch_fasta,
        ch_intervals,
        ch_targets,
        ch_exons
    )

    ch_multiqc_files = ch_multiqc_files.mix(VCF_QC_WGS.out.vcf_metrics)
    ch_multiqc_files = ch_multiqc_files.mix(VCF_QC_TARGET.out.vcf_metrics)

    ch_multiqc_files = ch_multiqc_files.mix(VCF_QC_WGS.out.vcf_stats)
    ch_multiqc_files = ch_multiqc_files.mix(VCF_QC_TARGET.out.vcf_stats)

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

        PED_FAMILY_CHECK(ch_somalier_input_ped, ch_samplesheet)

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
        DRAGEN METRICS (alternative to BAM_QC / VCF_QC, per sample)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    */
    // One row per DRAGEN-flagged sample. A sample can have several samplesheet rows
    // (e.g. a CRAM row and a VCF row), all repeating the same dragenMetrics value --
    // dedupe to one UNTAR/glob per sample.
    ch_dragen_samples = ch_samplesheet
        .filter { meta, _files -> meta.dragenMetrics }
        .map { meta, _files -> [ meta.sample, meta - meta.subMap('lane', 'runId', 'fileType') ] }
        .unique { it[0] }
        .map { _sample, meta -> meta }

    // meta.dragenMetrics accepts either a directory of that sample's DRAGEN metric
    // files, or a tarball of them (e.g. a CQGC-style per-sample "extra" bundle). A
    // tarball is extracted via UNTAR -- one task per sample needing it -- before the
    // same glob logic below applies uniformly to either shape via ch_dragen_dir.
    ch_dragen_input = ch_dragen_samples
        .map { meta -> [ meta, file(meta.dragenMetrics, checkIfExists: true) ] }
        .branch { meta, p ->
            tarball: p.name ==~ /.*\.(tar|tar\.gz|tgz|tar\.bz2|tbz2)$/
            directory: true
        }

    UNTAR ( ch_dragen_input.tarball )
    // UNTAR emits its version via the topic channel (topic: versions) -- picked up
    // automatically by the pipeline-wide channel.topic("versions") collection below.
    // if ch_dragen_input.tarball does not exist: zero UNTAR tasks get scheduled. Nextflow is purely reactive here: a process just does nothing when its input channel is empty. No task runs, no container gets pulled for it, no error.

    ch_dragen_dir = UNTAR.out.untar.mix(ch_dragen_input.directory) // mix both channels so that ch_dragen_dir ends up with one dir

    ch_dragen_files = ch_dragen_dir
        .flatMap { meta, dir ->
            // *.CollectHsMetrics.out: DRAGEN's own bundled Picard HsMetrics output
            // for exomes (no WGS equivalent). Not a DRAGEN-format CSV -- MultiQC's
            // own picard module content-detects it (by "HsMetrics" in the file, not
            // by filename) and parse_picard() already has the MEAN_TARGET_COVERAGE/
            // PCT_TARGET_BASES_20X fallback this needs; staging it is the whole fix.
            [
                file("${dir}/*.csv", checkIfExists: false),
                file("${dir}/**/*.csv", checkIfExists: false),
                file("${dir}/*.CollectHsMetrics.out", checkIfExists: false),
                file("${dir}/**/*.CollectHsMetrics.out", checkIfExists: false),
            ]
                .flatten()
                .findAll { f -> dragenFileMatchesSample(meta, f) }
                .collect { f -> [ meta, f ] }
        }

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_dragen_files.map { meta, f -> [ [id: meta.sample, sample: meta.sample, familyId: meta.familyId], f ] }
    )

    // DRAGEN per-region per-interval coverage → per-gene coverage TSV.
    ch_dragen_cov = ch_dragen_dir
        .flatMap { meta, dir ->
            [ file("${dir}/*_cov_report.bed", checkIfExists: false), file("${dir}/**/*_cov_report.bed", checkIfExists: false) ]
                .flatten()
                .findAll { f -> !f.name.endsWith('_read_cov_report.bed') }
                .findAll { f -> dragenFileMatchesSample(meta, f) }
                .collect { f ->
                    def region = f.name.tokenize('.')[1].replaceFirst(/_cov_report$/, '')
                    [ "${meta.sample}|${region}".toString(), meta, region, f ]
                }
        }

    ch_dragen_read_cov = ch_dragen_dir
        .flatMap { meta, dir ->
            [ file("${dir}/*_read_cov_report.bed", checkIfExists: false), file("${dir}/**/*_read_cov_report.bed", checkIfExists: false) ]
                .flatten()
                .findAll { f -> dragenFileMatchesSample(meta, f) }
                .collect { f ->
                    def region = f.name.tokenize('.')[1].replaceFirst(/_read_cov_report$/, '')
                    [ "${meta.sample}|${region}".toString(), f ]
                }
        }

    ch_dragen_gene_input = ch_dragen_cov
        .join(ch_dragen_read_cov)
        .map { _key, meta, region, cov, read_cov ->
            [ [id: meta.sample, sample: meta.sample, region: region, familyId: meta.familyId], cov, read_cov ]
        }

    DRAGEN_COVERAGE_BY_GENE(ch_dragen_gene_input)
    ch_multiqc_files = ch_multiqc_files.mix(DRAGEN_COVERAGE_BY_GENE.out.report)

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
    PED_MULTIQC_JOIN(ch_multiqc_input, ch_ped_by_id)
    ch_multiqc_input_with_ped = PED_MULTIQC_JOIN.out.multiqc_input_with_ped

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
