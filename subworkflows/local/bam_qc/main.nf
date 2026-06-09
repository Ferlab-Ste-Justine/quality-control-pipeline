include { SAMTOOLS_STATS } from '../../../modules/nf-core/samtools/stats/main'
include { PICARD_COLLECTWGSMETRICS } from '../../../modules/nf-core/picard/collectwgsmetrics/main'
include { PICARD_COLLECTHSMETRICS } from '../../../modules/nf-core/picard/collecthsmetrics/main'
include { PICARD_COLLECTQUALITYYIELDMETRICS } from '../../../modules/local/picard/collectqualityyieldmetrics/main'
include { VERIFYBAMID_VERIFYBAMID2 } from '../../../modules/nf-core/verifybamid/verifybamid2/main'
include { MOSDEPTH } from '../../../modules/nf-core/mosdepth/main'
include { QC_COVERAGE_REGIONS as QC_COVERAGE_REGIONS_R1 } from '../qc_coverage_regions/main'
include { QC_COVERAGE_REGIONS as QC_COVERAGE_REGIONS_R2 } from '../qc_coverage_regions/main'

workflow BAM_QC {

    take:

    ch_bam_bai           // channel: [ val(meta), path(bam/cram), path(bai/crai) ]
    ch_fasta             // channel: [optional] [ path(fasta) ]
    ch_fai               // channel: [optional] [ path(fai) ]
    ch_dict              // channel: [optional] [ path(dict) ]  - needed for HsMetrics BedToIntervalList
    ch_intervals         // channel: [optional] [ path(intervals) ]
    qc_regions_1         // channel: [optional] path to first qc regions file
    qc_regions_2         // channel: [optional] path to second qc regions file
    ch_intervals_list    // channel: [optional] [ path(interval_list) ]  - used by CollectWgsMetrics
    ch_bait_intervals    // channel: [optional] [ path(bed or interval_list) ] - HsMetrics bait
    ch_target_intervals  // channel: [optional] [ path(bed or interval_list) ] - HsMetrics target
    ch_svd_in            // channel: [optional] [ path(svd_ud), path(svd_mu), path(svd_bed) ]
    val_seq_type         // string:  'WGS' triggers CollectWgsMetrics, anything else uses CollectHsMetrics

    main:

    ch_versions = channel.empty()
    ch_reports = channel.empty()

    //
    // ----- SAMTOOLS STATS -----
    //

    SAMTOOLS_STATS( ch_bam_bai,
        ch_fasta.map { it -> [ [id:"fasta"], it] }
    )
    ch_reports = ch_reports.mix(SAMTOOLS_STATS.out.stats)
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions)

    //
    // ----- PICARD COVERAGE METRICS: WgsMetrics for WGS, HsMetrics for targeted/WXS -----
    //
    if (val_seq_type == 'WGS') {
        PICARD_COLLECTWGSMETRICS( ch_bam_bai,
            ch_fasta.map { it -> [ [id:"fasta"], it] },
            ch_fai.map { it -> [ [id:"fai"], it] },
            ch_intervals_list
        )
        ch_reports = ch_reports.mix(PICARD_COLLECTWGSMETRICS.out.metrics)
        ch_versions = ch_versions.mix(PICARD_COLLECTWGSMETRICS.out.versions)
    }
    else {
        // HsMetrics input expects [meta, bam, bai, bait_intervals, target_intervals]
        // ch_bait/target_intervals may be a BED; the module converts to interval_list internally.
        ch_hsmetrics_in = ch_bam_bai
            .combine(ch_bait_intervals.toList())
            .combine(ch_target_intervals.toList())
        PICARD_COLLECTHSMETRICS( ch_hsmetrics_in,
            ch_fasta.map { it -> [ [id:"fasta"], it] },
            ch_fai.map   { it -> [ [id:"fai"],   it] },
            ch_dict.map  { it -> [ [id:"dict"],  it] },
            [ [id:"gzi"], [] ]
        )
        ch_reports = ch_reports.mix(PICARD_COLLECTHSMETRICS.out.metrics)
    }

    //
    // ----- PICARD COLLECTQUALITYYIELDMETRICS (GA4GH yield_bp_q30) -----
    //
    PICARD_COLLECTQUALITYYIELDMETRICS( ch_bam_bai,
        ch_fasta.map { it -> [ [id:"fasta"], it] },
        ch_fai.map { it -> [ [id:"fai"], it] }
    )
    ch_reports = ch_reports.mix(PICARD_COLLECTQUALITYYIELDMETRICS.out.metrics)

    //
    // ----- MOSDEPTH -----
    //
    ch_input_mosdepth = ch_bam_bai.combine(ch_intervals.toList())
    MOSDEPTH(
        ch_input_mosdepth,
        ch_fasta.map { it -> [ [id:"fasta"], it] }
    )

    ch_reports = ch_reports.mix(MOSDEPTH.out.global_txt)
    ch_reports = ch_reports.mix(MOSDEPTH.out.regions_txt)

    // Coverage by gene for specified regions (if provided)
    if (params.qc_coverage_region_1) {
        QC_COVERAGE_REGIONS_R1(
            ch_bam_bai,
            qc_regions_1,
            ch_fasta
        )
        ch_reports = ch_reports.mix(QC_COVERAGE_REGIONS_R1.out.reports)
    }

    if (params.qc_coverage_region_2) {
        QC_COVERAGE_REGIONS_R2(
            ch_bam_bai,
            qc_regions_2,
            ch_fasta
        )
        ch_reports = ch_reports.mix(QC_COVERAGE_REGIONS_R2.out.reports)
    }

    //
    // ----- VERIFYBAMID2 - Contamination -----
    //
    VERIFYBAMID_VERIFYBAMID2(
        ch_bam_bai, ch_svd_in, [], ch_fasta)

    ch_reports = ch_reports.mix(VERIFYBAMID_VERIFYBAMID2.out.self_sm)
    ch_versions = ch_versions.mix(VERIFYBAMID_VERIFYBAMID2.out.versions)

    emit:
    reports = ch_reports                     // channel: [ val(meta), path(report) ]
    versions = ch_versions                   // channel: [ versions.yml ]
}
