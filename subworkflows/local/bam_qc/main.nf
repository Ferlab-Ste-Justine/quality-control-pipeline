include { SAMTOOLS_STATS } from '../../../modules/nf-core/samtools/stats/main'
include { PICARD_COLLECTWGSMETRICS } from '../../../modules/nf-core/picard/collectwgsmetrics/main'
include { VERIFYBAMID_VERIFYBAMID2 } from '../../../modules/nf-core/verifybamid/verifybamid2/main'
include { QC_COVERAGE_REGIONS } from '../qc_coverage_regions/main'
workflow BAM_QC {

    take:

    ch_bam_bai      // channel: [ val(meta), path(bam/cram), path(bai/crai) ]
    ch_fasta        // channel: [optional] [ path(fasta) ]
    ch_fai          // channel: [optional] [ path(fai) ]
    ch_intervals    // channel: [optional] [ path(intervals) ]
    qc_regions_1 // channel: [optional] path to first qc regions file
    qc_regions_2 // channel: [optional] path to second qc regions file
    ch_intervals_list // channel: [optional] [ path(interval_list) ]
    ch_svd_in       // channel: [optional] [ path(svd_ud), path(svd_mu), path(svd_bed) ]

    main:

    ch_versions = channel.empty()
    ch_reports = channel.empty()

    //
    // ----- SAMTOOLS STATS -----
    //

    SAMTOOLS_STATS( ch_bam_bai,
        ch_fasta.map { it -> [ [id:"fasta"], it] }
    )
    ch_reports = ch_reports.mix(SAMTOOLS_STATS.out.stats.map{it[1]}.collect())

    //
    // ----- PICARD COLLECTWGMETRICS -----
    //

    PICARD_COLLECTWGSMETRICS( ch_bam_bai,
        ch_fasta.map { it -> [ [id:"fasta"], it] },
        ch_fai.map { it -> [ [id:"fai"], it] },
        ch_intervals_list
    )

    ch_reports = ch_reports.mix(PICARD_COLLECTWGSMETRICS.out.metrics.map{it[1]}.collect())

    //
    // ----- MOSDEPTH -----
    //
    QC_COVERAGE_REGIONS (
        ch_bam_bai,
        ch_fasta,
        ch_intervals,
        qc_regions_1,
        qc_regions_2
    )

    ch_reports = ch_reports.mix(QC_COVERAGE_REGIONS.out.reports)

    //
    // ----- VERIFYBAMID2 - Contamination -----
    //
    VERIFYBAMID_VERIFYBAMID2(
        ch_bam_bai, ch_svd_in, [], ch_fasta)

    ch_reports = ch_reports.mix(VERIFYBAMID_VERIFYBAMID2.out.self_sm.map{it[1]}.collect())

    // Collect versions
    ch_versions = ch_versions.mix(QC_COVERAGE_REGIONS.out.versions)
    ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions)
    ch_versions = ch_versions.mix(PICARD_COLLECTWGSMETRICS.out.versions)
    ch_versions = ch_versions.mix(VERIFYBAMID_VERIFYBAMID2.out.versions)

    emit:
    reports = ch_reports                     // channel: [ path(report1), path(report2), ... ]
    versions = ch_versions                     // channel: [ versions.yml ]
}
