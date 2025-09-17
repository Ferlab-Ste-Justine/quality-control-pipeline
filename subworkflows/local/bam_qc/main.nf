include { QC_COVERAGE_REGIONS           } from '../qc_coverage_regions/main'
include { QUALIMAP_BAMQC     } from '../../../modules/nf-core/qualimap/bamqc/main'
include { QUALIMAP_BAMQCCRAM } from '../../../modules/nf-core/qualimap/bamqccram/main'
include { VERIFYBAMID_VERIFYBAMID2 } from '../../../modules/nf-core/verifybamid/verifybamid2/main'

workflow BAM_QC {

    take:

    ch_bam_bai      // channel: [ val(meta), path(bam/cram), path(bai/crai) ]
    ch_fasta        // channel: [optional] [ path(fasta) ]
    ch_fai          // channel: [optional] [ path(fai) ]
    ch_intervals    // channel: [optional] [ path(intervals) ]
    qc_regions_1 // channel: [optional] path to first qc regions file
    qc_regions_2 // channel: [optional] path to second qc regions file
    ch_svd_in       // channel: [optional] [ path(svd_ud), path(svd_mu), path(svd_bed) ]

    main:

    ch_versions = channel.empty()
    ch_reports = channel.empty()

    //
    // ----- COVERAGE -----
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
    // ----- QUALIMAP -----
    //

    // check if data_type in meta is cram or bam
    ch_input_qualimap = ch_bam_bai
        .branch { meta, aln, idx ->
            cram: aln.getExtension() == 'cram'
            bam: aln.getExtension() == 'bam'
            }

    QUALIMAP_BAMQCCRAM (
        ch_input_qualimap.cram,
        ch_intervals,
        ch_fasta, ch_fai)

    ch_reports = ch_reports.mix(QUALIMAP_BAMQCCRAM.out.results.map{it[1]}.collect().ifEmpty([]))

    QUALIMAP_BAMQC (
        ch_input_qualimap.bam.map { meta, bam, _idx -> [meta, bam] },
        ch_intervals )

    ch_reports = ch_reports.mix(QUALIMAP_BAMQC.out.results.map{it[1]}.collect().ifEmpty([]))

    //
    // ----- VERIFYBAMID2 - Contamination -----
    //
    VERIFYBAMID_VERIFYBAMID2(
        ch_bam_bai, ch_svd_in, [], ch_fasta)

    ch_reports = ch_reports.mix(VERIFYBAMID_VERIFYBAMID2.out.self_sm.map{it[1]}.collect().ifEmpty([]))

    // Collect versions
    ch_versions = ch_versions.mix(QC_COVERAGE_REGIONS.out.versions)
    ch_versions = ch_versions.mix(QUALIMAP_BAMQCCRAM.out.versions)
    ch_versions = ch_versions.mix(QUALIMAP_BAMQC.out.versions)
    ch_versions = ch_versions.mix(VERIFYBAMID_VERIFYBAMID2.out.versions)

    emit:
    reports = ch_reports                     // channel: [ path(report1), path(report2), ... ]
    versions = ch_versions                     // channel: [ versions.yml ]
}
