include { MOSDEPTH           } from '../../../modules/nf-core/mosdepth/main'
include { D4_COVERAGE_STATS as D4_COVERAGE_STATS_R1  } from '../d4_coverage_stats/main'
include { D4_COVERAGE_STATS as D4_COVERAGE_STATS_R2  } from '../d4_coverage_stats/main'

workflow QC_COVERAGE_REGIONS {

    take:
    ch_bam_bai // channel: [ val(meta), bam, bai ]
    fasta  // [optional] path to reference fasta file
    regions_target // [optional] path to target regions file if WXS
    qc_regions_1   // [optional] path to first qc regions file
    qc_regions_2   // [optional] path to second qc regions file

    main:

    ch_versions = channel.empty()
    ch_reports = channel.empty()

    ch_input_mosdepth = ch_bam_bai.combine(regions_target)

    MOSDEPTH(
        ch_input_mosdepth,
        fasta.map { it -> [ [id:"fasta"], it] }
    )
    ch_versions = ch_versions.mix(MOSDEPTH.out.versions)
    ch_reports = ch_reports.mix(MOSDEPTH.out.global_txt.map{it -> it[1]}.collect().ifEmpty([]))
    ch_reports = ch_reports.mix(MOSDEPTH.out.regions_txt.map{it -> it[1]}.collect().ifEmpty([]))

    D4_COVERAGE_STATS_R1(
        MOSDEPTH.out.per_base_d4,
        qc_regions_1.map { it -> [ [id:"qc_regions_1"], it] }
    )
    ch_reports = ch_reports.mix(D4_COVERAGE_STATS_R1.out.coverage_stats.map{it -> it[1]}.collect().ifEmpty([]))

    D4_COVERAGE_STATS_R2(
        MOSDEPTH.out.per_base_d4,
        qc_regions_2.map { it -> [ [id:"qc_regions_2"], it] }
    )
    ch_reports = ch_reports.mix(D4_COVERAGE_STATS_R2.out.coverage_stats.map{it -> it[1]}.collect().ifEmpty([]))

    emit:
    reports  = ch_reports                     // channel: [ path(report files) ]
    versions = ch_versions                     // channel: [ versions.yml ]
}
