include { FASTQC                 } from '../../../modules/nf-core/fastqc/main'
include { FASTQ_NGSCHECKMATE     } from '../../nf-core/fastq_ngscheckmate/main'

workflow FASTQ_QC {

    take:
    ch_fastq // channel: [ val(meta), [fastq1, fastq2?] ]
    ncm_snp_pt // channel: [ val(meta), path(snp_pt) ]

    main:
    ch_reports = channel.empty()
    ch_versions = channel.empty()

    FASTQ_NGSCHECKMATE (ch_fastq, ncm_snp_pt)
    ch_versions = ch_versions.mix(FASTQ_NGSCHECKMATE.out.versions)
    ch_reports = ch_reports.mix(FASTQ_NGSCHECKMATE.out.corr_matrix)
    ch_reports = ch_reports.mix(FASTQ_NGSCHECKMATE.out.matched)
    ch_reports = ch_reports.mix(FASTQ_NGSCHECKMATE.out.all)
    ch_reports = ch_reports.mix(FASTQ_NGSCHECKMATE.out.vaf)

    FASTQC (ch_fastq)
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    ch_reports = ch_reports.mix(FASTQC.out.zip)
    ch_reports = ch_reports.mix(FASTQC.out.html)

    emit:
    reports  = ch_reports                       // channel: [ val(meta), path(report) ]
    versions = ch_versions                      // channel: [ versions.yml ]
}
