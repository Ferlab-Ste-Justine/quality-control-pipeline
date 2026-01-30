include { FASTQC                 } from '../../../modules/nf-core/fastqc/main'
include { FASTQ_NGSCHECKMATE     } from '../../nf-core/fastq_ngscheckmate/main'

workflow FASTQ_QC {

    take:
    ch_fastq // channel: [ val(meta), [fastq1, fastq2?] ]
    ncm_snp_pt // channel: [ val(meta), path(snp_pt) ]

    main:
    ch_reports = Channel.empty()
    ch_versions = Channel.empty()

    FASTQ_NGSCHECKMATE (ch_fastq, ncm_snp_pt)
    ch_versions = ch_versions.mix(FASTQ_NGSCHECKMATE.out.versions)
    ch_reports = ch_reports.mix(FASTQ_NGSCHECKMATE.out.corr_matrix.map { _meta, report -> report })
    ch_reports = ch_reports.mix(FASTQ_NGSCHECKMATE.out.matched.map { _meta, report -> report })
    ch_reports = ch_reports.mix(FASTQ_NGSCHECKMATE.out.all.map { _meta, report -> report })
    ch_reports = ch_reports.mix(FASTQ_NGSCHECKMATE.out.vaf.map { _meta, report -> report })
    ch_reports = ch_reports.mix(FASTQ_NGSCHECKMATE.out.pdf.map { _meta, report -> report }.ifEmpty([]))

    FASTQC (ch_fastq)
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())
    ch_reports = ch_reports.mix(FASTQC.out.zip.collect{it[1]})
    ch_reports = ch_reports.mix(FASTQC.out.html.collect{it[1]})

    emit:
    reports  = ch_reports                       // channel: [ *.tsv, *.html, *.zip ]
    versions = ch_versions                     // channel: [ versions.yml ]
}
