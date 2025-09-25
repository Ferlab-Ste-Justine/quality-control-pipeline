include { FASTQC                 } from '../../../modules/nf-core/fastqc/main'
include { FQ_LINT                } from '../../../modules/nf-core/fq/lint/main'
include { SEQFU_CHECK            } from '../../../modules/local/seqfu/check/main'

workflow FASTQ_QC {

    take:
    ch_fastq // channel: [ val(meta), [fastq1, fastq2?] ]

    main:
    ch_reports = Channel.empty()
    ch_versions = Channel.empty()

    SEQFU_CHECK (
        ch_fastq
    )

    ch_reports = ch_reports.mix(SEQFU_CHECK.out.check.collect{it[1]})

    FQ_LINT (
        ch_fastq
    )

    ch_reports = ch_reports.mix(FQ_LINT.out.lint.collect{it[1]})

    FASTQC (
        ch_fastq
    )
    ch_reports = ch_reports.mix(FASTQC.out.zip.collect{it[1]})
    ch_reports = ch_reports.mix(FASTQC.out.html.collect{it[1]})


    ch_versions = ch_versions.mix(SEQFU_CHECK.out.versions.first())
    ch_versions = ch_versions.mix(FQ_LINT.out.versions.first())
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    emit:
    reports  = ch_reports                       // channel: [ *.tsv, *.html, *.zip ]
    versions = ch_versions                     // channel: [ versions.yml ]
}
