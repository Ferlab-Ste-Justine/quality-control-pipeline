include { MOSDEPTH } from '../../../modules/nf-core/mosdepth/main'
include { COVERAGE_BY_GENE } from '../../../modules/local/coverage_by_gene/main'

workflow QC_COVERAGE_REGIONS {

    take:
    ch_bam_bai // channel: [ val(meta), bam, bai ]
    qc_regions // channel: path to qc regions file (BED format)
    fasta  // [optional]  val: path to reference fasta file

    main:

    ch_input_mosdepth = ch_bam_bai.combine(qc_regions)

    MOSDEPTH(
        ch_input_mosdepth,
        fasta.map { it -> [ [id:"fasta"], it] }
    )
    ch_mean_thresholds = MOSDEPTH.out.regions_bed.join(MOSDEPTH.out.thresholds_bed)

    COVERAGE_BY_GENE(ch_mean_thresholds)

    emit:
    reports  = COVERAGE_BY_GENE.out.report      // channel: [ path(report files) ]
}
