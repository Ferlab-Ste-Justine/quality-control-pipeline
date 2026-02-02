include { BCFTOOLS_VIEW as COUNT_DEL	} from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as COUNT_INS	} from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as COUNT_SNV	} from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as HET_SNV	} from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as HOM_SNV	} from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as HOM_INDEL	} from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as HET_INDEL	} from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_STATS as STATS_TSTV  } from '../../../modules/nf-core/bcftools/stats/main'

include { VCF_METRICS                       } from '../../../modules/local/vcf_metrics/main'

workflow VCF_QC {

    take:

    ch_vcf_tbi      // channel: [ val(meta), path(bam/cram), path(bai/crai) ]
    ch_fasta        // channel: [optional] [ path(fasta) ]
    ch_intervals    // channel: [optional] [ path(intervals) ]
    ch_targets     // channel: [optional] path to targets file ]
    ch_exons       // channel: [optional] path to exons file ]

    main:

    COUNT_DEL ( ch_vcf_tbi, ch_intervals, ch_targets, [] )
    ch_dels = COUNT_DEL.out.vcf
        .map { meta, vcf ->
            def dels = vcf.countLines()
            [meta, [count_dels: dels]  ]
        }

    COUNT_INS ( ch_vcf_tbi, ch_intervals, ch_targets, [] )
    ch_ins = COUNT_INS.out.vcf
        .map { meta, vcf ->
            def ins = vcf.countLines()
            [meta, [count_ins: ins] ]
        }

    COUNT_SNV ( ch_vcf_tbi, ch_intervals, ch_targets, [] )
    ch_snvs = COUNT_SNV.out.vcf
        .map { meta, vcf ->
            def snvs = vcf.countLines()
            [meta, [count_snvs: snvs] ]
        }

    HET_SNV ( ch_vcf_tbi, ch_intervals, ch_targets, [] )
    ch_het_snvs = HET_SNV.out.vcf
        .map { meta, vcf ->
            def het_snvs = vcf.countLines()
            [meta, [count_het_snvs: het_snvs] ]
        }

    HOM_SNV ( ch_vcf_tbi, ch_intervals, ch_targets, [] )
    ch_hom_snvs = HOM_SNV.out.vcf
        .map { meta, vcf ->
            def hom_snvs = vcf.countLines()
            [meta, [count_hom_snvs: hom_snvs] ]
        }

    HET_INDEL ( ch_vcf_tbi, ch_intervals, ch_targets, [] )
    ch_het_indels = HET_INDEL.out.vcf
        .map { meta, vcf ->
            def het_indels = vcf.countLines()
            [meta, [count_het_indels: het_indels] ]
        }

    HOM_INDEL ( ch_vcf_tbi, ch_intervals, ch_targets, [] )
    ch_hom_indels = HOM_INDEL.out.vcf
        .map { meta, vcf ->
            def hom_indels = vcf.countLines()
            [meta, [count_hom_indels: hom_indels] ]
        }

    STATS_TSTV ( ch_vcf_tbi,
                ch_intervals.map { it -> [[id: "intervals"],it]},
                ch_targets.map { it -> [[id: "targets"],it]},
                [[],[]],
                ch_exons.map { it -> [[id: "exons"],it]},
                ch_fasta.map { it -> [[id: "fasta"],it]}
    )
    ch_tstv = STATS_TSTV.out.stats
        .map { meta, stats ->
            def tstv_line = stats.readLines().find { it -> it.startsWith('TSTV') }
            def tstv_fields = tstv_line.tokenize('\t')
            def tstv_ratio = tstv_fields[4]
            [meta, [tstv_ratio: tstv_ratio]]
        }

    ch_vcf_metrics = ch_dels
    .join(ch_ins)
    .join(ch_snvs)
    .join(ch_het_snvs)
    .join(ch_hom_snvs)
    .join(ch_hom_indels)
    .join(ch_het_indels)
    .join(ch_tstv)
    .map { meta, dels, ins, snvs, het_snv, hom_snv, hom_indel, het_indel, tstv_ratio ->
        [meta, dels + ins + snvs + het_snv + hom_snv + hom_indel + het_indel + tstv_ratio ]
    }

    VCF_METRICS ( ch_vcf_metrics )

    emit:
    vcf_metrics = VCF_METRICS.out.json   // channel: [ val(meta), val(metrics), path(metrics_file) ]
    vcf_stats = STATS_TSTV.out.stats   // channel: [ val(meta), path(stats) ]
}
