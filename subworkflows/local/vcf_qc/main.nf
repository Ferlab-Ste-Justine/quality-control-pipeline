include { BCFTOOLS_VIEW as COUNT_DEL	} from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as COUNT_INS	} from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as COUNT_SNV	} from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as HET_SNV	} from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as HOM_SNV	} from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as HOM_INDEL	} from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_VIEW as HET_INDEL	} from '../../../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_STATS as STATS_TSTV  } from '../../../modules/nf-core/bcftools/stats/main'

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
        def het_hom_ratio_snvs = het_snv.count_het_snvs > 0 ? het_snv.count_het_snvs / hom_snv.count_hom_snvs : 0
        def het_hom_ratio_indels = het_indel.count_het_indels > 0 ? het_indel.count_het_indels / hom_indel.count_hom_indels : 0
        def ins_dels_ratio = ins.count_ins > 0 ? ins.count_ins / dels.count_dels : 0
        [meta, dels + ins + snvs + het_snv + hom_snv + hom_indel + het_indel + [ het_hom_ratio_snvs: het_hom_ratio_snvs, het_hom_ratio_indels: het_hom_ratio_indels, ins_dels_ratio: ins_dels_ratio ] + tstv_ratio ]
    }

    emit:
    vcf_metrics = ch_vcf_metrics   // channel: [ val(meta), val(metrics), path(tstv_stats) ]
    vcf_stats = STATS_TSTV.out.stats   // channel: [ val(meta), path(stats) ]
}
