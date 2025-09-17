//
// Workflow that merges multiple BAM/CRAM files per sample into a single file
//
include { SAMTOOLS_MERGE     } from '../../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index/main'

workflow BAM_MERGE_REHEADER {

    take:
    ch_bam      // channel: [mandatory] meta, [bam/cram], [bai/crai]
    fasta       // channel: [mandatory for cram] fasta
    fasta_fai   // channel: [mandatory for cram] fai

    main:
    ch_versions = Channel.empty()

    // branch bams that need merging or not. bam is a list of files
    bam_to_merge = ch_bam.branch { meta, bam, bai ->
        merge: bam.size() > 1
            return [meta , bam]
        direct: bam.size() <= 1
            return [meta, bam, bai]
    }

    // merge mapped files by samplex
    SAMTOOLS_MERGE( bam_to_merge.merge,
                    fasta.map{ it -> [ [ id:'fasta' ], it ] },
                    fasta_fai.map{ it -> [ [ id:'fasta_fai' ], it ] } )

    ch_out_merge = SAMTOOLS_MERGE.out.cram ?: SAMTOOLS_MERGE.out.bam
    ch_out_merge = SAMTOOLS_MERGE.out.bam
        .mix(SAMTOOLS_MERGE.out.cram)

    SAMTOOLS_INDEX ( ch_out_merge )

    ch_out_index = SAMTOOLS_INDEX.out.bai
        .mix(SAMTOOLS_INDEX.out.crai)
        .mix(SAMTOOLS_INDEX.out.csi)

    bam_bai = ch_out_merge
        .join(ch_out_index, failOnMismatch:true, failOnDuplicate:true)
        .mix(bam_to_merge.direct)

    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions.first())
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    emit:

    bam_bai      = bam_bai // channel: [ val(meta), bam, bai ] or [ val(meta), cram, crai ]
    versions = ch_versions                     // channel: [ versions.yml ]
}

