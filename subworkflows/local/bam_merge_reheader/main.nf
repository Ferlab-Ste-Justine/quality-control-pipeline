include { SAMTOOLS_SAMPLES     } from '../../../modules/local/samtools/samples/main'
include { SAMTOOLS_REHEADER  } from '../../../modules/local/samtools/reheader/main'
include { SAMTOOLS_MERGE     } from '../../../modules/nf-core/samtools/merge/main'
include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index/main'

workflow BAM_MERGE_REHEADER {

    take:
    // TODO nf-core: edit input (take) channels
    ch_bam      // channel: [mandatory] meta, bam/cram
    fasta       // channel: [mandatory for cram] fasta
    fasta_fai   // channel: [mandatory for cram] fai

    main:
    ch_versions = Channel.empty()

    // TODO: define channel of files grouped by sample
    bam_to_merge = ch_bam.branch { meta, bam ->
        merge: bam.size() > 1
        direct: bam.size() <= 1
    }

    // merge mapped files by samplex
    SAMTOOLS_MERGE( bam_to_merge,
                    fasta.map{ it -> [ [ id:'fasta' ], it ] },
                    fasta_fai.map{ it -> [ [ id:'fasta_fai' ], it ] } )

    // TODO : allow cram or bam output
    ch_out_merge = SAMTOOLS_MERGE.out.cram ? SAMTOOLS_MERGE.out.cram : SAMTOOLS_MERGE.out.bam
    ch_versions = ch_versions.mix(SAMTOOLS_MERGE.out.versions.first())

    // TODO nf-core: substitute modules here for the modules of your subworkflow
    SAMTOOLS_REHEADER ( ch_out_merge )
    ch_versions = ch_versions.mix(SAMTOOLS_REHEADER.out.versions.first())

    SAMTOOLS_INDEX ( SAMTOOLS_REHEADER.out.bam )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    emit:
    // TODO nf-core: edit emitted channels
    bam      = SAMTOOLS_REHEADER.out.bam           // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    csi      = SAMTOOLS_INDEX.out.csi          // channel: [ val(meta), [ csi ] ]

    versions = ch_versions                     // channel: [ versions.yml ]
}

