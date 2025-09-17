//
// Workflow that checks if internal sampleID in BAM/CRAM matches with sample_registration id and renames sampleid if not
//
include { SAMTOOLS_SAMPLES     } from '../../../modules/local/samtools/samples/main'
include { SAMTOOLS_INDEX     } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_REHEADER as SAMTOOLS_REHEADER_SAMPLE   } from '../../../modules/local/samtools/reheader/main'

workflow BAM_ID_REPAIR {

    take:
    ch_input // channel: [ mandatory ] meta, bam/cram, bai/crai
    ch_fasta // channel: [ optional ] 

    main:
    ch_versions = Channel.empty()

    SAMTOOLS_SAMPLES(ch_input, ch_fasta, [])

    branched_bams = ch_input
        .join(SAMTOOLS_SAMPLES.out.output)
            .map{ meta, bam, bai, query_result ->
                assert query_result.countLines() < 2 : 'More than one sample in bam file'
                def sample_name = query_result.text.split().first() // txt file with the result
                [ meta, bam, bai, sample_name ]
            }
            .branch { meta, bam, bai, sample_name ->
                reheader: (sample_name != meta.participant_sample) && params.skip_reheader == false
                direct: (sample_name == meta.participant_sample) || params.skip_reheader == true
            }
    
    // create input to vcf reheader option --samples
    ch_reheader_input =  branched_bams.reheader
                            .map { meta, bam, _bai, _sample_name ->
                                [ meta, bam ]
                            }

    // Edit Sample ID in bam/cram
    SAMTOOLS_REHEADER_SAMPLE(ch_reheader_input)

    // index new reheaded file
    SAMTOOLS_INDEX(SAMTOOLS_REHEADER_SAMPLE.out.bam)

    index_ch = SAMTOOLS_INDEX.out.bai
        .mix(SAMTOOLS_INDEX.out.crai)
        .mix(SAMTOOLS_INDEX.out.csi)

    bam_bai = SAMTOOLS_REHEADER_SAMPLE.out.bam
                .join( index_ch )
                .mix( branched_bams.direct
                        .map { meta, bam, idx, _sample_name ->
                        [ meta, bam, idx ] } )

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(SAMTOOLS_REHEADER_SAMPLE.out.versions.first())
    ch_versions = ch_versions.mix(SAMTOOLS_SAMPLES.out.versions.first())
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    emit:
    bam_bai                     // channel: [ val(meta), path(bam/cram), path(bai/crai) ]
    versions = ch_versions       // channel: [ path(versions.yml) ]
}