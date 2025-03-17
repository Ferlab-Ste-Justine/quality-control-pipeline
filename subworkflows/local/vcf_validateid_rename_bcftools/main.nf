include { BCFTOOLS_QUERY as BCFTOOLS_QUERY_SAMPLE     } from '../../../modules/nf-core/bcftools/query/main'
include { BCFTOOLS_REHEADER    } from '../../../modules/nf-core/bcftools/reheader/main'

workflow VCF_VALIDATEID_RENAME_BCFTOOLS {

    take:
    ch_input // channel: [ mandatory ] meta, vcf, tbi

    main:
    ch_versions = Channel.empty()

    BCFTOOLS_QUERY_SAMPLE(ch_input, [], [], []) // --list-samples in task.ext.args
    branched_vcfs = ch_input
        .join(BCFTOOLS_QUERY_SAMPLE.out.output)
            .map{ meta, vcf, tbi, query_result ->
                def sample_name = query_result.text.trim() // txt file with the result
                def needs_reheader = (sample_name != meta.id)
                [ meta, vcf, tbi, needs_reheader, sample_name ]
            }
            .branch { meta, vcf, tbi, needs_reheader, sample_name ->
                reheader: needs_reheader
                direct: !needs_reheader
            }
    // create input to vcf reheader option --samples
    ch_reheader_input =  branched_vcfs.reheader
                            .map { meta, vcf, _tbi, _needs_reheader, sample_name ->
                                def rename_tsv = file("${meta.id}_renameVCF.tsv")
                                rename_tsv.text = "${sample_name}\t${meta.id}"
                                [ meta, vcf , [], rename_tsv ]
                            }
    // Edit Sample ID in vcf
    BCFTOOLS_REHEADER(ch_reheader_input, [[:],[]])
    vcf_tbi = BCFTOOLS_REHEADER.out.vcf
                .join(BCFTOOLS_REHEADER.out.index)
                .mix( branched_vcfs.direct
                        .map { meta, vcf, tbi, _needs_reheader, _sample_name ->
                        [ meta, vcf, tbi ] } )

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(BCFTOOLS_QUERY_SAMPLE.out.versions.first())
    ch_versions = ch_versions.mix(BCFTOOLS_REHEADER.out.versions.first())

    emit:
    vcf_tbi                      // channel: [ val(meta), path(vcf), path(tbi) ]
    versions = ch_versions       // channel: [ path(versions.yml) ]
}