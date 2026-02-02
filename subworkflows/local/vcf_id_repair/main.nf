//
// Workflow that checks if internal sampleID in VCF matches with sample_registration id and renames sampleid if not
//
include { BCFTOOLS_QUERY as BCFTOOLS_QUERY_SAMPLE     } from '../../../modules/nf-core/bcftools/query/main'
include { BCFTOOLS_REHEADER as BCFTOOLS_REHEADER_SAMPLE   } from '../../../modules/nf-core/bcftools/reheader/main'

process createReheaderSampleInput{
    input:
        tuple val(meta), path(vcf), val(sample_name)

    output:
        tuple val(meta), path(vcf), [], path("${meta.id}_renameVCF.tsv")

    exec:
    f = file(["${task.workDir}","${meta.id}_renameVCF.tsv"].join(File.separator))
    f.text = ["${sample_name}","${meta.id}"].join("\t")
}

workflow VCF_ID_REPAIR {

    take:
    ch_input // channel: [ mandatory ] meta, vcf, tbi

    main:
    ch_versions = Channel.empty()

    BCFTOOLS_QUERY_SAMPLE(ch_input, [], [], []) // --list-samples in task.ext.args
    branched_vcfs = ch_input
        .join(BCFTOOLS_QUERY_SAMPLE.out.output)
            .map{ meta, vcf, tbi, query_result ->
                def sample_name = query_result.text.trim() // txt file with the result
                [ meta, vcf, tbi, sample_name ]
            }
            .branch { meta, vcf, tbi, sample_name ->
                reheader: (sample_name != meta.id)
                direct: (sample_name == meta.id)
            }

    branched_vcfs.reheader
        | map { meta, vcf, _tbi, sample_name -> [ meta, vcf, sample_name ] }
        | createReheaderSampleInput

    // Edit Sample ID in vcf
    BCFTOOLS_REHEADER_SAMPLE(createReheaderSampleInput, [[:],[]])
    vcf_tbi = BCFTOOLS_REHEADER_SAMPLE.out.vcf
                .join(BCFTOOLS_REHEADER_SAMPLE.out.index)
                .mix( branched_vcfs.direct
                        .map { meta, vcf, tbi, _sample_name ->
                        [ meta, vcf, tbi ]
                        })

    // Gather versions of all tools used
    ch_versions = ch_versions.mix(
        BCFTOOLS_QUERY_SAMPLE.out.versions.first(),
        BCFTOOLS_REHEADER_SAMPLE.out.versions.first()
    )

    emit:
    vcf_tbi                      // channel: [ val(meta), path(vcf), path(tbi) ]
    versions = ch_versions       // channel: [ path(versions.yml) ]
}
