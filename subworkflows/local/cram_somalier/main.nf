include { SOMALIER_EXTRACT } from '../../../modules/nf-core/somalier/extract/main'
include { SOMALIER_RELATE  } from '../../../modules/nf-core/somalier/relate/main'
include { SAMTOOLS_INDEX      } from '../../../modules/nf-core/samtools/index/main'

// BASED OF nf-core/vcf_extract_relate_somalier
// This workflow is adapted to work with CRAM files instead of VCF files

workflow CRAM_SOMALIER {

    take:
        ch_crams                // channel: [mandatory] [ val(meta), path(cram), path(crai), val(count) ]
        ch_fasta                // channel: [mandatory] [ val(meta), path(fasta) ]
        ch_fasta_fai            // channel: [mandatory] [ val(meta), path(fai) ]
        ch_somalier_sites       // channel: [mandatory] [ path(somalier_sites_vcf) ]
        ch_peds                 // channel: [mandatory] [ val(meta), path(ped) ]
        ch_sample_groups        // channel: [optional]  [ path(txt) ]
        val_common_id           // string:  [optional]  A common identifier for the samples that need to be related. - Family ID for example

    main:

    ch_versions = Channel.empty()

    ch_input = ch_crams
        .branch { meta, cram, crai, _count ->
            crai: crai != []
                return [ meta, cram, crai ]
            no_crai: crai == []
                return [ meta, cram ]
        }

    SAMTOOLS_INDEX ( ch_input.no_crai )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())

    ch_somalierextract_input = ch_input.no_crai
        .join(SAMTOOLS_INDEX.out.crai) //.bai
        .mix(ch_input.crai)
        .map { meta, cram, crai ->
            [ meta, cram, crai, meta.samplename_somalier ]
        }

    SOMALIER_EXTRACT(
        ch_somalierextract_input,
        ch_fasta,
        ch_fasta_fai,
        ch_somalier_sites
    )

    ch_versions = ch_versions.mix(SOMALIER_EXTRACT.out.versions)

     // Prepare input for CRAM Somalier subworkflow

    if (params.somalier_perfamily) {
        ch_somalierrelate_input = SOMALIER_EXTRACT.out.extract
            .join(ch_crams, failOnDuplicate: true, failOnMismatch: true)
            .map { meta, extract, _cram, _crai, count ->
                def new_meta = val_common_id ? meta + [id:meta[val_common_id]] : meta
                [ count ? groupKey(new_meta, count): new_meta, extract ]
            }
            .groupTuple()
            .join(ch_peds, failOnDuplicate: true, failOnMismatch: true)
            .map { meta, extract, ped ->
                def extract2 = extract[0] instanceof ArrayList ? extract[0] : extract
                def sorted_extract = extract2.sort { a, b -> file(a).name <=> file(b).name }
                def new_meta = meta instanceof nextflow.extension.GroupKey ? meta.target : meta
                [ new_meta, sorted_extract, ped ]
            } // Sort and flatten the extract list, remove the GroupKey wrapper if present
    }
    else {
        ch_somalierrelate_input = SOMALIER_EXTRACT.out.extract
            .map { _meta, extract_path -> extract_path}
            .collect()
            .map { paths_list -> [ [id: 'Cohort'], paths_list ] }
            .combine(ch_peds.map { _meta, ped -> ped }.toList())
    }

    SOMALIER_RELATE(
        ch_somalierrelate_input,
        ch_sample_groups
    )

    ch_versions = ch_versions.mix(SOMALIER_RELATE.out.versions)

    emit:
    extract        = SOMALIER_EXTRACT.out.extract       // channel: [ val(meta), path(extract) ]
    html           = SOMALIER_RELATE.out.html           // channel: [ val(meta), path(html) ]
    pairs_tsv      = SOMALIER_RELATE.out.pairs_tsv      // channel: [ val(meta), path(tsv) ]
    samples_tsv    = SOMALIER_RELATE.out.samples_tsv    // channel: [ val(meta), path(tsv) ]
    versions       = ch_versions                        // channel: [ path(versions.yml) ]
}
