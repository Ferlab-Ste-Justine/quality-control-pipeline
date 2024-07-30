#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ferlab/mypipeline-nfcore
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/ferlab/mypipeline
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
params.ascat_alleles           = getGenomeAttribute('ascat_alleles')
params.ascat_genome            = getGenomeAttribute('ascat_genome')
params.ascat_loci              = getGenomeAttribute('ascat_loci')
params.ascat_loci_gc           = getGenomeAttribute('ascat_loci_gc')
params.ascat_loci_rt           = getGenomeAttribute('ascat_loci_rt')
params.bwa                     = getGenomeAttribute('bwa')
params.bwamem2                 = getGenomeAttribute('bwamem2')
params.cf_chrom_len            = getGenomeAttribute('cf_chrom_len')
params.chr_dir                 = getGenomeAttribute('chr_dir')
params.dbsnp                   = getGenomeAttribute('dbsnp')
params.dbsnp_tbi               = getGenomeAttribute('dbsnp_tbi')
params.dbsnp_vqsr              = getGenomeAttribute('dbsnp_vqsr')
params.dict                    = getGenomeAttribute('dict')
params.dragmap                 = getGenomeAttribute('dragmap')
params.fasta                   = getGenomeAttribute('fasta')
params.fasta_fai               = getGenomeAttribute('fasta_fai')
params.germline_resource       = getGenomeAttribute('germline_resource')
params.germline_resource_tbi   = getGenomeAttribute('germline_resource_tbi')
params.intervals               = getGenomeAttribute('intervals')
params.known_indels            = getGenomeAttribute('known_indels')
params.known_indels_tbi        = getGenomeAttribute('known_indels_tbi')
params.known_indels_vqsr       = getGenomeAttribute('known_indels_vqsr')
params.known_snps              = getGenomeAttribute('known_snps')
params.known_snps_tbi          = getGenomeAttribute('known_snps_tbi')
params.known_snps_vqsr         = getGenomeAttribute('known_snps_vqsr')
params.mappability             = getGenomeAttribute('mappability')
params.ngscheckmate_bed        = getGenomeAttribute('ngscheckmate_bed')
params.pon                     = getGenomeAttribute('pon')
params.pon_tbi                 = getGenomeAttribute('pon_tbi')
params.sentieon_dnascope_model = getGenomeAttribute('sentieon_dnascope_model')
params.snpeff_db               = getGenomeAttribute('snpeff_db')
params.snpeff_genome           = getGenomeAttribute('snpeff_genome')
params.vep_cache_version       = getGenomeAttribute('vep_cache_version')
params.vep_genome              = getGenomeAttribute('vep_genome')
params.vep_species             = getGenomeAttribute('vep_species')

aligner = params.aligner

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { QUALITY_CONTROL                   } from './workflows/quality_control'
include { PIPELINE_INITIALISATION           } from './subworkflows/local/utils_nfcore_template_pipeline'
include { PIPELINE_COMPLETION               } from './subworkflows/local/utils_nfcore_template_pipeline'
include { ANNOTATION_CACHE_INITIALISATION   } from './subworkflows/local/annotation_cache_initialisation'
include { DOWNLOAD_CACHE_SNPEFF_VEP         } from './subworkflows/local/download_cache_snpeff_vep'
include { PREPARE_GENOME                    } from './subworkflows/local/prepare_genome'
include { PREPARE_INTERVALS                 } from './subworkflows/local/prepare_intervals'

// Initialize fasta file with meta map:
fasta = params.fasta ? Channel.fromPath(params.fasta).map{ it -> [ [id:it.baseName], it ] }.collect() : Channel.empty()

// Initialize file channels based on params, defined in the params.genomes[params.genome] scope
bcftools_annotations    = params.bcftools_annotations    ? Channel.fromPath(params.bcftools_annotations).collect()      : Channel.empty()
dbsnp                   = params.dbsnp                   ? Channel.fromPath(params.dbsnp).collect()                     : Channel.value([])
fasta_fai               = params.fasta_fai               ? Channel.fromPath(params.fasta_fai).collect()                 : Channel.empty()
germline_resource       = params.germline_resource       ? Channel.fromPath(params.germline_resource).collect()         : Channel.value([]) // Mutect2 does not require a germline resource, so set to optional input
known_indels            = params.known_indels            ? Channel.fromPath(params.known_indels).collect()              : Channel.value([])
known_snps              = params.known_snps              ? Channel.fromPath(params.known_snps).collect()                : Channel.value([])
pon                     = params.pon                     ? Channel.fromPath(params.pon).collect()                       : Channel.value([]) // PON is optional for Mutect2 (but highly recommended)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow QUALITY_CONTROL_PIPELINE {

    take:
    samplesheet // channel: samplesheet read in from --input

    main:
    versions = Channel.empty()

    // build indexes if needed
    // PREPARE_GENOME(
    //     params.ascat_alleles,
    //     params.ascat_loci,
    //     params.ascat_loci_gc,
    //     params.ascat_loci_rt,
    //     bcftools_annotations,
    //     params.chr_dir,
    //     dbsnp,
    //     fasta,
    //     germline_resource,
    //     known_indels,
    //     known_snps,
    //     pon)

    // fasta_fai   = params.fasta_fai  ? Channel.fromPath(params.fasta_fai).map{ it -> [ [id:'fai'], it ] }.collect()
    //                                 : PREPARE_GENOME.out.fasta_fai

    // Build intervals if needed
    // PREPARE_INTERVALS(fasta_fai, params.intervals, params.no_intervals, params.nucleotides_per_second, params.outdir, params.step)

        // Intervals for speed up preprocessing/variant calling by spread/gather
    // [interval.bed] all intervals in one file
    // intervals_bed_combined         = params.no_intervals ? Channel.value([]) : PREPARE_INTERVALS.out.intervals_bed_combined
    // intervals_bed_gz_tbi_combined  = params.no_intervals ? Channel.value([]) : PREPARE_INTERVALS.out.intervals_bed_gz_tbi_combined
    // intervals_bed_combined_for_variant_calling = PREPARE_INTERVALS.out.intervals_bed_combined

    // For QC during preprocessing, we don't need any intervals (MOSDEPTH doesn't take them for WGS)
    intervals_for_preprocessing = Channel.value([ [ id:'null' ], [] ])
    // intervals            = PREPARE_INTERVALS.out.intervals_bed        // [ interval, num_intervals ] multiple interval.bed files, divided by useful intervals for scatter/gather
    // intervals_bed_gz_tbi = PREPARE_INTERVALS.out.intervals_bed_gz_tbi // [ interval_bed, tbi, num_intervals ] multiple interval.bed.gz/.tbi files, divided by useful intervals for scatter/gather
    // intervals_and_num_intervals = intervals.map{ interval, num_intervals ->
    //     if ( num_intervals < 1 ) [ [], num_intervals ]
    //     else [ interval, num_intervals ]
    // }
    // intervals_bed_gz_tbi_and_num_intervals = intervals_bed_gz_tbi.map{ intervals, num_intervals ->
    //     if ( num_intervals < 1 ) [ [], [], num_intervals ]
    //     else [ intervals[0], intervals[1], num_intervals ]
    // }
    // if (params.tools && params.tools.split(',').contains('cnvkit')) {
    //     if (params.cnvkit_reference) {
    //         cnvkit_reference = Channel.fromPath(params.cnvkit_reference).collect()
    //     } else {
    //         PREPARE_REFERENCE_CNVKIT(fasta, intervals_bed_combined)
    //         cnvkit_reference = PREPARE_REFERENCE_CNVKIT.out.cnvkit_reference
    //         versions = versions.mix(PREPARE_REFERENCE_CNVKIT.out.versions)
    //     }
    // } else {
    //     cnvkit_reference = Channel.value([])
    // }
    // // Gather used softwares versions
    // versions = versions.mix(PREPARE_GENOME.out.versions)
    // versions = versions.mix(PREPARE_INTERVALS.out.versions)

    // vep_fasta = (params.vep_include_fasta) ? fasta.map{ fasta -> [ [ id:fasta.baseName ], fasta ] } : [[id: 'null'], []]

    // // Download cache
    // if (params.download_cache) {
    //     // Assuming that even if the cache is provided, if the user specify download_cache, sarek will download the cache
    //     ensemblvep_info = Channel.of([ [ id:"${params.vep_cache_version}_${params.vep_genome}" ], params.vep_genome, params.vep_species, params.vep_cache_version ])
    //     snpeff_info     = Channel.of([ [ id:"${params.snpeff_genome}.${params.snpeff_db}" ], params.snpeff_genome, params.snpeff_db ])
    //     DOWNLOAD_CACHE_SNPEFF_VEP(ensemblvep_info, snpeff_info)
    //     snpeff_cache = DOWNLOAD_CACHE_SNPEFF_VEP.out.snpeff_cache
    //     vep_cache    = DOWNLOAD_CACHE_SNPEFF_VEP.out.ensemblvep_cache.map{ meta, cache -> [ cache ] }

    //     versions = versions.mix(DOWNLOAD_CACHE_SNPEFF_VEP.out.versions)
    // } else {
    //     // Looks for cache information either locally or on the cloud
    //     ANNOTATION_CACHE_INITIALISATION(
    //         (params.snpeff_cache && params.tools && (params.tools.split(',').contains("snpeff") || params.tools.split(',').contains('merge'))),
    //         params.snpeff_cache,
    //         params.snpeff_genome,
    //         params.snpeff_db,
    //         (params.vep_cache && params.tools && (params.tools.split(',').contains("vep") || params.tools.split(',').contains('merge'))),
    //         params.vep_cache,
    //         params.vep_species,
    //         params.vep_cache_version,
    //         params.vep_genome,
    //         "Please refer to https://nf-co.re/sarek/docs/usage/#how-to-customise-snpeff-and-vep-annotation for more information.")

    //         snpeff_cache = ANNOTATION_CACHE_INITIALISATION.out.snpeff_cache
    //         vep_cache    = ANNOTATION_CACHE_INITIALISATION.out.ensemblvep_cache
    // }

    QUALITY_CONTROL (
        samplesheet,
        fasta,
        intervals_for_preprocessing)

    // emit:
    // multiqc_report = QUALITY_CONTROL_PIPELINE.out.multiqc_report // channel: /path/to/multiqc_report.html

}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:

    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.help,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // WORKFLOW: Run main workflow
    //
    QUALITY_CONTROL_PIPELINE (
        PIPELINE_INITIALISATION.out.samplesheet
    )

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.outdir,
        params.monochrome_logs,
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Get attribute from genome config file e.g. fasta
//

def getGenomeAttribute(attribute) {
    if (params.genomes && params.genome && params.genomes.containsKey(params.genome)) {
        if (params.genomes[ params.genome ].containsKey(attribute)) {
            return params.genomes[ params.genome ][ attribute ]
        }
    }
    return null
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
