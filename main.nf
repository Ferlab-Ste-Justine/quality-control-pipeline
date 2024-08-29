#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Ferlab-Ste-Justine/quality-control-pipeline
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

    intervals_for_preprocessing = Channel.value([ [ id:'null' ], [] ])

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
