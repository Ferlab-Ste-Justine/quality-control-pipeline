//
// Subworkflow with functionality specific to the Ferlab-Ste-Justine/quality-control-pipeline pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input samplesheet

    main:

    ch_versions = channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    //
    // Create channel from input file provided through params.input
    //

    ch_samplesheet = channel
        .fromList(samplesheetToList(input, "${projectDir}/assets/schema_input.json"))
        .map { meta, file1, file2 ->
            def fileType = inferFileTypeFromExtension(file1, meta.fileType)
            [ meta + [ participant_sample: "${meta.participant}_${meta.sample}", fileType: fileType ], [file1, file2] ]
        }
        .tap { ch_participant_sample } // save raw input channel
        .map { meta, files -> [meta.participant, meta.sequencingType, meta, files] }
        .groupTuple()
        .map { participant, seqtypes, metas, files ->
            def sequencingTypes = seqtypes.unique()
            [ sequencingTypes.size(), metas, files ]
        }
        .transpose()
        .map { n_sequencingTypes, meta, files -> [meta + [n_seqTypes:n_sequencingTypes ], files] }
        .map { meta, files -> [ meta - meta.subMap('lane'), meta.lane, files ] }
        .groupTuple() // group by meta
        .map { meta, lanes, files  ->  [meta + [n_lanes:lanes.size()], lanes.withIndex(), files] }
        .transpose()
        .map { meta, lane, files  ->  [meta + [lane:lane[0], lane_idx:lane[1]]] + files }
        .map { meta, file1, file2 ->
            if (meta.fileType == "FASTQ") {
                if (file2) {
                    assert (file2.name.endsWith('.fastq.gz') || file2.name.endsWith('.fq.gz')) : log.error("File 2 for sample ${meta.sample} does not have a valid FASTQ extension.")
                }
                def new_id = meta.n_seqTypes > 1 ? "${meta.sample}_${meta.sequencingType}_${meta.lane}" : "${meta.sample}_${meta.lane}"
                return [ meta + [ id:new_id, paired_end:file2 ? true : false ], file2 ? [ file1, file2 ] : [ file1 ] ]
            }
            else {
                if (!file2) {
                    def index_file = findIndex(meta.fileType, file1)
                    return [ meta, [file1, index_file]]
                }
            }
            [ meta, [file1, file2] ]
        }

    emit:
    samplesheet = ch_samplesheet
    versions    = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    multiqc_report  //  string: Path to MultiQC report

    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    //
    // Completion email and summary
    //
    workflow.onComplete {

        completionSummary(monochrome_logs)
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// Validate and Infer fileType from file extension
//
def inferFileTypeFromExtension(file, fileType=null) {
    def name = file.getFileName().toString() - '.gz'

    // Define mappings from file type to a list of possible extensions.
    def fileTypeMappings = [
        "GVCF" : ['.gvcf', '.g.vcf'],
        "VCF"  : ['.vcf'],
        "FASTQ": ['.fastq', '.fq'],
        "BAM"  : ['.bam'],
        "CRAM" : ['.cram']
    ]

    // Find the first map entry where any of its extensions match the end of the filename.
    def matchedEntry = fileTypeMappings.find { entry ->
        entry.value.any { extension -> name.endsWith(extension) }
    }

    if (matchedEntry) {
        // if fileType exists, check it matches inferred fileType - Validation
        if (fileType && matchedEntry?.key != fileType) {
            error("Inferred fileType '${matchedEntry.key}' from file extension does not match provided fileType '${fileType}' for file: ${name}. Please check the input samplesheet.")
        }
        return matchedEntry.key
    } else {
        log.warn("Unsupported fileType or file extension for file '${name}'. Supported fileTypes are: fastq, bam, cram, gvcf.")
    }
}

//
// Find index file for alignment or variant files
//
def findIndex(fileType, dataFile) {
    def index = dataFile.toString() + (fileType in ["BAM","CRAM"] ? (fileType == "BAM" ? '.bai' : '.crai') : '.tbi')
    if(!file(index).exists()) {
        log.debug("Index file not found for file: ${dataFile}. Expected index at: ${index}")
        return []
    }
    return file(index)
}

//
// Build PED rows for one family from samplesheet meta. Returns a list of
// [filename, tsv_line] tuples suitable for collectFile. In a prenatal trio
// (Mother/Father/Fetus, no Proband row) the Mother row is the proband biologically;
// Mother/Father still resolve to parents=0, and Fetus links to the family's
// Father/Mother rows — same as any other child relationship.
//
def buildPedRowsForFamily(familyId, metas, perFamily) {
    // Anchor parent/proband lookups to the primary sample row (sample_idx == 0),
    // whose samplename_somalier equals meta.sample.
    def fatherRow  = metas.find { m -> m.relationship_to_proband == 'Father'  && (m.sample_idx ?: 0) == 0 }
    def motherRow  = metas.find { m -> m.relationship_to_proband == 'Mother'  && (m.sample_idx ?: 0) == 0 }
    def probandRow = metas.find { m -> m.relationship_to_proband == 'Proband' && (m.sample_idx ?: 0) == 0 }

    def fatherId  = fatherRow?.samplename_somalier  ?: '0'
    def motherId  = motherRow?.samplename_somalier  ?: '0'
    def probandId = probandRow?.samplename_somalier ?: '0'
    def probandSex = probandRow?.sex

    def filename = perFamily ? "${familyId}.ped".toString() : 'Cohort.ped'

    return metas.collect { m ->
        def rel = m.relationship_to_proband ?: 'Proband'
        def parents = pedParentsFor(rel, fatherId, motherId, probandId, probandSex)
        def sex = sexForPed(m, rel)
        def phenotype = phenotypeForPed(m.affected_status)
        def line = [familyId, m.samplename_somalier, parents[0], parents[1], sex, phenotype].join('\t')
        [ filename, line ]
    }
}

//
// Map (relationship, family context) to [paternal_id, maternal_id].
//
def pedParentsFor(rel, fatherId, motherId, probandId, probandSex) {
    // These all link to the family's Father/Mother rows (proband, full siblings and twins
    // share both parents; the fetus links to its Mother/Father rows).
    def is_child = ['Proband', 'Fetus', 'Brother', 'Sister',
        'Identical twin', 'Fraternal twin brother', 'Fraternal twin sister']
    if (rel in is_child) {
        return [fatherId, motherId]
    }
    if (rel == 'Son' || rel == 'Daughter') {
        if (probandSex == 'Male')   { return [probandId, '0'] }
        if (probandSex == 'Female') { return ['0', probandId] }
        return ['0', '0']
    }
    // Father, Mother, Half-brother, Half-sister, Other, null
    return ['0', '0']
}

//
// PED sex code from samplesheet, warning if it disagrees with the relationship label.
//
def sexForPed(meta, rel) {
    def relationshipsBySex = [
        'Male'   : ['Father', 'Brother', 'Half-brother', 'Fraternal twin brother', 'Son'],
        'Female' : ['Mother', 'Sister', 'Half-sister', 'Fraternal twin sister', 'Daughter'],
    ]
    // Relationships not listed (Proband, Fetus, Identical twin, Other) imply no fixed sex.
    def expected = relationshipsBySex.find { entry -> rel in entry.value }?.key
    if (expected && meta.sex && meta.sex != 'NA' && meta.sex != expected) {
        log.warn("Sample ${meta.samplename_somalier}: relationship '${rel}' expects sex=${expected} but samplesheet has '${meta.sex}'. Using samplesheet sex.")
    }
    return meta.sex == 'Female' ? 2 : (meta.sex == 'Male' ? 1 : 0)
}

//
// PED phenotype code from affected_status.
//
def phenotypeForPed(status) {
    if (status == 'Affected')   { return 2 }
    if (status == 'Unaffected') { return 1 }
    return 0
}

//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    genomeExistsError()
}

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}
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

//
// Exit pipeline if incorrect --genome key provided
//
def genomeExistsError() {
    if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
        def error_string = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
            "  Genome '${params.genome}' not found in any config files provided to the pipeline.\n" +
            "  Currently, the available genome keys are:\n" +
            "  ${params.genomes.keySet().join(", ")}\n" +
            "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        error(error_string)
    }
}
//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "FastQC (Andrews 2010),",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}

