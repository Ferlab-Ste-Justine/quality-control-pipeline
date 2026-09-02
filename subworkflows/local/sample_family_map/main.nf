//
// Workflow that builds a sample -> familyId map from the samplesheet, keyed by
// a String sample name. Used to attach familyId to files that are matched by
// sample name parsed from a filename (always a String), e.g. DRAGEN metrics
// files, so the join works even when meta.sample is an Integer (samplesheet
// supplied a numeric ID).
//
workflow SAMPLE_FAMILY_MAP {

    take:
    ch_samplesheet // channel: [ meta, files ]

    main:
    sample_family = ch_samplesheet
        .map { meta, _files -> [ meta.sample.toString(), meta.familyId ] }
        .unique()

    emit:
    sample_family // channel: [ sample:String, familyId ]
}
