//
// Workflow that attaches the per-family ped file to the per-family MultiQC
// input, tolerating families present on only one side.
//
// meta and files are wrapped as one nested element (not two positional ones)
// before the join so that a right-only remainder (a ped family absent from
// the samplesheet) pads with a single null instead of one null per left
// field — with two positional fields, remainder:true's single-null padding
// leaves the downstream closure with too few arguments and crashes.
//
workflow PED_MULTIQC_JOIN {

    take:
    ch_multiqc_input // channel: [ meta, files ]  meta.id = familyId (or 'Cohort')
    ch_ped_by_id     // channel: [ id, ped ]

    main:
    multiqc_input_with_ped = ch_multiqc_input
        .map { meta, files -> [ meta.id, [meta, files] ] }
        .join(ch_ped_by_id, remainder: true)
        .filter { _id, meta_files, _ped -> meta_files != null }
        .map { _id, meta_files, ped -> [ meta_files[0], meta_files[1], ped ?: [] ] }

    emit:
    multiqc_input_with_ped
}
