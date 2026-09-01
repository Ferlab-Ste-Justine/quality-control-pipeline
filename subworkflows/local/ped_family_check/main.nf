//
// Workflow that fails fast when the ped file references a family that has no
// matching samples in the samplesheet. Left unchecked, such a family survives
// downstream as a dangling entry in the ped/MultiQC join and either crashes
// there or (once that join is hardened) silently drops the pedigree data with
// no explanation.
//
workflow PED_FAMILY_CHECK {

    take:
    ch_somalier_input_ped // channel: [ meta, ped ]     meta.id = ped family id
    ch_samplesheet        // channel: [ meta, files ]   meta.familyId = samplesheet family id

    main:
    ch_samplesheet_families = ch_samplesheet
        .map { meta, _files -> meta.familyId }
        .unique()

    ch_somalier_input_ped
        .map { meta, _ped -> meta.id }
        .collect()
        .map { it -> [it] }
        .combine(ch_samplesheet_families.collect().map { it -> [it] })
        .subscribe { pedFamilies, samplesheetFamilies ->
            def missing = pedFamilies - samplesheetFamilies
            if (missing) {
                error("Family ID(s) [${missing.join(', ')}] found in the pedigree but absent from the samplesheet's familyId column. Please check that the ped file's family_id matches the samplesheet.")
            }
        }

    emit:
    ch_somalier_input_ped
}
