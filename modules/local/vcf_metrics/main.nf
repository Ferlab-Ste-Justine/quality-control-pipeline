process VCF_METRICS {
    tag "${meta.id}"
    label 'process_single'

    input:
    tuple val(meta), val(metrics)

    output:
    tuple val(meta), path("*.vcf_metrics.json"), emit: json
    tuple val("${task.process}"), val('cat'), eval("cat --version | head -n 1 | cut -f4 -d' '"), topic: versions, emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def het_hom_ratio_snvs = metrics.count_hom_snvs > 0 ? metrics.count_het_snvs / metrics.count_hom_snvs : 0
    def het_hom_ratio_indels = metrics.count_hom_indels > 0 ? metrics.count_het_indels / metrics.count_hom_indels : 0
    def ins_dels_ratio = metrics.count_dels > 0 ? metrics.count_ins / metrics.count_dels : 0
    """
    cat <<-END_JSON > ${prefix}.vcf_metrics.json
    {
        "sample_id": "${meta.id}",
        "total_deletions": ${metrics.count_dels},
        "total_insertions": ${metrics.count_ins},
        "total_snvs": ${metrics.count_snvs},
        "heterozygous_snvs": ${metrics.count_het_snvs},
        "homozygous_snvs": ${metrics.count_hom_snvs},
        "heterozygous_indels": ${metrics.count_het_indels},
        "homozygous_indels": ${metrics.count_hom_indels},
        "het_hom_ratio_snvs": ${het_hom_ratio_snvs},
        "het_hom_ratio_indels": ${het_hom_ratio_indels},
        "ins_dels_ratio": ${ins_dels_ratio},
        "ti_tv_ratio": ${metrics.tstv_ratio}
    }
    END_JSON
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    cat > ${prefix}.vcf_metrics.json <<'END'
    {"sample_id":"stub"}
    END
    """
}
