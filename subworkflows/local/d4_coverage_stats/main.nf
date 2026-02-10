include { D4TOOLS_STAT as D4TOOLS_STAT_MEAN } from '../../../modules/local/d4tools/stat/main'
include { D4TOOLS_STAT as D4TOOLS_STAT_MEDIAN } from '../../../modules/local/d4tools/stat/main'
include { D4TOOLS_STAT as D4TOOLS_STAT_PCOV } from '../../../modules/local/d4tools/stat/main'
include { D4TOOLS_STAT as D4TOOLS_STAT_COUNT } from '../../../modules/local/d4tools/stat/main'
// include { MERGE_REGIONS_COV } from '../../../modules/local/merge_regions_cov/main'

workflow D4_COVERAGE_STATS {

    take:
    d4 // channel: [ val(meta), d4 ]
    ch_region_bed // [optional] [ val(meta2), bed] path to bed file with regions to calculate coverage stats over

    main:

    // mean - default
    D4TOOLS_STAT_MEAN(
        d4,
        ch_region_bed
    )

    // run with option --stat median
    D4TOOLS_STAT_MEDIAN(
        d4,
        ch_region_bed
    )

    // run with --stat perc_cov=5,15,20,50,100,200
    D4TOOLS_STAT_PCOV(
        d4,
        ch_region_bed
    )

    // run with --stat count
    D4TOOLS_STAT_COUNT(
        d4,
        ch_region_bed
    )

    ch_stats = D4TOOLS_STAT_MEAN.out.stat
        .join(D4TOOLS_STAT_MEDIAN.out.stat)
        .join(D4TOOLS_STAT_PCOV.out.stat)
        .join(D4TOOLS_STAT_COUNT.out.stat) // [ meta, mean, median, pcov, count ]

    // MERGE_REGIONS_COV(
    //     ch_stats
    // )

    emit:
    coverage_stats = ch_stats                     // channel: [ val(meta), path("*.stat.txt") ]

}
