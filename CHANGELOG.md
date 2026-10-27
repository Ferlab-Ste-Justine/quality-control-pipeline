# Ferlab-Ste-Justine/quality-control-pipeline: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### `Added`

[#46](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/46) Reject a sample that provides both `BAM` and `CRAM` alignment inputs in the same run, with a clear error instead of a downstream MultiQC file-collision crash.
[#48](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/48) Added missing `meta.yml` files for local subworkflows and modules that lacked them, and rewrote stale ones that still described leftover scaffold subworkflows instead of the real pipeline.
[#48](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/48) Added `scripts/run-smoke-tests.sh` and `scripts/run-test-suite.sh`, local helper scripts for manual smoke-testing and running the full nf-test + lint suite before pushing.
[#49](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/49) Added `CLAUDE.md`, giving Claude Code repo context: pipeline flow, DRAGEN-mode gating, the topic-channel version-collection pattern, and known lint debt.
[#49](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/49) Added `.github/workflows/ci-full-run.yml` (full `-profile test`/`test_dragen` pipeline runs), `.github/workflows/ci-pr-title-lint.yml`, and a `.github/actions/copy-test-data` composite action.
[#50](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/50) Added `conf/containers_*.config`, the container-manifest files nf-core/tools 4.1.0's new `container_configs` lint check expects; regenerated automatically from module `meta.yml` files going forward.
[#51](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/51) Added a `dragenMetrics` samplesheet column: a per-sample directory or tarball of DRAGEN's pre-computed metrics. A tarball (`.tar`, `.tar.gz`/`.tgz`, `.tar.bz2`/`.tbz2`) is extracted via the newly-installed `untar` module (`modules/nf-core/untar`), one task per sample needing it, before the existing per-sample metric globbing applies — flat single-sample bundles and subdirectory trees inside the archive both work.
[#55](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/55) Extended the `--json-out` per-sample JSON sidecar (`multiqc_report.py`): added `schema_version`, `qc_metrics.cnv_metrics`, `pedigree.relatedness` (per-pair somalier comparisons), and a raw, unnormalized passthrough of each source tool's own native fields under `dragen_metrics.{mapping,wgs_coverage,vc,ploidy_estimation,cnv,gvcf,capture_coverage,time}_metrics` and `picard_metrics.{HsMetrics,WgsMetrics,QualityYieldMetrics}` — intended to let downstream consumers (e.g. CQDG) eventually source this JSON instead of generating their own equivalent report. Also added parsing for three DRAGEN file types not previously staged at all: `*.gvcf_metrics.csv`, `*.qc-coverage-region-1_coverage_metrics.csv` (`capture_coverage_metrics`), and `*.time_metrics.csv`.

### `Removed`

[#50](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/50) Removed the `fastp` and `samtools/idxstats` modules: both installed but never invoked anywhere in the pipeline.

### `Fixed`

[#46](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/46) Accept integer sample IDs, and (for `participant`) decimal-number IDs, in the input samplesheet schema, without silently corrupting IDs with a leading zero before a decimal point.
[#46](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/46) Fixed a DRAGEN-mode crash when a sample's numeric ID didn't type-match the sample name parsed from DRAGEN metric filenames.
[#46](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/46) Fixed a pipeline crash — and, separately, silently dropped pedigree data — when a `--ped_file` family had no matching samples in the samplesheet; now fails fast with a clear error instead.
[#46](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/46) Fixed the bundled example `assets/samplesheet.csv`, which used an outdated column format that no longer matched the current schema.
[#46](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/46) Fixed a flaky nf-test snapshot comparison caused by FastQC embedding the current date in its HTML report.
[#46](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/46) Fixed a false-positive validation error (and a mangled BAM/CRAM reheader sample name) for a decimal `participant` ID with a magnitude small enough that `BigDecimal.toString()` renders it in scientific notation.
[#48](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/48) Fixed a typo in `.nf-core.yml` (`include_configs` → `included_configs`) that had silently disabled that lint check.
[#48](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/48) Fixed `linting.yml`'s `master`/`main` branch-name check.
[#48](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/48) Fixed the pre-commit CI job's prettier/editorconfig failures (mostly YAML list-indentation style), and added a `[LICENSE]` exemption to `.editorconfig` instead of reflowing the verbatim Apache License 2.0 text.
[#48](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/48) Raised the declared `nextflowVersion` floor from `23.10.1` to `24.10.5`: the pipeline's use of topic channels for version collection (see `main.nf`) requires Nextflow >=24.04.0, so `23.10.1` was never actually runnable, just untested until this PR added it to the CI matrix.
[#49](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/49) Fixed a stale/typo'd CI reference in `.github/CONTRIBUTING.md` (`.github/workflow/ci.yml` → `.github/workflows/ci-full-run.yml`).
[#50](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/50) Fixed the remaining `check_local_copy` lint failures (module/subworkflow version drift deferred since #48) by updating `bcftools/{index,query,reheader}`, `samtools/{index,merge,stats}`, `fastqc`, `multiqc`, and `fastq_ngscheckmate` (+ its two child modules) to the newest versions compatible with our declared Nextflow floor, re-pointing containers to biocontainers/depot.galaxyproject.org and recording the substitution via `nf-core modules patch`. `nf-core pipelines lint --release` goes from 14 failures to 1 (the already-documented `multiqc_config: export_plots` nf-core/tools bug).
[#50](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/50) Fixed two real, pre-existing bugs these module updates surfaced: `bam_merge`'s `SAMTOOLS_MERGE` call was passing `fasta`/`fai` as two separate tuples instead of the one combined tuple every upstream version expects, and `bam_qc`'s `SAMTOOLS_STATS` call wasn't passing `fai` at all despite the module requiring it — both latent since the local module content had silently diverged from its recorded `modules.json` pin long before this PR.
[#51](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/51) Fixed a silent misattribution risk in DRAGEN metrics ingestion: a sample's `dragenMetrics` directory/tarball is expected to hold only that sample's own files, but nothing previously checked this — any file found there was unconditionally attributed to that samplesheet row. A mismatched `dragenMetrics` value (e.g. two samples pointed at the same bundle by mistake) would have silently attached another sample's QC data with no error. Added `dragenFileMatchesSample()`, which compares each matched file's filename-embedded sample ID against `meta.sample` and drops (with a `log.warn`) anything that doesn't match.
[#52](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/52) Fixed the custom "Quality Control" MultiQC tables (General Status, Alignment Metrics, Variant Calling Metrics, Pedigree) getting a random `pconfig` id on every run, which broke any downstream consumer (e.g. Radiant's ETL) trying to find a given table by name in the exported `<familyId>_multiqc_report_data.zip`. Each table now gets a fixed id (`qc-general-status`, `qc-alignment-metrics`, `qc-variant-calling-metrics`, `qc-pedigree`), documented in `docs/output.md` as a stable contract for consumers.
[#53](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/53) Fixed DRAGEN mode never surfacing CNV metrics in the MultiQC report: `*.cnv_metrics.csv` was already staged into `MULTIQC_PYTHON` by the generic DRAGEN `*.csv` glob, but `multiqc_report.py` never parsed it. Added a "CNV Metrics" table (amplification/deletion counts, passing counts, segment count, coverage MAD — from DRAGEN's `CNV SUMMARY` section), published with the fixed id `qc-cnv-metrics` (see #52). A sample with no `cnv_metrics.csv` is simply omitted from this table, matching how every other table in this report already handles a sample with no data for that category.
[#54](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/54) Fixed DRAGEN mode never surfacing exome capture-enrichment metrics: DRAGEN's own bundled Picard `*.CollectHsMetrics.out` wasn't staged into `MULTIQC_PYTHON`, so General Status/Alignment Metrics showed no coverage for DRAGEN-sourced exome samples. MultiQC's `picard` module already content-detects and parses this file natively — staging it (added to the existing `ch_dragen_files` glob) is the whole fix. Also fixed `assets/multiqc_config.yml`'s `extra_fn_clean_exts`, which had an entry for our own pipeline's HsMetrics filename but not DRAGEN's, so the sample name was resolving uncleaned (e.g. `NA19238.CollectHsMetrics` instead of `NA19238`). Separately found and worked around a genuine MultiQC 1.30 bug: a malformed Picard metrics file crashes its `HsMetrics` parser with `dictionary changed size during iteration`, silently taking down the _entire_ picard module for that run — including other, unrelated samples' otherwise-valid data in the same family report. Excluded `tests/fixtures/**/*.out` from pre-commit's `trailing-whitespace` hook, which had been silently corrupting the test fixture's column count and triggering exactly this crash.
[#55](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/55) Fixed `tests/dragen.nf.test` silently running under the wrong Nextflow profile in both `scripts/run-test-suite.sh` and CI (`.github/workflows/nf-test.yml`): nf-test's `--profile` flag applies to every test an invocation runs, not per test file, so a single blanket `test,docker` call left this test's `params.ped_file` set (from `conf/test.config`'s default) instead of unset — flipping `workflows/qualitycontrol.nf` onto the provided-PED-file code path instead of the derive-from-samplesheet path this test actually needs, silently dropping derived `.ped` output. The test now declares `profile "test_dragen"` itself, so it gets the right profile in any invocation that passes extra profiles with a `+` prefix (CI's `--profile=+docker`, and `scripts/run-test-suite.sh`'s single `--profile +docker` run). Splitting it into its own CI job by test path doesn't work: `nf-test --changed-since` replaces any paths given with the changed files (and a changed `triggers` file selects every test), so each job ended up running the other's tests under the wrong profile.
[#55](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/55) Fixed a data-correctness bug in the new `dragen_metrics` raw passthrough: several DRAGEN CSVs carry multiple sections that aren't safe to merge blindly — `vc_metrics.csv`'s PREFILTER/POSTFILTER/SUMMARY sections can report different values under the same metric name, `mapping_metrics.csv`'s per-read-group section repeats SUMMARY's own metric names, and `cnv_metrics.csv`'s `SEX GENOTYPER` section lists every reference-panel sample's own inferred sex keyed by that sample's ID. The raw dump previously flattened every section in the file with no collision protection; it's now restricted to the single canonical section the curated `SampleMetrics` fields already read from.
[#55](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/55) Fixed non-deterministic JSON sidecar content: `pedigree.relatedness`'s array order depended on somalier/upstream channel ordering, which isn't guaranteed stable between runs, and `dragen_metrics`/`picard_metrics`'s key order depended on file-staging order — both now use a fixed, deterministic order.
[#55](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/55) Fixed CI breaking when GitHub moved `ubuntu-latest` to Ubuntu 26.04: the `nf-test-changes` job ran nf-test 0.9.2 on the image's default Java 25, which nf-test can't parse (`Unsupported class file major version 69`). All workflows now pin `runs-on: ubuntu-24.04`, and `nf-test-changes` installs Java 17 (Zulu) explicitly. Also fixed `confirm-pass` passing when `nf-test-changes` failed (the skipped test matrix didn't count as a failure); it now depends on both jobs.

### `Changed`

[#47](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/47) Updated documentation related to dragen-metrics mode and usage
[#48](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/48) Bumped the pinned `nf-core/tools` lint version to 4.1.0 and resolved the resulting lint failures/warnings, including regrouping `nextflow_schema.json` into `$defs` sections and updating `.prettierrc.yml` to match the current template.
[#48](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/48) Updated the README Nextflow badge and the nf-test CI matrix's pinned Nextflow versions to match the corrected `nextflowVersion` floor.
[#49](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/49) Harmonized `.github/` CI conventions with the sibling Post-processing-Pipeline: pinned `nf-core/setup-nextflow` (previously unpinned), replaced `commit_lint.yml` with a PR-title lint (PRs are squash-merged, so the PR title becomes the commit message on `main`), and added `run_all`/disk-cleanup support to `nf-test.yml`. `nf-test.config`'s `nft-utils` plugin deliberately stays at `1.0.0` rather than the sibling's `0.0.3`, which lacks a `removeFromYamlMap` overload this pipeline's tests rely on.
[#51](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/51) Replaced the pipeline-wide `--dragen_metrics_dir` parameter with the per-sample `dragenMetrics` samplesheet column — DRAGEN sourcing is no longer an all-or-nothing whole-run mode; standard-recompute and DRAGEN-sourced samples can now be mixed in the same run, including within the same family. `qualitycontrol.nf`'s fastq/aln/vcf channels are each branched again on `meta.dragenMetrics`, independently per row, instead of gating the whole `FASTQ_QC`/`BAM_MERGE`/`BAM_QC`/`VCF_QC` block behind a single pipeline flag. Existing invocations/config files using `--dragen_metrics_dir`/`--dragen_metrics` need updating to the samplesheet column instead.
[#49](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/49) Bumped prettier to 3.6.2 and added `trailing-whitespace`/`end-of-file-fixer` pre-commit hooks, matching the fleet.
[#50](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/50) `multiqc` deliberately stays at 1.27 rather than upstream's 1.31 — it's driven through a custom script calling MultiQC's Python API directly, so a version bump is a real compatibility risk, not a drive-by container swap. Patched to resolve its lint failure instead of bumping.
[#50](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/50) `samtools/merge` deliberately stays at 1.23.1 rather than upstream's latest 1.24, which introduces a variable-scoping pattern that fails to compile under our declared Nextflow floor (24.10.5).

## v2.0.0 - [11/06/2026]

### `Fixed`

[#35](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/35) Fixed issue with MultiQC report title argument.
[#34](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/34) Fixed issue with ped file channel initialization.
[#36](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/36) Fixed using MultiQC title in report.
[#37](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/37) Fixed issue with somalier input ped file generation.

### `Added`

[#36](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/36) Added test profile
[#38](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/38) Added new samplesheet fields: `relationship_to_proband` and `affected_status` and feature to generate pedigree file from samplesheet.
[#39](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/39) Added `cohort_mode`. This will allow the reports to be one per family (or sample if a solo), or per cohort/batch.
[#41](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/41) Added `--dragen_metrics_dir`: when set, BAM_QC and VCF_QC are skipped and the report's alignment / variant metrics come from DRAGEN's pre-computed CSVs. Somalier still runs against any BAM/CRAM provided in the samplesheet for pedigree validation.
[#42](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/42) Added `DRAGEN_COVERAGE_BY_GENE` module to aggregate dragen region metrics to `coverage_by_gene.tsv`.

### `Changed`

[#40](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/40) Refactoring QC report to support family vs cohort reports and fix some parsing errors.
[#44](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/44) Updated targeted sequencing stream to take baits and targets.
[#43](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/43) Updated nf-tests and added snapshots.

### `Deprecated`

[#39](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/39) Deprecated `somalier_perfamily`

## v0.0.3dev - [10/04/2026]

### `Fixed`

[#33](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/33) Fixed sample name parsing in picard and verifybamID multiqc modules. Removed default of 0 when value is None.

## v0.0.2dev - [10/04/2026]

### `Added`

[#30](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/30) Added pipeline documentation and refactor configs.

### `Changed`

[#29](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/29) Refactor coverage per gene workflow.
[#32](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/32) Use biocontainers picard instead of broadinstitue's picard docker image.

### `Deprecated`

[#29](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/29) Deprecated d4tools workflow and modules.

### `Fixed`

[#31](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/31) Fixed bug when processing without specified coverage regions.

## v0.0.1dev - [26/09/2025]

### `Added`

- [#28](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/28) Added first version of multiqc report and gathering of VCF stats.
- [#26](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/26) Added workflow to perform post-variant calling quality control according to the GA4GH QC standards.
- [#21](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/21) Added option to calculate coverage metrics on multiple BED files using d4tools
- [#20](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/20) Added FastQ validation and QC workflow
- [#19](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/19) Added somalier workflow and modules
- [#18](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/18) Added bam/cram merge workflow
- [#16](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/16) Added bam/cram QC workflow and modules
- [#12](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/12) Added vcf id repair subworkflow
- [#11](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/11) Added nf-core modules

### `Fixed`

- [#13](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/13) Fixed naming of vcf id repair wf
- [#10](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/10) Cleaned-up template

### `Dependencies`

- [#15](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/15) Edited utils_nfcore - Parsing samplesheet input into 3 channels
- [#14](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/14) Added CODEOWNERS, Input schema and edited configs

### `Deprecated`
