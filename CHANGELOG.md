# Ferlab-Ste-Justine/quality-control-pipeline: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

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

### `Changed`
[#40](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/40) Refactoring QC report to support family vs cohort reports and fix some parsing errors.

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
