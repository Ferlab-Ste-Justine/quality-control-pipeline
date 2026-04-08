# Ferlab-Ste-Justine/quality-control-pipeline: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased - [03/04/2026]

### `Changed`
[#29](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/29) Refactor coverage per gene workflow.

### `Deprecated`
[#29](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/pull/29) Deprecated d4tools workflow and modules.

## v1.0.0dev - [26/09/2025]

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
