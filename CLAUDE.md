# CLAUDE.md

This file gives Claude Code the context it needs to work effectively in this repository.

## Project overview

`Ferlab-Ste-Justine/quality-control-pipeline` is a Nextflow DSL2 pipeline for comprehensive quality control of genomic sequencing data. It accepts FASTQ reads, BAM/CRAM alignments, and VCF/GVCF variant files — running a full suite of QC tools appropriate to each data type — and aggregates all results into a single (or per-family) MultiQC report. It can also ingest pre-computed DRAGEN per-sample metrics in place of recomputing them.

The repo is structured following nf-core conventions. It is _not_ a published nf-core pipeline.

Nextflow version range: `>=24.10.5, <26.0.0`. The floor is 24.10.5, not a rounder/lower number: this pipeline's version collection uses topic channels (`topic: versions` process outputs, `channel.topic("versions")` in `workflows/qualitycontrol.nf`), which don't exist before Nextflow 24.04.0 — 24.10.5 is simply the lowest version actually exercised by CI, so the declared floor is never again an untested claim (the previous floor, 23.10.1, was never actually runnable; nobody caught it until it was finally added to the CI matrix). Pipeline version is tracked in `nextflow.config` (`manifest.version`) and `.nf-core.yml` (`template.version`).

## High-level pipeline flow

The single entry workflow is in `main.nf`, which calls `QUALITYCONTROL` in `workflows/qualitycontrol.nf`. The samplesheet is branched by `meta.fileType` into `fastq` / `aln` (BAM or CRAM) / `vcf` (VCF or GVCF) channels. `--dragen_metrics_dir` gates two mutually exclusive modes:

1. **Standard mode** (default, no `--dragen_metrics_dir`):
   - `FASTQ_QC` — FastQC + ngsCheckMate (sample identity) on FASTQ rows.
   - `BAM_MERGE` — merges multi-lane/multi-run BAM/CRAM per sample, then `BAM_QC` (instantiated twice, as `BAM_QC_WGS` / `BAM_QC_TARGET`, split by `meta.sequencingType`) runs Samtools stats, Picard `CollectWgsMetrics`/`CollectHsMetrics`, quality-yield metrics, Mosdepth, per-gene coverage, and VerifyBamID2 contamination estimation.
   - `VCF_QC` (also instantiated as `VCF_QC_WGS` / `VCF_QC_TARGET`) — BCFtools variant counts by type/zygosity and Ts/Tv ratio on VCF/GVCF rows.
2. **DRAGEN mode** (`--dragen_metrics_dir` set): skips `BAM_MERGE`/`BAM_QC`/`VCF_QC` entirely. `SAMPLE_FAMILY_MAP` builds a sample→familyId map from the samplesheet; DRAGEN's own pre-computed per-sample CSV/BED metrics are picked up directly by glob from `--dragen_metrics_dir`, and `DRAGEN_COVERAGE_BY_GENE` aggregates DRAGEN's per-region coverage into the same per-gene coverage report the standard mode produces via Mosdepth.

Both modes always run `CRAM_SOMALIER` for pedigree/relatedness checking, either per-family or cohort-wide (`--cohort_mode`). The PED file comes from `--ped_file` if set, otherwise it's derived from samplesheet columns (`relationship_to_proband`, `affected_status`, `sex`) via `buildPedRowsForFamily`. `PED_FAMILY_CHECK` fails fast if a `--ped_file` family has no matching samples in the samplesheet. All QC outputs converge into `MULTIQC_PYTHON`, joined with PED data via `PED_MULTIQC_JOIN` — one report per family in per-family mode, one cohort-wide report in `--cohort_mode`.

The WGS/targeted distinction (`experimentalStrategy` column → `meta.sequencingType`) is load-bearing: `BAM_QC` and `VCF_QC` each run as two separate instances with different Picard metrics modules and bait/target BED inputs depending on it.

## Repository layout

```
main.nf                       # Entry point — calls PIPELINE_INITIALISATION, QUALITYCONTROL, PIPELINE_COMPLETION
nextflow.config                # Params, profiles, per-process resources, manifest
nextflow_schema.json           # Authoritative parameter schema (use this, not the README), grouped into $defs
nf-test.config                 # nf-test runner config (profile "test")
workflows/qualitycontrol.nf     # Main QUALITYCONTROL workflow — standard vs DRAGEN mode gating, somalier, MultiQC assembly
subworkflows/local/             # fastq_qc, bam_qc, vcf_qc, bam_merge, bam_id_repair, vcf_id_repair, cram_somalier,
                                #   sample_family_map, ped_family_check, ped_multiqc_join, qc_coverage_regions,
                                #   utils_nfcore_quality-control-pipeline_pipeline
subworkflows/nf-core/           # utils_nextflow_pipeline, utils_nfcore_pipeline, utils_nfschema_plugin,
                                #   fastq_ngscheckmate, vcf_extract_relate_somalier
modules/local/                  # samtools/{reheader,samples}, picard, d4tools, coverage_by_gene,
                                #   dragen_coverage_by_gene, vcf_metrics, multiqc_python
modules/nf-core/                # bcftools, fastp, fastqc, gatk4/bedtointervallist, mosdepth, multiqc, ngscheckmate,
                                #   picard, samtools, somalier, tabix, verifybamid
conf/                           # base.config, modules.config, test.config, test_dragen.config, test_full.config,
                                #   igenomes.config, modules/{bam_qc,fastq_qc,qc_coverage,somalier,vcf_qc}.config
assets/                         # samplesheet.csv, schema_input.json, qc_thresholds.yml, multiqc_config.yml,
                                #   BED files (capture/primary targets, autosomes non-gap regions, CDS canonical)
docs/                           # usage.md, output.md
```

## How to run

Typical invocation (from the README):

```bash
nextflow run Ferlab-Ste-Justine/quality-control-pipeline \
   -profile docker \
   --input samplesheet.csv \
   --outdir ./results \
   --fasta /path/to/reference.fa \
   --somalier_sites /path/to/sites.vcf.gz
```

Important conventions:

- Pass parameters via CLI flags or `-params-file` (JSON/YAML). **Do not** put params in a `-c` config file — `-c` is reserved for resource/infrastructure tuning. `docs/usage.md` and the README both call this out.
- `--dragen_metrics_dir` switches the pipeline into DRAGEN mode (see above) — mutually exclusive in practice with running the standard BAM_QC/VCF_QC path, though somalier still runs against any BAM/CRAM in the samplesheet either way.
- `--cohort_mode` (default `false`): per-family MultiQC reports and PED handling, vs. one cohort-wide report/PED.
- A given `sample` cannot provide both `BAM` and `CRAM` in the same run — validated at startup (`checkSingleAlignmentFileType` in `subworkflows/local/utils_nfcore_quality-control-pipeline_pipeline/main.nf`); mixing formats across _different_ samples is fine.

### Test dataset

The test data is expected to be accessible locally under the launch directory. Before testing the pipeline, verify that the `data-test/` directory exists.
The data lives in a private AWS S3 bucket: `s3://ferlab-public-dataset/nextflow/quality-control-pipeline/V1/data-test`.
In CI, `nf-test.yml` and `ci-full-run.yml` both download it through the `.github/actions/copy-test-data` composite action, the only place CI defines the S3 path. When the dataset version changes, update it there, in `scripts/run-smoke-tests.sh`, in `tests/nextflow.config` (`pipelines_testdata_base_path`, which nf-core lint requires but no test reads), and here.

### Quick smoke test

Unlike some sibling pipelines, a full `-stub` run does **not** work end-to-end here — only a handful of modules (`picard/collectqualityyieldmetrics`, `vcf_metrics`, `multiqc_python`, `d4tools/stat`, `coverage_by_gene`) define a `stub:` block, and running `-stub` against the rest of the pipeline fails partway through with an unrelated crash. For a real Docker-free sanity check that the pipeline actually launches and compiles under a given Nextflow version, use `-preview` instead — it fully executes/compiles the DSL2 script (including top-level statements) while skipping real task execution:

```bash
nextflow run . -profile test,docker -preview --outdir ./results_preview
```

### Test profile

If running locally:

- Make sure Docker Desktop is installed and running.

```bash
nextflow run Ferlab-Ste-Justine/quality-control-pipeline -profile test,docker
```

A separate `test_dragen` profile (`conf/test_dragen.config`) exercises DRAGEN mode against the same `data-test/` dataset:

```bash
nextflow run Ferlab-Ste-Justine/quality-control-pipeline -profile test_dragen,docker
```

To clean up outputs, run `nextflow clean -f`.

For the manual smoke-test routine (debug+test profile, plain test profile) as a single command, run `scripts/run-smoke-tests.sh` — it verifies `data-test/` is synced and Docker is running before starting, and leaves existing output directories in place rather than wiping them.

### Tests (nf-test)

`nf-test` is the testing framework. Config in `nf-test.config` sets `profile "test"` (add `docker` on the command line, e.g. `--profile test,docker`). nf-core upstream tests are excluded via the `ignore` glob.

If running locally:

- Make sure to run `export NXF_FILE_ROOT=$PWD` to allow nf-core's `nf-test` framework to find test files.

```bash
nf-test test                                  # run all tests
nf-test test tests/default.nf.test --tag pipeline --profile test,docker   # run the full end-to-end pipeline tests
nf-test test subworkflows/local/bam_qc        # target one module/subworkflow
```

The three full end-to-end pipeline tests are `tests/default.nf.test`, `tests/default_fam.nf.test`, and `tests/dragen.nf.test`. Test snapshots live alongside each module/subworkflow as `tests/main.nf.test.snap`.

To clean up test outputs, run `nf-test clean` or manually delete the `.nf-test/` directory.

`scripts/run-test-suite.sh` runs the full nf-test suite as a single pre-push gate, along with nf-core lint, pre-commit, an installed-nf-core-CLI-version check, and a launch check (via `-preview`) under the pipeline's declared minimum Nextflow version.

### Linting

CI workflows live in `.github/workflows/`: `linting.yml` (pre-commit + nf-core lint), `nf-test.yml` (sharded nf-test, using the composite actions in `.github/actions/`), `ci-full-run.yml` (full `-profile test`/`test_dragen` pipeline runs on Nextflow 24.10.5 and 25.10.4), and `ci-pr-title-lint.yml` (the PR title must look like `<type>: <TICKET-123> <description>`, e.g. `fix: BIOINFO-230 pin actions/checkout`; PRs are squash-merged, so the title becomes the commit message on `main`). Run lint locally with:

```bash
nf-core pipelines lint --release
```

`.nf-core.yml` carries lint overrides — several nf-core-template files are deliberately not present (e.g. `CODE_OF_CONDUCT.md`, nf-core logos, AWS CI workflows) because this is a Ferlab workflow, not a published nf-core pipeline. Don't reintroduce those files; instead update `.nf-core.yml` if you need to change lint behavior. As of this writing, `nf-core pipelines lint --release` still reports 14 failures — all pre-existing `check_local_copy` module/subworkflow version drift (deliberately deferred to a separate ticket, since some of it looks like intentional customization, e.g. container-registry choice) plus one known nf-core/tools 4.1.0 tool limitation (`multiqc_config: export_plots`, not actually gated by its ignore-list entry).

To format the files before committing, run:

```bash
pre-commit run --all-files
```

`nf-core pipelines lint --release` and `pre-commit run --all-files` are also bundled into `scripts/run-test-suite.sh`, alongside the full nf-test suite — duplicated on purpose with `.github/workflows/linting.yml` so formatting issues surface locally before CI does.

## Samplesheet format

See `assets/schema_input.json` for the authoritative schema; `docs/usage.md` has the full column reference and worked examples. Summary:

- Required: `participant`, `sample`, `fileType` (`FASTQ`/`BAM`/`CRAM`/`VCF`/`GVCF`), `file1`.
- Optional: `file2` (R2/index), `familyId` (required for `--ped_file` or per-family mode), `experimentalStrategy` (defaults `WGS`), `sex`, `status`, `relationship_to_proband`, `affected_status`, `lane`, `runId`.
- Rows sharing `participant` + `sample` + `experimentalStrategy` (e.g. multiple `lane`s) are merged before QC.
- A `sample` cannot mix `BAM` and `CRAM` fileTypes across its own rows (validated at startup); different samples may use different formats freely.

## Working in this codebase

A few patterns worth knowing before editing:

- **Channel shape convention.** Alignment channels carry `[meta, bam/cram, bai/crai]`; VCF channels carry `[meta, vcf, tbi]`. `meta` is trimmed of run-specific keys (`lane`, `runId`) once past per-run processing (e.g. right after the initial samplesheet branch in `qualitycontrol.nf`).
- **Startup validation** lives in `subworkflows/local/utils_nfcore_quality-control-pipeline_pipeline/main.nf`: `checkParticipantNotCorrupted` and `checkSingleAlignmentFileType` are plain `def` functions operating on parsed Lists (not Channels) — nf-test's `nextflow_function` test type can't `await` a Channel-returning function, which is why subworkflows that need to validate _channel_ data (e.g. `ped_family_check`, `sample_family_map`) stay as full `workflow`s instead.
- **Per-process resources** live in `conf/base.config`, gated by `check_max(...)` and the `max_cpus`/`max_memory`/`max_time` params.
- **DRAGEN mode is a genuinely different code path**, not a variation of the standard one — it skips `BAM_MERGE`/`BAM_QC`/`VCF_QC` entirely (see High-level pipeline flow above). When touching version-collection or MultiQC-input wiring, check both branches of the `if (params.dragen_metrics_dir)` gate in `qualitycontrol.nf`.
- **Version collection uses topic channels** (`topic: versions` process outputs, gathered via `channel.topic("versions")` in `qualitycontrol.nf`), gated by the `nextflow.preview.topic` flag in `main.nf` for Nextflow versions before it left preview. This is _not_ the classic nf-core `ch_versions.mix(...)` pattern used by the sibling pipelines — see "Project overview" above for why the Nextflow version floor is what it is.
- **Adding an nf-core module:** use `nf-core modules install <tool>` so `modules.json` stays consistent. Local-only logic goes under `modules/local/`. Several installed `modules/nf-core/**` modules have diverged from their recorded `modules.json` pin (deliberately, in at least some cases — e.g. container registry choice) without being recorded via `nf-core modules patch`; that's the pre-existing lint debt mentioned above, not something to "fix" by blindly running `nf-core modules update`.
- **Schema and params stay in sync.** `nextflow_schema.json` is the source of truth for parameter validation (via the `nf-schema` plugin). When adding a param, update both `nextflow.config` defaults and the schema.

## Outputs

`--outdir` is required. Per-sample QC files are published under `reports/QC/{sample}/`; per-family (or cohort-wide) pedigree/relatedness outputs under `reports/pedigree/`; merged/reheadered BAM/CRAM files under `results/{sequencingType}/`; the aggregated MultiQC report(s) under `multiqc/`. Pipeline run metadata, configs, timeline/report/trace/dag are written to `${outdir}/pipeline_info/`.

See `docs/output.md` for the full output layout.

## Pointers

- Parameter documentation: `nextflow_schema.json` (authoritative). README/usage docs intentionally avoid duplicating parameter details.
- Samplesheet reference: `docs/usage.md`.
- Output layout: `docs/output.md`.
- Changelog and version history: `CHANGELOG.md`.
- PR templates: `.github/PULL_REQUEST_TEMPLATE.md`.
