# Ferlab-Ste-Justine/quality-control-pipeline: Usage

## Introduction

This pipeline performs quality control on genomic sequencing data (FASTQ, BAM/CRAM, VCF). It accepts a CSV samplesheet listing samples and their associated files, runs QC tools appropriate to each data type, and produces per-sample output directories alongside a single aggregated MultiQC report.

## Samplesheet

Prepare a CSV samplesheet with one row per file (or per sequencing run for multi-run samples). Pass it to the pipeline with `--input`.

```bash
--input '[path to samplesheet.csv]'
```

### Column reference

| Column                 | Required | Description                                                                                                                                                                               |
| ---------------------- | -------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `participant`          | Yes      | Participant/individual identifier. Used to group samples from the same person across sequencing types. May be a plain string, an integer, or (in a JSON/YAML samplesheet) a decimal number — avoid a leading zero directly before a decimal point (e.g. `01.5`), which validation rejects rather than silently truncate. |
| `sample`               | Yes      | Sample identifier. Used as the output directory name (`reports/QC/{sample}`). Must be unique per participant + sequencing type combination. May be a plain string or an integer.           |
| `fileType`             | Yes      | Type of input file. One of: `FASTQ`, `BAM`, `CRAM`, `VCF`, `GVCF` (case-insensitive). A given `sample` cannot provide both `BAM` and `CRAM` in the same run — pick one alignment format per sample. |
| `file1`                | Yes      | Path to the primary file: FASTQ R1, BAM, CRAM, or VCF/GVCF.                                                                                                                               |
| `file2`                | No       | Path to the secondary file: FASTQ R2, BAI index, CRAI index, or VCF TBI/CSI index.                                                                                                        |
| `familyId`             | No       | Family identifier. Required when running per-family (`--cohort_mode false`) or with `--ped_file`. May be a plain string or an integer. When using `--ped_file`, every family referenced in the ped file must also appear in the samplesheet's `familyId` column, or the pipeline fails fast at startup. |
| `experimentalStrategy` | No       | Sequencing strategy. One of: `WGS`, `WXS`, `TARS`, `RNAS`, `ATACS`, `BIS`, `TMS`, `CHIPS`. Defaults to `WGS`. Controls which QC modules run (e.g. Picard WGS metrics only run for `WGS`). |
| `sex`                  | No       | Sample sex. One of: `Female`, `Male`, `Other`, `NA`. Defaults to `NA`.                                                                                                                    |
| `status`               | No       | Sample status. `0` = normal, `1` = tumor. Defaults to `0`.                                                                                                                                |
| `relationship_to_proband` | No    | Relationship of the sample to the family proband. One of: `Proband`, `Mother`, `Father`, `Brother`, `Sister`, `Half-brother`, `Half-sister`, `Identical twin`, `Fraternal twin brother`, `Fraternal twin sister`, `Son`, `Daughter`, `Other`, `Fetus`. Defaults to `Proband`. Used when generating a pedigree from the samplesheet. |
| `affected_status`      | No       | Affected status for pedigree generation. One of: `Affected`, `Unaffected`, `Unknown`. Defaults to `Unknown`.                                                                              |
| `lane`                 | No       | Lane or run identifier. When multiple rows share the same `participant`, `sample`, and `experimentalStrategy`, they are merged before QC.                                                 |
| `runId`                | No       | Run identifier. Used together with `lane` to uniquely identify a sequencing run.                                                                                                          |

### Example samplesheet

```csv
participant,sample,familyId,experimentalStrategy,sex,lane,fileType,file1,file2
P001,S001,FAM1,WGS,Female,L001,FASTQ,/data/S001_L001_R1.fastq.gz,/data/S001_L001_R2.fastq.gz
P001,S001,FAM1,WGS,Female,L002,FASTQ,/data/S001_L002_R1.fastq.gz,/data/S001_L002_R2.fastq.gz
P002,S002,FAM1,WGS,Male,,BAM,/data/S002.bam,/data/S002.bam.bai
P003,S003,FAM2,WGS,Female,,CRAM,/data/S003.cram,/data/S003.cram.crai
P004,S004,FAM2,WGS,Male,,VCF,/data/S004.vcf.gz,/data/S004.vcf.gz.tbi
```

> `S001` above has two rows for two sequencing lanes — they will be merged before alignment QC.

## Running the pipeline

### Minimal command

```bash
nextflow run Ferlab-Ste-Justine/quality-control-pipeline \
   -profile docker \
   --input samplesheet.csv \
   --outdir ./results
```

### Recommended command (with reference files)

```bash
nextflow run Ferlab-Ste-Justine/quality-control-pipeline \
   -profile docker \
   --input samplesheet.csv \
   --outdir ./results \
   --fasta /path/to/GRCh38.fa \
   --fai /path/to/GRCh38.fa.fai \
   --fasta_dict /path/to/GRCh38.dict \
   --somalier_sites /path/to/sites.hg38.vcf.gz \
   --ped_file /path/to/cohort.ped \
   --cohort_mode false \
   --verifybamid_svd_prefix /path/to/1000g.phase3 \
   --ngscheckmate_snp_pt /path/to/SNP_GRCh38_hg38_wChr.bed
```

### Using a params file

Pipeline parameters can be specified in a YAML or JSON file via `-params-file`:

```bash
nextflow run Ferlab-Ste-Justine/quality-control-pipeline \
   -profile docker \
   -params-file params.yaml
```

```yaml title="params.yaml"
input: ./samplesheet.csv
outdir: ./results
fasta: /path/to/GRCh38.fa
fai: /path/to/GRCh38.fa.fai
fasta_dict: /path/to/GRCh38.dict
somalier_sites: /path/to/sites.hg38.vcf.gz
ped_file: /path/to/cohort.ped
cohort_mode: false
verifybamid_svd_prefix: /path/to/1000g.phase3
ngscheckmate_snp_pt: /path/to/SNP_GRCh38_hg38_wChr.bed
```

> [!WARNING]
> Do not use `-c <file>` to specify pipeline parameters — custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) and infrastructural settings.

## Pipeline parameters

### Input / output

| Parameter         | Description                                                                  |
| ----------------- | ---------------------------------------------------------------------------- |
| `--input`         | Path to the input samplesheet CSV.                                           |
| `--outdir`        | Directory where results will be saved. Use absolute paths for cloud storage. |
| `--multiqc_title` | Title string printed in the MultiQC report header.                           |

### Reference files

| Parameter       | Description                                                                         |
| --------------- | ----------------------------------------------------------------------------------- |
| `--fasta`       | Path to the reference genome FASTA file. Required for alignment QC and Somalier.    |
| `--fai`         | Path to the FASTA index (`.fai`).                                                   |
| `--fasta_dict`  | Path to the sequence dictionary (`.dict`). Required for Picard.                     |
| `--regions_bed` | BED file of capture or analysis regions. Used by GATK BedToIntervalList and Picard. |

### Sample identity & relatedness (Somalier)

| Parameter          | Description                                                                                                                                                                                                                                                                                                                                                                                                                                               |
| ------------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--somalier_sites` | VCF of known variant sites used by Somalier for fingerprinting. Refer to the [brentp/somalier github release page](https://github.com/brentp/somalier/releases#release-v0.2.19). The release corresponding to the container version (v0.2.19) is provided, however a more recent version may also be used. Choose the reference that maps your sample's genome build and ensure samples within the same `SOMALIER_RELATE` run uses the somalier reference |
| `--ped_file`       | Pedigree file (PED format) describing family relationships. Optional — when omitted, a pedigree is generated automatically from the samplesheet (`familyId`, `sex`, `relationship_to_proband`, `affected_status`).                                                                                                                                                                                                                                        |
| `--cohort_mode`    | Boolean. When `true` the pipeline runs as a single cohort: one Somalier relate across all samples, one MultiQC report. When `false` (default), samples are grouped by `familyId` — Somalier runs once per family and one MultiQC report is produced per family.                                                                                                                                                                                           |

### Contamination (VerifyBamID2)

| Parameter                  | Description                                                                   |
| -------------------------- | ----------------------------------------------------------------------------- |
| `--verifybamid_svd_prefix` | File path prefix for VerifyBamID2 SVD reference files (`.UD`, `.mu`, `.bed`). |

### Sample identity (ngsCheckMate)

| Parameter               | Description                                                        |
| ----------------------- | ------------------------------------------------------------------ |
| `--ngscheckmate_snp_pt` | SNP panel BED file for ngsCheckMate FASTQ-level identity checking. |

### Coverage QC regions

Up to two custom BED region sets can be provided for per-region coverage analysis (e.g. a gene panel BED and a clinically relevant regions BED):

| Parameter                | Description                                                                            |
| ------------------------ | -------------------------------------------------------------------------------------- |
| `--qc_coverage_region_1` | BED file for the first QC coverage region set.                                         |
| `--region_1_name`        | Label for the first region set, used in output file names. Defaults to `qc_regions_1`. |
| `--qc_coverage_region_2` | BED file for the second QC coverage region set.                                        |
| `--region_2_name`        | Label for the second region set. Defaults to `qc_regions_2`.                           |

### VCF QC

| Parameter       | Description                                                 |
| --------------- | ----------------------------------------------------------- |
| `--targets_bed` | BED file of target regions for filtering VCF QC statistics. |
| `--exons_bed`   | BED file of exon regions, passed to `bcftools stats`.       |

### DRAGEN metrics input

When samples have already been processed by DRAGEN, the pipeline can build the report from DRAGEN's per-sample metric CSVs instead of recomputing them with BAM_QC / VCF_QC.

| Parameter              | Description
| ---------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--dragen_metrics_dir` | Directory containing DRAGEN per-sample metric CSVs. When set, BAM_QC and VCF_QC are skipped. Somalier still runs against any BAM/CRAM in the samplesheet for pedigree validation. |

The directory is globbed at any depth for these filenames (both `<sample>.<type>.csv` and `<sample>.final.<type>.csv` are recognised, therefore ensure the file prefix corresponds to the sample ID):

| Filename suffix                              | Skips                               | Used for                                                                                                                                                                                                                                                           |
| -------------------------------------------- | ----------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| `*.mapping_metrics.csv`                      | FASTQ_QC, BAM_MERGE, BAM_QC, VCF_QC | `MULTIQC_PYTHON`: Alignment metrics + Q30 yield + estimated contamination                                                                                                                                                                                                            |
| `*.wgs_coverage_metrics.csv`                 | FASTQ_QC, BAM_MERGE, BAM_QC, VCF_QC | `MULTIQC_PYTHON`: Mean autosome coverage, uniformity (MAD proxy), %15x                                                                                                                                                                                                               |
| `*.vc_metrics.csv`                           | FASTQ_QC, BAM_MERGE, BAM_QC, VCF_QC | `MULTIQC_PYTHON`: SNV / insertion / deletion counts, Ti/Tv, het:hom ratios                                                                                                                                                                                                           |
| `*.ploidy_estimation_metrics.csv`            |                                     | `CRAM_SOMALIER` (optional): XX / XY ploidy → predicted-sex fallback for `sex_check` when somalier didn't run for that sample                                                                                                                                                                   |
| `*_cov_report.bed` + `*_read_cov_report.bed` |                                     | `DRAGEN_COVERAGE_BY_GENE` (optional): Paired per region (`qc-coverage-region-[1,2,3]`); the `DRAGEN_COVERAGE_BY_GENE` module aggregates them into a `*.coverage_by_gene.tsv`. The per-gene report tab uses region 1 by default and falls back to region 2 when region 1 has no gene-annotated intervals. |

`--dragen_metrics_dir` accepts local paths and remote URIs (`s3://`, `gs://`, `az://`) provided the appropriate Nextflow plugin/credentials are configured. Files whose sample prefix doesn't match a row in the samplesheet are silently ignored.

Example — DRAGEN-only run with a small mixed samplesheet (some BAM for somalier, the rest GVCF-only):

```bash
nextflow run Ferlab-Ste-Justine/quality-control-pipeline \
   -profile docker \
   --input samplesheet.csv \
   --outdir ./results \
   --fasta /path/to/GRCh38.fa \
   --fai /path/to/GRCh38.fa.fai \
   --somalier_sites /path/to/sites.hg38.vcf.gz \
   --dragen_metrics_dir s3://my-bucket/dragen-outputs/
```

The following organization for the DRAGEN directory `s3://my-bucket/dragen-outputs/` are valid:

<details> <summary>DRAGEN directory organized by samples</summary>

```bash
my-bucket/dragen-outputs/
├── S001/
│   ├── S001(\.final)?.mapping_metrics.csv
│   ├── S001(\.final)?.wgs_coverage_metrics.csv
│   ├── S001(\.final)?.vc_metrics.csv
│   ├── S001(\.final)?.ploidy_estimation_metrics.csv
│   ├── S001_cov_report.bed
│   └── S001_read_cov_report.bed
├── S002/
│   ├── S002(\.final)?.mapping_metrics.csv
│   ├── S002(\.final)?.wgs_coverage_metrics.csv
│   ├── S002(\.final)?.vc_metrics.csv
│   ├── S002(\.final)?.ploidy_estimation_metrics.csv
│   ├── S002_cov_report.bed
│   └── S002_read_cov_report.bed
└── S003/
    ├── S003(\.final)?.mapping_metrics.csv
    ├── S003(\.final)?.wgs_coverage_metrics.csv
    ├── S003(\.final)?.vc_metrics.csv
    ├── S003(\.final)?.ploidy_estimation_metrics.csv
    ├── S003_cov_report.bed
    └── S003_read_cov_report.bed
```

</details>

<details> <summary>DRAGEN directory organized by runs</summary>

```bash
my-bucket/dragen_outputs/
├── Run001/
│   ├── S001/
│   │   ├── S001(\.final)?.mapping_metrics.csv
│   │   ├── S001(\.final)?.wgs_coverage_metrics.csv
│   │   ├── S001(\.final)?.vc_metrics.csv
│   │   ├── S001(\.final)?.ploidy_estimation_metrics.csv
│   │   ├── S001_cov_report.bed
│   │   └── S001_read_cov_report.bed
│   └── S002/
│       ├── S002(\.final)?.mapping_metrics.csv
│       ├── S002(\.final)?.wgs_coverage_metrics.csv
│       ├── S002(\.final)?.vc_metrics.csv
│       ├── S002(\.final)?.ploidy_estimation_metrics.csv
│       ├── S002_cov_report.bed
│       └── S002_read_cov_report.bed
└── Run002/
    ├── S003/
    │   ├── S003(\.final)?.mapping_metrics.csv
    │   ├── S003(\.final)?.wgs_coverage_metrics.csv
    │   ├── S003(\.final)?.vc_metrics.csv
    │   ├── S003(\.final)?.ploidy_estimation_metrics.csv
    │   ├── S003_cov_report.bed
    │   └── S003_read_cov_report.bed
    └── S004/
        ├── S004(\.final)?.mapping_metrics.csv
        ├── S004(\.final)?.wgs_coverage_metrics.csv
        ├── S004(\.final)?.vc_metrics.csv
        ├── S004(\.final)?.ploidy_estimation_metrics.csv
        ├── S004_cov_report.bed
        └── S004_read_cov_report.bed
```

</details>


### Pipeline behaviour

| Parameter         | Description                                                                                                      |
| ----------------- | ---------------------------------------------------------------------------------------------------------------- |
| `--skip_merge`    | Skip BAM/CRAM merging even when a sample has multiple lanes. Default: `true`.                                    |
| `--skip_reheader` | Skip reheadering BAM/CRAM files when the read group sample name does not match the samplesheet. Default: `true`. |

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Several generic profiles are bundled with the pipeline to configure the container engine:

- `docker` — run with [Docker](https://docker.com/)
- `singularity` — run with [Singularity](https://sylabs.io/docs/)
- `apptainer` — run with [Apptainer](https://apptainer.org/)
- `podman` — run with [Podman](https://podman.io/)
- `conda` — run with [Conda](https://conda.io/docs/) (last resort; prefer containers)
- `test` — complete configuration for automated testing with bundled test data

Multiple profiles can be combined: `-profile test,docker`. Profiles are loaded in order, so later profiles overwrite earlier ones.

### `-resume`

Restart a pipeline from where it left off, reusing cached results for any steps whose inputs have not changed:

```bash
nextflow run ... -resume
```

Supply a specific run name to resume a particular run: `-resume [run-name]`. Use `nextflow log` to list previous run names.

### `-r` (revision)

Pin the pipeline to a specific release for reproducibility:

```bash
nextflow run Ferlab-Ste-Justine/quality-control-pipeline -r 1.0.0 ...
```

Find available releases on the [GitHub releases page](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/releases).

## Resource configuration

Default resource requirements are defined in [`conf/base.config`](../conf/base.config). If a process exits with a retryable error code, Nextflow will automatically resubmit it with increased resources (2×, then 3× original). If it fails after three attempts the pipeline stops.

To override resources for a specific process, use a custom config file passed with `-c`:

```groovy title="custom.config"
process {
    withName: 'PICARD_COLLECTWGSMETRICS' {
        memory = '32.GB'
        time   = '12.h'
    }
}
```

## Running in the background

Use Nextflow's `-bg` flag to detach the pipeline from your terminal session:

```bash
nextflow run ... -bg
```

Alternatively, use `screen`, `tmux`, or submit the Nextflow process itself as a cluster job.

## Nextflow memory requirements

If the Nextflow JVM requests excessive memory, limit it by adding the following to your shell profile (`~/.bashrc` or `~/.bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
