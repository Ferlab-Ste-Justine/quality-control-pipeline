# Ferlab-Ste-Justine/quality-control-pipeline: Output

## Introduction

This document describes the output produced by the pipeline. Most plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

All paths below are relative to the top-level output directory (`--outdir`).

## Output directory structure

```
outdir/
├── reports/
│   ├── QC/
│   │   └── {sample.id}/             # Per-sample QC outputs (one directory per sample)
│   │       ├── *_fastqc.html/.zip   # FastQC
│   │       ├── *.vaf                # ngsCheckMate per-sample VAF
│   │       ├── *.txt                # Samtools stats / Mosdepth summaries
│   │       ├── *.wgs_metrics        # Picard CollectWgsMetrics
│   │       ├── *.selfSM             # VerifyBamID2 contamination
│   │       ├── *.somalier           # Somalier extract
│   │       └── *.tsv                # Coverage-by-gene reports
│   │   └── snp_pt/                  # ngsCheckMate batch-level outputs
│   └── pedigree/
│       ├── {family.id}/             # Per-family somalier + PED (cohort_mode=false)
│       │   ├── *.ped
│       │   ├── *.html
│       │   └── *.tsv
│       └── *.html / *.tsv           # Cohort-level somalier outputs (cohort_mode=true)
├── multiqc/
│   ├── multiqc_report.html
│   └── multiqc_data/
├── results/
│   └── {sequencingType}/
│       ├── Merged/                  # Multi-run merged BAM/CRAM files
│       └── ID_rename/               # Reheadered BAM/CRAM files
└── pipeline_info/
```

## Pipeline overview

The pipeline processes FASTQ, BAM/CRAM, and VCF files through the following steps:

- [FastQC](#fastqc) — Raw read quality metrics
- [ngsCheckMate](#ngscheckmate) — FASTQ-level sample identity / cross-contamination
- [Samtools](#samtools) — Alignment statistics and BAM header inspection
- [Picard CollectWgsMetrics](#picard-collectwgsmetrics) — WGS coverage and quality metrics
- [Mosdepth](#mosdepth) — Sequencing depth and coverage
- [Coverage by gene](#coverage-by-gene) — Per-gene coverage summaries for defined regions
- [VerifyBamID2](#verifybamid2) — Alignment-level contamination estimation
- [Somalier](#somalier) — Sample identity and genetic relatedness
- [VCF QC](#vcf-qc) — Variant-level quality metrics
- [DRAGEN metrics input](#dragen-metrics-input) — Pre-computed DRAGEN CSVs used in place of BAM_QC / VCF_QC
- [MultiQC](#multiqc) — Aggregate report summarising all QC results
- [Pipeline information](#pipeline-information) — Nextflow execution reports

---

## FASTQ QC

### FastQC

<details markdown="1">
<summary>Output files</summary>

- `reports/QC/{sample.id}/`
  - `*_fastqc.html`: FastQC report with per-base quality, GC content, adapter content, and overrepresented sequences.
  - `*_fastqc.zip`: Zip archive containing the FastQC report, tab-delimited data, and plot images.

</details>

[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) provides general quality metrics about sequenced reads including quality score distributions, per-base sequence content, and adapter contamination.

### ngsCheckMate

<details markdown="1">
<summary>Output files</summary>

- `reports/QC/{sample.id}/`
  - `*.vaf`: Per-sample VAF (variant allele frequency) file at SNP positions used for identity checking.
- `reports/QC/snp_pt/`
  - `output_all.txt`: All pairwise sample comparisons with correlation scores.
  - `output_matched.txt`: Sample pairs determined to be matched (same individual).
  - `corr_matrix.txt`: Pairwise correlation matrix across all samples.
  - `output.pdf` *(optional)*: Heatmap visualisation of the correlation matrix.

</details>

[ngsCheckMate](https://github.com/parklab/NGSCheckMate) detects sample swaps and cross-contamination by comparing VAF profiles at known SNP positions across all FASTQ inputs.

---

## BAM/CRAM QC

### Samtools

<details markdown="1">
<summary>Output files</summary>

- `reports/QC/{sample.id}/`
  - `*.stats`: Comprehensive alignment statistics (flagstat, base quality, insert size, etc.) from `samtools stats`.
  - `*.txt` *(samples)*: Sample names extracted from the BAM/CRAM read group headers.

</details>

[Samtools](http://www.htslib.org/) is used to merge multi-run BAM/CRAM files, inspect read group headers, and compute alignment statistics.

> **Note:** Merged BAM/CRAM files are published under `results/{sequencingType}/Merged/`. Reheadered files (sample ID repair) are published under `results/{sequencingType}/ID_rename/`.

### Picard CollectWgsMetrics

<details markdown="1">
<summary>Output files</summary>

- `reports/QC/{sample.id}/`
  - `*.wgs_metrics`: WGS coverage metrics including mean coverage, PCT_EXC_* exclusion fractions, and median insert size. Only produced for WGS samples.

</details>

[Picard CollectWgsMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360037269351-CollectWgsMetrics-Picard) summarises whole-genome sequencing coverage and base quality metrics.

### Mosdepth

<details markdown="1">
<summary>Output files</summary>

- `reports/QC/{sample.id}/`
  - `*.mosdepth.global.dist.txt`: Cumulative coverage distribution across the genome.
  - `*.mosdepth.region.dist.txt` *(if regions BED provided)*: Coverage distribution per region.
  - `*.mosdepth.summary.txt`: Mean coverage per chromosome and total.
  - `*.regions.bed.gz` *(QC coverage regions)*: Per-interval mean coverage for defined QC regions.
  - `*.thresholds.bed.gz` *(QC coverage regions)*: Fraction of bases at coverage thresholds (5x, 15x, 20x, 30x, 50x, 100x, 200x, 300x, 400x, 500x, 1000x).

</details>

[Mosdepth](https://github.com/brentp/mosdepth) computes fast per-base and per-region sequencing depth. It is run twice: once for overall alignment QC and once per QC coverage region set (up to two region BED files).

### Coverage by gene

<details markdown="1">
<summary>Output files</summary>

- `reports/QC/{sample.id}/`
  - `*.<region_name>.tsv`: Per-gene mean coverage and coverage-at-threshold summary, derived from Mosdepth output for each QC coverage region.

</details>

A local module aggregates Mosdepth region and threshold outputs into a per-gene coverage table. Output file names reflect the `--region_1_name` / `--region_2_name` parameters (defaults: `qc_regions_1`, `qc_regions_2`).

### VerifyBamID2

<details markdown="1">
<summary>Output files</summary>

- `reports/QC/{sample.id}/`
  - `*.selfSM`: Per-sample contamination estimate. Key field: `FREEMIX` (estimated contamination fraction).

</details>

[VerifyBamID2](https://github.com/Griffan/VerifyBamID) estimates DNA contamination from aligned reads by comparing allele frequencies at known SNP positions against population reference SVD files.

---

## Somalier

<details markdown="1">
<summary>Output files — per-family mode (`--cohort_mode false`)</summary>

- `reports/QC/{sample.id}/`
  - `*.somalier`: Somalier extract file per sample (fingerprint at sites VCF positions).
- `reports/pedigree/{family.id}/`
  - `*.ped`: Per-family PED file generated from the input pedigree.
  - `*.html`: Interactive HTML relatedness report for the family.
  - `*.pairs.tsv`: Pairwise relatedness metrics (IBS0, IBS2, kinship coefficient).
  - `*.samples.tsv`: Per-sample summary (het rate, depth, ancestry PCs).

</details>

<details markdown="1">
<summary>Output files — cohort mode (default)</summary>

- `reports/QC/{sample.id}/`
  - `*.somalier`: Somalier extract file per sample.
- `reports/pedigree/`
  - `*.html`: Interactive HTML relatedness report for the entire cohort.
  - `*.pairs.tsv`: Pairwise relatedness metrics across all samples.
  - `*.samples.tsv`: Per-sample summary across the cohort.

</details>

[Somalier](https://github.com/brentp/somalier) checks sample identity and genetic relatedness by extracting genotype-like information at known sites. It can run per-family (using a pedigree file split by family) or across the entire cohort.

---

## VCF QC

<details markdown="1">
<summary>Output files</summary>

- `reports/QC/{sample.id}/`
  - `*.metrics.json`: Variant counts including SNVs, insertions, deletions, het/hom breakdowns, and Ts/Tv ratio.
  - `*.stats`: Full bcftools stats output.

</details>

VCF QC uses [BCFtools](https://samtools.github.io/bcftools/) to compute per-sample variant statistics. Counts are broken down by variant class (SNV, insertion, deletion) and zygosity (het/hom), and the Ts/Tv ratio is extracted from `bcftools stats`.

---

## DRAGEN metrics input

When `--dragen_metrics_dir` is set, BAM_QC and VCF_QC are skipped and the report's alignment / variant metrics are populated from DRAGEN's pre-computed CSVs. The pipeline does not republish these files — they're staged into the MultiQC work directory and surfaced through the MultiQC report and the per-sample JSON sidecars below.

Files are picked up at any depth under `--dragen_metrics_dir` (local path or `s3://` / `gs://` / `az://` URI). Both `<sample>.<type>.csv` and `<sample>.final.<type>.csv` filenames are recognised.

| Filename suffix | Used for |
|---|---|
| `*.mapping_metrics.csv` | Alignment metrics + Q30 yield + estimated contamination |
| `*.wgs_coverage_metrics.csv` | Mean autosome coverage, uniformity (MAD proxy), %15x |
| `*.vc_metrics.csv` | SNV / insertion / deletion counts, Ti/Tv, het:hom ratios |
| `*.ploidy_estimation_metrics.csv` | XX / XY ploidy → predicted-sex fallback for `sex_check` when somalier did not run for that sample |

Somalier still runs against any BAM/CRAM provided in the samplesheet for pedigree validation. Samples without BAM/CRAM (e.g. GVCF-only with DRAGEN metrics) still appear in the per-family report with their `pedigree_sex` from the samplesheet-derived PED, and pick up `sex_check` from the DRAGEN ploidy estimate.

---

## MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: Standalone HTML report summarising QC results from all tools and samples.
  - `multiqc_data/`: Directory of parsed, tab-delimited statistics from each tool.

</details>

[MultiQC](http://multiqc.info) aggregates QC outputs from all tools in the pipeline into a single interactive report. The report is configured via `assets/multiqc_config.yml`.

---

## Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - `execution_report.html`: Nextflow execution report with task-level resource usage.
  - `execution_timeline.html`: Timeline of task execution across the run.
  - `execution_trace.txt`: Tab-delimited trace of every task (CPU, memory, wall time).
  - `pipeline_dag.dot` / `pipeline_dag.svg`: Directed acyclic graph of the workflow.
  - `params.json`: Parameters used for the pipeline run.
  - `quality-control-pipeline_software_mqc_versions.yml`: Software versions for all tools.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) generates execution reports, timelines, and traces that are useful for troubleshooting and auditing pipeline runs.
