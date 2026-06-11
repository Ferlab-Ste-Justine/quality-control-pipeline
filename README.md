# Ferlab-Ste-Justine/quality-control-pipeline

[![GitHub Actions CI Status](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/actions/workflows/nf-test.yml/badge.svg)](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/actions/workflows/linting.yml/badge.svg)](https://github.com/Ferlab-Ste-Justine/quality-control-pipeline/actions/workflows/linting.yml)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-≥23.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)

## Introduction

**Ferlab-Ste-Justine/quality-control-pipeline** is a Nextflow (DSL2) bioinformatics pipeline for comprehensive quality control of genomic sequencing data. It accepts FASTQ reads, BAM/CRAM alignments, and VCF variant files — running a full suite of QC tools across each data type — and aggregates all results into a single MultiQC HTML report. It can also ingest pre-computed DRAGEN per-sample metrics in place of recomputing them. It is designed to handle multi-sample cohorts with mixed sequencing strategies (WGS, WES, targeted panels) and supports pedigree-based sample identity checks.

### Pipeline steps

1. **BAM/CRAM merging** — merge multi-run alignment files per sample ([Samtools](http://www.htslib.org/))
2. **FASTQ QC** — raw read quality metrics ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)) and sample identity checking ([ngsCheckMate](https://github.com/parklab/NGSCheckMate))
3. **Alignment QC** — alignment statistics (Samtools stats), WGS metrics ([Picard CollectWgsMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360037269351)), sequencing depth ([Mosdepth](https://github.com/brentp/mosdepth)), and DNA contamination estimation ([VerifyBamID2](https://github.com/Griffan/VerifyBamID))
4. **Per-region coverage** — per-gene coverage summaries for up to two custom BED region sets
5. **Sample identity & relatedness** — genetic relatedness checking, per-family or cohort-wide ([Somalier](https://github.com/brentp/somalier))
6. **VCF QC** — variant counts by type and zygosity, Ts/Tv ratio ([BCFtools](https://samtools.github.io/bcftools/))
7. **Report aggregation** — all QC results consolidated into a single interactive report ([MultiQC](http://multiqc.info))

Alignment and variant metrics can alternatively be sourced from pre-computed DRAGEN metrics CSVs via `--dragen_metrics_dir`.

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

Prepare a samplesheet CSV describing your samples (see [usage docs](docs/usage.md) for full column reference):

```csv
participant,sample,familyId,fileType,file1,file2
P001,S001,FAM1,FASTQ,/data/S001_R1.fastq.gz,/data/S001_R2.fastq.gz
P002,S002,FAM1,CRAM,/data/S002.cram,/data/S002.crai
P003,S003,FAM2,VCF,/data/S003.vcf.gz,/data/S003.vcf.gz.tbi
```

Then run the pipeline:

```bash
nextflow run Ferlab-Ste-Justine/quality-control-pipeline \
   -profile docker \
   --input samplesheet.csv \
   --outdir ./results \
   --fasta /path/to/reference.fa \
   --somalier_sites /path/to/sites.vcf.gz
```

For all available parameters, see [docs/usage.md](docs/usage.md). For a description of the output files, see [docs/output.md](docs/output.md).

> [!WARNING]
> Provide pipeline parameters via the CLI or a Nextflow `-params-file`. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) or infrastructural tweaks, not for pipeline parameters.

## Credits

Ferlab-Ste-Justine/quality-control-pipeline was originally written by Georgette Femerling, Lysiane Bouchard.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
