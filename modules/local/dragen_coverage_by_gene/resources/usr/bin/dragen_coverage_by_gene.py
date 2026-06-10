#!/usr/bin/env python3
"""
Aggregate DRAGEN per-interval coverage to per-gene coverage.

Joins DRAGEN's `qc-coverage-region-N_cov_report.bed` (per-interval coverage +
threshold percentages) with `qc-coverage-region-N_read_cov_report.bed`
(per-interval gene assignment) on (chrom, start, end), then groups by gene
and emits the same schema as `aggregate_mosdepth_by_gene.py`:

  #gene  size  totalcvg  count  average_coverage  coverage5  coverage15  ...

Gene names in the read-coverage BED are stored as Python-list-literal strings
(e.g. "['OR4F5']", "['A', 'B']", "[]"). Multi-gene rows contribute to each
listed gene; empty rows are dropped.
"""
import argparse
import ast
import sys
import pandas as pd

__version__ = "1.0.0"


COV_REPORT_COLS = [
    "chrom", "start", "end",
    "total_cvg", "mean_cvg", "Q1_cvg", "median_cvg", "Q3_cvg", "min_cvg", "max_cvg",
    "pct_above_5", "pct_above_15", "pct_above_20", "pct_above_30",
    "pct_above_50", "pct_above_100", "pct_above_200", "pct_above_300",
    "pct_above_400", "pct_above_500", "pct_above_1000",
]

READ_COV_REPORT_COLS = [
    "chrom", "start", "end", "name", "gene_id",
    "total_cvg", "read1_cvg", "read2_cvg",
]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    p.add_argument("-v", "--version", action="version", version=f"%(prog)s {__version__}")
    p.add_argument("--cov-report",      required=True, help="DRAGEN *_cov_report.bed")
    p.add_argument("--read-cov-report", required=True, help="DRAGEN *_read_cov_report.bed")
    p.add_argument("--output",          required=True, help="Output TSV path")
    return p.parse_args()


def read_bed(path, names):
    return pd.read_csv(path, sep="\t", comment="#", header=None, names=names)


def explode_gene(name):
    """Parse the Python-list-literal in DRAGEN's `name` column to a list of genes.

    Drops empty strings and rows with no genes. Returns [] for intergenic so the
    caller can drop those.
    """
    if not isinstance(name, str) or not name.strip():
        return []
    try:
        parsed = ast.literal_eval(name)
    except (ValueError, SyntaxError):
        return [name.strip()]
    if not isinstance(parsed, (list, tuple)):
        parsed = [parsed]
    return [str(g).strip() for g in parsed if str(g).strip()]


def aggregate(cov, read_cov):
    merged = pd.merge(cov, read_cov, on=["chrom", "start", "end"], how="inner",
                      suffixes=("", "_read"))
    merged["size"] = merged["end"] - merged["start"]
    merged["genes"] = merged["name"].apply(explode_gene)
    merged = merged[merged["genes"].map(len) > 0].explode("genes").rename(
        columns={"genes": "gene"})

    threshold_cols = [c for c in cov.columns if c.startswith("pct_above_")]
    # DRAGEN's pct is 0-100; aggregate as length-weighted then convert to 0-1.
    for col in threshold_cols:
        merged[col] = merged[col] * merged["size"] / 100.0  # bases above N reconstructed
    merged["totalcvg"] = merged["total_cvg"].astype(int)

    keep = ["gene", "size", "totalcvg"] + threshold_cols
    grouped = merged[keep].groupby("gene", as_index=False).sum()
    grouped["count"] = merged.groupby("gene", as_index=False).size()["size"]
    grouped["average_coverage"] = grouped["totalcvg"] / grouped["size"]
    for col in threshold_cols:
        # bases-above-N → fraction by dividing by total size for the gene
        new = col.replace("pct_above_", "coverage")
        grouped[new] = grouped[col] / grouped["size"]
    grouped = grouped.drop(columns=threshold_cols)
    grouped = grouped.rename(columns={"gene": "#gene"})

    ordered = ["#gene", "size", "totalcvg", "count", "average_coverage"] + [
        c.replace("pct_above_", "coverage") for c in threshold_cols
    ]
    return grouped[ordered]


def main():
    args = parse_args()
    cov = read_bed(args.cov_report, COV_REPORT_COLS)
    read_cov = read_bed(args.read_cov_report, READ_COV_REPORT_COLS)
    result = aggregate(cov, read_cov)
    if result.empty:
        print("No gene-annotated intervals — nothing to write.", file=sys.stderr)
        sys.exit(0)
    result.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
