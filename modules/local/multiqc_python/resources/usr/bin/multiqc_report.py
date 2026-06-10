#!/usr/bin/env python3
"""MultiQC report generation for the quality-control-pipeline.

Produces a GA4GH QC compliant report

Inputs (all CLI args):
    files                  positional, list of QC files staged by Nextflow
    --config <yaml>        MultiQC config (passed through to parse_logs)
    --title <str>          report title (passed through to write_report)
    --filename <str>       output report filename
    --mode {cohort,family} affects which sections are shown:
                             - family: pedigree table; somalier built-in plots
                               are removed (see assets/multiqc_config.yml).
                             - cohort: skip pedigree table; keep somalier plots.
    --ped <ped>            PED file driving the canonical sample list — every
                           sample in the PED appears in every table even if it
                           has no QC data for that section (cells render as na).
    --thresholds <yaml>    QC thresholds; defaults baked in if omitted.
    --json-out <dir>       directory to write per-sample JSON files.
"""

import argparse
import json
import logging
import os
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import multiqc
import yaml
from multiqc.plots import table

log = logging.getLogger("multiqc")


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

@dataclass
class SampleMetrics:
    """Per-sample QC metrics.
    """
    # aln_metrics
    total_reads: Optional[float] = None
    mapped_reads: Optional[float] = None
    duplicate_reads: Optional[float] = None
    insert_size_std_deviation: Optional[float] = None
    mad_autosome_coverage: Optional[float] = None  # GA4GH coverage_uniformity
    mean_autosome_coverage: Optional[float] = None
    mean_insert_size: Optional[float] = None
    pct_autosomes_15x: Optional[float] = None
    pct_reads_mapped: Optional[float] = None
    pct_reads_properly_paired: Optional[float] = None
    yield_bp_q30: Optional[int] = None
    cross_contamination_rate: Optional[float] = None
    # variant_metrics
    count_deletions: Optional[int] = None
    count_insertions: Optional[int] = None
    count_snvs: Optional[int] = None
    ratio_heterozygous_homozygous_indel: Optional[float] = None
    ratio_heterozygous_homozygous_snv: Optional[float] = None
    ratio_insertion_deletion: Optional[float] = None
    ratio_transitions_transversions_snv: Optional[float] = None
    count_het_snvs: Optional[int] = None
    count_hom_snvs: Optional[int] = None
    count_het_indels: Optional[int] = None
    count_hom_indels: Optional[int] = None
    # pedigree.
    family_id: Optional[str] = None
    pedigree_sex: Optional[str] = None  # 1=male, 2=female (PED convention)
    somalier_sex: Optional[str] = None
    sex_check: str = "na"               # pass | fail | na; set by parse_somalier




@dataclass
class PedRow:
    family_id: str
    sample_id: str
    paternal_id: str
    maternal_id: str
    sex: str
    phenotype: str


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _as_float(v):
    if v is None or v == "" or v == "?":
        return None
    try:
        return float(v)
    except (ValueError, TypeError):
        return None


def _as_int(v):
    f = _as_float(v)
    return int(f) if f is not None else None


_PASSFAIL_OPS = {
    ">=": lambda v, t: v >= t,
    "<=": lambda v, t: v <= t,
    ">":  lambda v, t: v > t,
    "<":  lambda v, t: v < t,
}


def _passfail(*checks):
    """Combine one or more (value, threshold, op) checks into 'pass'|'fail'|'na'.

    Returns 'na' if any value/threshold is missing, 'fail' if any check fails,
    otherwise 'pass'. Use a single tuple for a one-shot threshold check, or
    multiple tuples for combined criteria (all must pass).
    """
    for value, threshold, op in checks:
        if value is None or threshold is None:
            return "na"
        if not _PASSFAIL_OPS[op](value, threshold):
            return "fail"
    return "pass"


def _sample_name_from_filename(path):
    """Derive a sample name by stripping everything from the first dot.

    Used only for gene-coverage CSVs (no sample column inside the file).
    """
    bn = os.path.basename(path)
    if "." in bn:
        bn = bn.split(".", 1)[0]
    return bn.strip("._-")


# ---------------------------------------------------------------------------
# Configuration loading
# ---------------------------------------------------------------------------

def parse_thresholds(path):
    """Load QC thresholds from a YAML file. Nextflow always stages
    assets/qc_thresholds.yml when params.qc_thresholds is unset, so this path
    is expected to exist in normal pipeline use.
    """
    if not path:
        raise ValueError(
            "No QC thresholds file provided. Pass --thresholds <yaml> or set params.qc_thresholds."
        )
    with open(path) as fh:
        loaded = yaml.safe_load(fh) or {}
    if not isinstance(loaded, dict):
        raise ValueError(f"Thresholds file {path} must be a mapping, got {type(loaded).__name__}")
    return loaded


def parse_ped(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                log.warning("PED line skipped (need 6 columns) got %i",len(parts))
                continue
            rows.append(PedRow(*parts[:6]))
    return rows


# ---------------------------------------------------------------------------
# Per-module parsers
# ---------------------------------------------------------------------------

def _get_module(module):
    try:
        return multiqc.get_module_data(module=module) or {}
    except ValueError:
        return {}


# Map samtools stats fields -> SampleMetrics attributes (all coerced to float).
_SAMTOOLS_FIELDS = {
    "total_reads":                "raw_total_sequences",
    "mapped_reads":               "reads_mapped",
    "duplicate_reads":            "reads_duplicated",
    "pct_reads_mapped":           "reads_mapped_percent",
    "pct_reads_properly_paired":  "reads_properly_paired_percent",
    "mean_insert_size":           "insert_size_average",
    "insert_size_std_deviation":  "insert_size_standard_deviation",
}


def parse_samtools(samples):
    """Read per-sample samtools stats. MultiQC shape: {sample: {field: value}}."""
    for sample, row in _get_module("samtools").items():
        if not isinstance(row, dict) or "raw_total_sequences" not in row:
            continue
        s = samples.setdefault(sample, SampleMetrics())
        for attr, src in _SAMTOOLS_FIELDS.items():
            setattr(s, attr, _as_float(row.get(src)))


def parse_picard(samples):
    """Read picard sub-tools: CollectWgsMetrics, CollectHsMetrics, CollectQualityYieldMetrics.

    MultiQC shape: {sub_tool: {sample: {field: value}}}.
    """
    data = _get_module("picard")
    if not isinstance(data, dict):
        return
    for sub_data in data.values():
        if not isinstance(sub_data, dict):
            continue
        for sample, row in sub_data.items():
            if not isinstance(row, dict):
                continue
            s = samples.setdefault(sample, SampleMetrics())
            # WgsMetrics
            if "MEAN_COVERAGE" in row:
                s.mean_autosome_coverage = _as_float(row.get("MEAN_COVERAGE"))
                s.mad_autosome_coverage = _as_float(row.get("MAD_COVERAGE"))
                s.pct_autosomes_15x = _as_float(row.get("PCT_15X"))
            # HsMetrics (targeted/WXS) — TARGET-based equivalents. HsMetrics has no
            # 15X bucket (jumps 10X -> 20X); use 20X as the closest proxy. HsMetrics
            # has no MAD coverage, so mad_autosome_coverage stays None for WXS.
            elif "MEAN_TARGET_COVERAGE" in row:
                s.mean_autosome_coverage = _as_float(row.get("MEAN_TARGET_COVERAGE"))
                s.pct_autosomes_15x = _as_float(row.get("PCT_TARGET_BASES_20X"))
            # CollectQualityYieldMetrics (GA4GH yield_bp_q30)
            if "PF_Q30_BASES" in row:
                s.yield_bp_q30 = _as_int(row.get("PF_Q30_BASES"))


def parse_verifybamid(samples):
    """Read verifybamid FREEMIX as cross_contamination."""
    for sample, row in _get_module("verifybamid").items():
        if not isinstance(row, dict) or "FREEMIX" not in row:
            continue
        s = samples.setdefault(sample, SampleMetrics())
        s.cross_contamination_rate = _as_float(row.get("FREEMIX"))


# Normalize somalier sex values (somalier writes `sex` as 1.0/2.0 and
# `original_pedigree_sex` as "male"/"female").
_SEX_NORM = {"1": "male", "1.0": "male", "male": "male",
             "2": "female", "2.0": "female", "female": "female"}


def parse_somalier(samples, thresholds):
    """Read somalier per-sample sex predictions and pair-level relatedness.

    MultiQC's somalier module keys entries by:
      - sample name for per-sample rows (containing `sex` + `original_pedigree_sex`)
      - "<sample_a>*<sample_b>" for pair rows (containing `relatedness`)

    For per-sample rows we set somalier_sex and precompute sex_check by comparing
    against `original_pedigree_sex` (which somalier copies from the input PED).
    Pair rows are pre-enriched with `expected_relationship` and `verdict` so the
    downstream report + per-sample aggregation just read those fields.
    """
    pairs = []
    for key, val in _get_module("somalier").items():
        if not isinstance(val, dict):
            continue
        if "sex" in val:
            s = samples.setdefault(key, SampleMetrics())
            s.somalier_sex = val["sex"]
            ped = _SEX_NORM.get(str(val.get("original_pedigree_sex", "")).lower())
            som = _SEX_NORM.get(str(val["sex"]).lower())
            s.sex_check = "na" if (ped is None or som is None) else ("pass" if ped == som else "fail")
        if "relatedness" in val and isinstance(key, str) and "*" in key:
            sa, sb = key.split("*", 1)
            rel = _as_float(val.get("relatedness"))
            exp = _as_float(val.get("expected_relatedness"))
            pairs.append({
                "sample_a": sa,
                "sample_b": sb,
                "relatedness": rel,
                "expected_relatedness": exp,
                "ibs0": val.get("ibs0"),
                "ibs2": val.get("ibs2"),
                "expected_relationship": (expected_relationship(exp) or EXPECTED_RELATIONSHIP_OTHER)["name"],
                "verdict": classify_pair_relatedness(rel, exp, thresholds),
            })
    return pairs


# Map VCF_QC's per-sample JSON fields -> SampleMetrics attribute (+ coercer).
_VCF_INT_FIELDS = {
    "count_snvs":        "total_snvs",
    "count_insertions":  "total_insertions",
    "count_deletions":   "total_deletions",
    "count_het_snvs":    "heterozygous_snvs",
    "count_hom_snvs":    "homozygous_snvs",
    "count_het_indels":  "heterozygous_indels",
    "count_hom_indels":  "homozygous_indels",
}
_VCF_FLOAT_FIELDS = {
    "ratio_transitions_transversions_snv": "ti_tv_ratio",
    "ratio_insertion_deletion":            "ins_dels_ratio",
    "ratio_heterozygous_homozygous_snv":   "het_hom_ratio_snvs",
    "ratio_heterozygous_homozygous_indel": "het_hom_ratio_indels",
}


def parse_vcf_metrics(files, samples):
    """Parse per-sample VCF metrics JSON files produced by VCF_QC."""
    for path in files:
        try:
            with open(path) as fh:
                j = json.load(fh)
        except Exception as e:
            log.warning(f"Could not parse VCF metrics file {path}: {e}")
            continue
        if not isinstance(j, dict) or "sample_id" not in j:
            continue
        s = samples.setdefault(j["sample_id"], SampleMetrics())
        for attr, src in _VCF_INT_FIELDS.items():
            setattr(s, attr, _as_int(j.get(src)))
        for attr, src in _VCF_FLOAT_FIELDS.items():
            setattr(s, attr, _as_float(j.get(src)))


def _read_dragen_csv(path):
    """Return {section: {metric: (value_str, pct_str_or_None)}} for a DRAGEN CSV.

    DRAGEN metric CSVs have rows: `SECTION,subgroup,metric,value[,pct]`. We
    ignore the subgroup column and key by (section, metric). When a metric is
    repeated under different subgroups (e.g. per read-group) the SUMMARY row
    (subgroup empty) wins -- which matches the GA4GH per-sample view.
    """
    out = defaultdict(dict)
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split(",")
            if len(parts) < 4:
                continue
            section, subgroup, metric = parts[0], parts[1], parts[2]
            value = parts[3]
            pct = parts[4] if len(parts) >= 5 else None
            if subgroup == "" or section not in out or metric not in out[section]:
                out[section][metric] = (value, pct)
    return out


def parse_dragen_csv_files(files, samples):
    """Populate SampleMetrics from DRAGEN per-sample CSVs.

    Recognised filenames (under any directory):
      *.mapping_metrics.csv           -> alignment metrics + Q30 + contamination
      *.wgs_coverage_metrics.csv      -> coverage + uniformity (proxy for MAD)
      *.vc_metrics.csv                -> variant calling counts and ratios
      *.ploidy_estimation_metrics.csv -> XX/XY -> somalier_sex fallback

    Sample names are taken from the filename prefix (everything before the
    first `.`), which matches DRAGEN's `<sample>.final.<metric>.csv` convention.
    """
    for path in files:
        name = os.path.basename(path)
        sample = name.split(".", 1)[0]
        s = samples.setdefault(sample, SampleMetrics())
        data = _read_dragen_csv(path)
        if name.endswith(".mapping_metrics.csv"):
            row = data.get("MAPPING/ALIGNING SUMMARY", {})
            s.total_reads               = _as_float(row.get("Total input reads", (None,))[0])
            s.mapped_reads              = _as_float(row.get("Mapped reads", (None,))[0])
            s.duplicate_reads           = _as_float(row.get("Number of duplicate marked reads", (None,))[0])
            s.pct_reads_mapped          = _as_float((row.get("Mapped reads", (None, None)) + (None,))[1])
            s.pct_reads_properly_paired = _as_float((row.get("Properly paired reads", (None, None)) + (None,))[1])
            s.mean_insert_size          = _as_float(row.get("Insert length: mean", (None,))[0])
            s.insert_size_std_deviation = _as_float(row.get("Insert length: standard deviation", (None,))[0])
            s.yield_bp_q30              = _as_int(row.get("Q30 bases", (None,))[0])
            s.cross_contamination_rate  = _as_float(row.get("Estimated sample contamination", (None,))[0])
        elif name.endswith(".wgs_coverage_metrics.csv"):
            row = data.get("COVERAGE SUMMARY", {})
            s.mean_autosome_coverage = _as_float(row.get("Average autosomal coverage over genome", (None,))[0])
            s.pct_autosomes_15x      = _as_float((row.get("PCT of genome with coverage [  15x: inf)", (None, None)) + (None,))[1]
                                                 or row.get("PCT of genome with coverage [  15x: inf)", (None,))[0])
            # DRAGEN has no MAD coverage; use the "Uniformity of coverage" proxy.
            s.mad_autosome_coverage  = _as_float(row.get("Uniformity of coverage (PCT > 0.2*mean) over genome", (None,))[0])
        elif name.endswith(".vc_metrics.csv"):
            row = data.get("VARIANT CALLER POSTFILTER", {})
            s.count_snvs       = _as_int(row.get("SNPs", (None,))[0])
            ins_hom = _as_int(row.get("Insertions (Hom)", (None,))[0]) or 0
            ins_het = _as_int(row.get("Insertions (Het)", (None,))[0]) or 0
            del_hom = _as_int(row.get("Deletions (Hom)", (None,))[0]) or 0
            del_het = _as_int(row.get("Deletions (Het)", (None,))[0]) or 0
            s.count_insertions = ins_hom + ins_het
            s.count_deletions  = del_hom + del_het
            het = _as_float(row.get("Heterozygous", (None,))[0])
            hom = _as_float(row.get("Homozygous", (None,))[0])
            s.ratio_heterozygous_homozygous_snv = (het / hom) if (het is not None and hom) else None
            s.ratio_transitions_transversions_snv = _as_float(row.get("Ti/Tv ratio", (None,))[0])
            s.ratio_insertion_deletion = (s.count_insertions / s.count_deletions) if s.count_deletions else None
        elif name.endswith(".ploidy_estimation_metrics.csv"):
            row = data.get("PLOIDY ESTIMATION", {})
            ploidy = row.get("Ploidy estimation", (None,))[0]
            if ploidy in ("XX", "XY") and not s.somalier_sex:
                # Use DRAGEN's ploidy as the somalier_sex fallback so sex_check
                # has a value when somalier didn't run (no BAM provided).
                s.somalier_sex = "female" if ploidy == "XX" else "male"
                # Refresh sex_check using PED's pedigree_sex (set by parse_ped).
                ped = _SEX_NORM.get(str(s.pedigree_sex).lower()) if s.pedigree_sex else None
                som = _SEX_NORM.get(s.somalier_sex)
                s.sex_check = "na" if (ped is None or som is None) else ("pass" if ped == som else "fail")


def parse_gene_coverage_files(files):
    """Parse the coverage_by_gene TSVs (from COVERAGE_BY_GENE module).

    Header columns: gene, size, totalcvg, count, average_coverage,
    coverage5, coverage15, coverage20, coverage30, coverage50, coverage100, ...
    """
    rows = []
    wanted = {"gene", "average_coverage", "coverage15", "coverage30", "coverage100"}
    for path in files:
        sample = _sample_name_from_filename(path)
        try:
            with open(path) as fh:
                header = None
                for line in fh:
                    if not line.strip():
                        continue
                    if line.startswith("#"):
                        # First # line is the header, e.g. "#gene\tsize\t..."
                        if header is None:
                            header = line.lstrip("#").rstrip("\n").split("\t")
                        continue
                    if header is None:
                        log.warning(f"{path}: data line before header, skipping: {line.rstrip()!r}")
                        continue
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) != len(header):
                        log.warning(f"{path}: column count mismatch, skipping line")
                        continue
                    record = dict(zip(header, fields))
                    if not wanted.issubset(record):
                        log.warning(f"{path}: expected columns missing ({wanted - set(record)}), skipping file")
                        break
                    try:
                        rows.append({
                            "sample": sample,
                            "gene": record["gene"],
                            "average_coverage": float(record["average_coverage"]),
                            "coverage15": float(record["coverage15"]),
                            "coverage30": float(record["coverage30"]),
                            "coverage100": float(record["coverage100"]),
                        })
                    except ValueError as e:
                        log.warning(f"Could not parse line in {path}: {line.rstrip()!r} - {e}")
        except Exception as e:
            log.warning(f"Could not read gene-coverage file {path}: {e}")
    return rows


# ---------------------------------------------------------------------------
# Pedigree analysis
# ---------------------------------------------------------------------------

# Matches the lookup used in MultiQC's somalier module
# (multiqc/modules/somalier/somalier.py — `relatedness_groups`).
# Somalier also writes expected_relatedness = -1 for unrelated samples from
# different families; that module normalises -1 -> 0 before this lookup.
EXPECTED_RELATIONSHIP_GROUPS = {
    0:    {"name": "Unrelated",    "color": "#cccccc"},
    0.49: {"name": "Sib-sib",      "color": "#f37b28"},
    0.5:  {"name": "Parent-child", "color": "#2f9f46"},
}
EXPECTED_RELATIONSHIP_OTHER = {"name": "Other", "color": "#4a7cb6"}
_ALL_RELATIONSHIP_GROUPS = list(EXPECTED_RELATIONSHIP_GROUPS.values()) + [EXPECTED_RELATIONSHIP_OTHER]


def expected_relationship(exp):
    """Bucket somalier `expected_relatedness` into a {'name', 'color'} dict.

    Somalier sets expected_relatedness = -1 for cross-family pairs and 0 for
    within-family unrelated pairs; both are normalised to 0.
    Returns None for missing values (callers should treat as 'na').
    """
    e = _as_float(exp)
    if e is None:
        return None
    if e == -1:
        e = 0
    return EXPECTED_RELATIONSHIP_GROUPS.get(e, EXPECTED_RELATIONSHIP_OTHER)


def classify_pair_relatedness(rel, exp, thresholds):
    """Bucket a somalier pair's relatedness into pass | warn | fail | na.

    Bucketing depends on the expected relationship implied by `exp`:

      Unrelated (exp = -1 or 0):
        rel <  pedigree_unrelated_pass_max               -> pass
        pedigree_unrelated_pass_max <= rel <= warn_max   -> warn (possible consanguinity)
        rel >  pedigree_unrelated_warn_max               -> fail

      Sib-sib (exp = 0.49) / Parent-child (exp = 0.50):
        pedigree_related_pass_min <= rel <= related_pass_max -> pass
        otherwise                                            -> fail

      Other expected values: pass if |rel - exp| <= 0.10, else fail.
    """
    if rel is None or exp is None:
        return "na"
    group = expected_relationship(exp)
    label = group["name"] if group else "na"
    if label == "Unrelated":
        # Negative observed relatedness is "more unrelated" than 0 — passes.
        if rel < thresholds["pedigree_unrelated_pass_max"]:
            return "pass"
        if rel <= thresholds["pedigree_unrelated_warn_max"]:
            return "warn"
        return "fail"
    if label in ("Sib-sib", "Parent-child"):
        if thresholds["pedigree_related_pass_min"] <= rel <= thresholds["pedigree_related_pass_max"]:
            return "pass"
        return "fail"
    # "Other" or "na": ±0.10 around expected when both numbers are present.
    return "pass" if abs(rel - exp) <= 0.10 else "fail"


def pedigree_validation_by_sample(pairs):
    """Aggregate per-pair verdicts into per-sample verdicts.

    Worst verdict wins (fail > warn > pass). Samples with no comparisons are absent
    from the result; callers default missing samples to 'na'.
    """
    SEVERITY = {"pass": 0, "warn": 1, "fail": 2}
    worst = {}
    for p in pairs:
        v = p["verdict"]
        if v == "na":
            continue
        for s in (p["sample_a"], p["sample_b"]):
            if SEVERITY[v] >= SEVERITY.get(worst.get(s, "pass"), 0):
                worst[s] = v
    return worst


# ---------------------------------------------------------------------------
# Section builders
# ---------------------------------------------------------------------------

# Reused for any column that renders pass | warn | fail | na badges with the
# usual green / yellow / red / grey palette.
PASS_WARN_FAIL_FORMATTING = {
    "cond_formatting_rules": {
        "pass": [{"s_eq": "pass"}],
        "warn": [{"s_eq": "warn"}],
        "fail": [{"s_eq": "fail"}],
        "na":   [{"s_eq": "na"}],
    },
    "cond_formatting_colours": [
        {"pass": "#5cb85c"},
        {"warn": "#f0ad4e"},
        {"fail": "#d9534f"},
        {"na":   "#cccccc"},
    ],
}


GA4GH_ALIGNMENT_HEADERS = {
    "yield_bp_q30":              {"title": "Yield ≥Q30 (bp)",          "description": "GA4GH yield_bp_q30; Picard CollectQualityYieldMetrics PF_Q30_BASES.",                                                        "format": "{:,.0f}"},
    "pct_reads_mapped":          {"title": "% Mapped",                  "description": "GA4GH pct_reads_mapped; samtools stats.",                                                                                  "format": "{:,.2f}%"},
    "pct_reads_properly_paired": {"title": "% Properly Paired",         "description": "GA4GH pct_reads_properly_paired; samtools stats.",                                                                         "format": "{:,.2f}%"},
    "mean_insert_size":          {"title": "Mean Insert Size",          "description": "GA4GH mean_insert_size; samtools stats.",                                                                                  "format": "{:,.2f}"},
    "insert_size_std_deviation": {"title": "Insert Size SD",            "description": "GA4GH insert_size_std_deviation; samtools stats.",                                                                         "format": "{:,.2f}"},
    "mean_autosome_coverage":    {"title": "Mean Autosome Coverage",    "description": "GA4GH mean_autosome_coverage; Picard CollectWgsMetrics MEAN_COVERAGE (WGS) or CollectHsMetrics MEAN_TARGET_COVERAGE (WXS).",  "format": "{:,.2f}X"},
    "pct_autosomes_15x":         {"title": "% Autosomes ≥15X",          "description": "GA4GH pct_autosomes_15x; Picard PCT_15X (WGS) or PCT_TARGET_BASES_20X used as proxy for WXS (HsMetrics has no 15X bucket).", "format": "{:,.2f}%"},
    "mad_autosome_coverage":     {"title": "Coverage Uniformity (MAD)", "description": "GA4GH coverage_uniformity (canonical ID mad_autosome_coverage); Picard MAD_COVERAGE. Not produced by HsMetrics.",            "format": "{:,.2f}"},
    "cross_contamination_rate":  {"title": "Cross-contamination",       "description": "GA4GH cross_contamination_rate; VerifyBamID2 FREEMIX.",                                                                     "format": "{:,.4f}"},
    "total_reads":               {"title": "Total Reads",               "description": "Total reads (raw_total_sequences).",  "format": "{:,.0f}", "hidden": True},
    "mapped_reads":              {"title": "Mapped Reads",              "description": "Mapped reads.",                       "format": "{:,.0f}", "hidden": True},
    "duplicate_reads":           {"title": "Duplicate Reads",           "description": "Duplicate reads.",                    "format": "{:,.0f}", "hidden": True},
}

GA4GH_VARIANT_HEADERS = {
    "count_snvs":                          {"title": "SNVs",         "description": "GA4GH count_snvs.",                                  "format": "{:,.0f}"},
    "count_het_snvs":                      {"title": "Het SNVs",     "description": "Heterozygous SNVs.",                                 "format": "{:,.0f}"},
    "count_hom_snvs":                      {"title": "Hom SNVs",     "description": "Homozygous SNVs.",                                   "format": "{:,.0f}"},
    "count_insertions":                    {"title": "Insertions",   "description": "GA4GH count_insertions.",                            "format": "{:,.0f}"},
    "count_deletions":                     {"title": "Deletions",    "description": "GA4GH count_deletions.",                             "format": "{:,.0f}"},
    "count_het_indels":                    {"title": "Het Indels",   "description": "Heterozygous indels.",                               "format": "{:,.0f}"},
    "count_hom_indels":                    {"title": "Hom Indels",   "description": "Homozygous indels.",                                 "format": "{:,.0f}"},
    "ratio_transitions_transversions_snv": {"title": "Ti/Tv",        "description": "GA4GH ratio_transitions_transversions_snv.",        "format": "{:,.2f}"},
    "ratio_insertion_deletion":            {"title": "Ins/Del",      "description": "GA4GH ratio_insertion_deletion.",                    "format": "{:,.2f}"},
    "ratio_heterozygous_homozygous_snv":   {"title": "Het/Hom SNV",  "description": "GA4GH ratio_heterozygous_homozygous_snv.",           "format": "{:,.2f}"},
    "ratio_heterozygous_homozygous_indel": {"title": "Het/Hom Indel","description": "GA4GH ratio_heterozygous_homozygous_indel.",         "format": "{:,.2f}"},
}


def add_general_status_section(module, samples, thresholds, pedigree_pass):
    """One pass/fail row per sample across the major QC dimensions."""
    rows = {}
    t = thresholds
    for sample in sorted(samples):
        s = samples[sample]
        rows[sample] = {
            "aln_quality":         _passfail((s.pct_reads_mapped,          t["pct_reads_mapped"],          ">="),
                                             (s.pct_reads_properly_paired, t["pct_reads_properly_paired"], ">=")),
            "coverage":            _passfail((s.mean_autosome_coverage,    t["mean_autosome_coverage"],    ">=")),
            "contamination":       _passfail((s.cross_contamination_rate,  t["cross_contamination_rate"], "<=")),
            "sex_check":           s.sex_check,
            "pedigree_validation": pedigree_pass.get(sample, "na"),
            "mean_coverage":       s.mean_autosome_coverage,
        }
    headers = {
        "aln_quality":         {"title": "Aln Quality",         "description": f"% mapped ≥ {thresholds['pct_reads_mapped']} and % properly paired ≥ {thresholds['pct_reads_properly_paired']}", "scale": "pass-fail"},
        "coverage":            {"title": "Coverage",            "description": f"Mean autosome coverage ≥ {thresholds['mean_autosome_coverage']}X",                                             "scale": "pass-fail"},
        "contamination":       {"title": "Contamination",       "description": f"FREEMIX ≤ {thresholds['cross_contamination_rate']}",                                                          "scale": "pass-fail"},
        "sex_check":           {"title": "Sex Check",           "description": "Inferred sex matches PED-recorded sex",                                                                "scale": "pass-fail"},
        "pedigree_validation": {
            "title": "Pedigree Validation",
            "description": "Aggregated across all somalier pair comparisons for the sample (any fail -> fail; any warn -> warn; otherwise pass).",
            **PASS_WARN_FAIL_FORMATTING,
        },
        "mean_coverage":       {"title": "Mean Coverage",       "description": "Mean autosome (WGS) / target (WXS) coverage",                                                                    "format": "{:,.2f}X"},
    }
    module.add_section(
        name="General Status",
        anchor="general-status",
        description="Per-sample pass/fail summary across QC dimensions.",
        plot=table.plot(rows, headers=headers),
    )


def add_metric_table(module, samples, *, name, anchor, headers):
    """Render a per-sample metric table; skip the section if every cell is None."""
    rows = {sample: {k: getattr(s, k) for k in headers} for sample, s in samples.items()}
    if not any(v is not None for row in rows.values() for v in row.values()):
        return
    module.add_section(
        name=name,
        anchor=anchor,
        description=f"Per-sample {name.lower()}.",
        plot=table.plot(rows, headers=headers),
    )


def add_pedigree_section(module, pairs, thresholds):
    """Pair-level pedigree validation table. Reads precomputed `verdict` and
    `expected_relationship` from each pair (set by parse_somalier).
    """
    if not pairs:
        return
    rows = {f"{p['sample_a']} | {p['sample_b']} [{i}]": p for i, p in enumerate(pairs)}
    unrelated_warn = f"{thresholds['pedigree_unrelated_pass_max']*100:.0f}-{thresholds['pedigree_unrelated_warn_max']*100:.0f}%"
    related_pass = f"{thresholds['pedigree_related_pass_min']*100:.0f}-{thresholds['pedigree_related_pass_max']*100:.0f}%"
    headers = {
        "sample_a":              {"title": "Sample A"},
        "sample_b":              {"title": "Sample B"},
        "expected_relatedness":  {"title": "Expected Relatedness", "format": "{:,.3f}"},
        "expected_relationship": {
            "title": "Expected Relationship",
            "description": "Pedigree-implied relationship derived from somalier's `expected_relatedness` (Unrelated / Sib-sib / Parent-child / Other).",
            # Colour each label with the same palette as the somalier relatedness plot.
            "cond_formatting_rules":   {g["name"]: [{"s_eq": g["name"]}] for g in _ALL_RELATIONSHIP_GROUPS},
            "cond_formatting_colours": [{g["name"]: g["color"]} for g in _ALL_RELATIONSHIP_GROUPS],
        },
        "relatedness":           {"title": "Relatedness",          "format": "{:,.3f}"},
        "ibs0":                  {"title": "IBS0",                 "format": "{:,.0f}"},
        "ibs2":                  {"title": "IBS2",                 "format": "{:,.0f}"},
        "verdict": {
            "title": "Validation",
            "description": (
                f"Unrelated: pass <{thresholds['pedigree_unrelated_pass_max']*100:.0f}%, "
                f"warn {unrelated_warn} (possible consanguinity), "
                f"fail >{thresholds['pedigree_unrelated_warn_max']*100:.0f}%. "
                f"Sib-sib / Parent-child: pass {related_pass}, else fail."
            ),
            **PASS_WARN_FAIL_FORMATTING,
        },
    }
    module.add_section(
        name="Pedigree",
        anchor="pedigree",
        description="Pair-level pedigree validation from Somalier (relatedness vs expected_relatedness).",
        plot=table.plot(rows, headers=headers),
    )


def add_gene_coverage_section(module, rows):
    """Add a per-sample tabbed gene coverage table to the QC module.

    One Bootstrap tab per sample; each tab has its own searchable + paginated
    table. The function returns silently when `rows` is empty (eg. cohort mode
    which skips this section entirely, or no coverage_by_gene files staged).
    """
    if not rows:
        return

    by_sample = defaultdict(list)
    for rec in rows:
        by_sample[rec["sample"]].append(rec)
    samples = sorted(by_sample)
    if not samples:
        return

    headers = [
        ("gene",             "Gene"),
        ("average_coverage", "Avg. Coverage"),
        ("coverage15",       "% ≥15X"),
        ("coverage30",       "% ≥30X"),
        ("coverage100",      "% ≥100X"),
    ]
    keys = [k for k, _ in headers]

    def fmt(v):
        if v is None or v == "":
            return ""
        try:
            return f"{float(v):,.2f}"
        except (ValueError, TypeError):
            return str(v)

    # nav-tabs header
    nav_html = ['<ul class="nav nav-tabs gc-nav" role="tablist">']
    for i, sample in enumerate(samples):
        active = " active" if i == 0 else ""
        nav_html.append(
            f'<li class="nav-item" role="presentation">'
            f'<a class="nav-link{active}" data-bs-toggle="tab" data-toggle="tab" '
            f'role="tab" href="#gc-pane-{sample}">{sample}</a>'
            f'</li>'
        )
    nav_html.append("</ul>")

    # one tab pane per sample
    pane_html = ['<div class="tab-content gc-tab-content">']
    headers_th = "".join(f"<th>{label}</th>" for _, label in headers)
    for i, sample in enumerate(samples):
        active = " show active" if i == 0 else ""
        tbody_rows = []
        for rec in by_sample[sample]:
            tbody_rows.append(
                "<tr>" + "".join(f"<td>{fmt(rec.get(k))}</td>" if k != 'gene' else f"<td>{rec.get('gene','')}</td>" for k in keys) + "</tr>"
            )
        pane_html.append(
            f'<div class="tab-pane fade{active}" id="gc-pane-{sample}" role="tabpanel">'
            f'<div class="gc-controls" data-sample="{sample}">'
            f'<input type="text" placeholder="Search gene" class="form-control gc-search" style="width:220px;display:inline-block;margin-right:8px;">'
            f'<button class="btn btn-sm btn-outline-secondary gc-prev" style="margin-right:4px;">Prev</button>'
            f'<button class="btn btn-sm btn-outline-secondary gc-next" style="margin-right:8px;">Next</button>'
            f'<span class="gc-info"></span>'
            f'</div>'
            f'<div class="mqc-table-wrapper">'
            f'<table class="mqc-table table table-striped table-condensed gc-table" data-sample="{sample}">'
            f'<thead><tr>{headers_th}</tr></thead>'
            f'<tbody>{"".join(tbody_rows)}</tbody>'
            f'</table></div>'
            f'</div>'
        )
    pane_html.append("</div>")

    # Single script that wires controls for ALL per-sample tables.
    script = """
<style>.gc-controls { margin: 0.5rem 0; } .gc-tab-content { padding-top: 0.75rem; }</style>
<script>
(function(){
  var PAGE_SIZE = 15;
  document.querySelectorAll('.gc-table').forEach(function(table){
    var sample = table.getAttribute('data-sample');
    var pane = table.closest('.tab-pane');
    var controls = pane.querySelector('.gc-controls');
    var searchEl = controls.querySelector('.gc-search');
    var prevEl = controls.querySelector('.gc-prev');
    var nextEl = controls.querySelector('.gc-next');
    var infoEl = controls.querySelector('.gc-info');
    var rows = Array.prototype.slice.call(table.tBodies[0].rows);
    try { rows.sort(function(a,b){ return (a.cells[0].textContent||'').localeCompare(b.cells[0].textContent||''); }); } catch(e){}
    var page = 0;
    function render(){
      var q = (searchEl.value || '').toLowerCase();
      var filt = rows.filter(function(r){ return r.cells[0].textContent.toLowerCase().indexOf(q) > -1; });
      var total = filt.length, pages = Math.max(1, Math.ceil(total/PAGE_SIZE));
      if(page >= pages) page = pages - 1;
      rows.forEach(function(r){ r.style.display = 'none'; });
      var start = page * PAGE_SIZE, end = Math.min(start + PAGE_SIZE, total);
      filt.slice(start, end).forEach(function(r){ r.style.display = 'table-row'; });
      infoEl.textContent = (total === 0 ? '0' : (start+1) + '-' + end) + ' of ' + total;
    }
    searchEl.addEventListener('input', function(){ page = 0; render(); });
    prevEl.addEventListener('click', function(){ if(page>0){ page--; render(); }});
    nextEl.addEventListener('click', function(){ page++; render(); });
    render();
  });
})();
</script>
"""
    module.add_section(
        name="Per-Gene Coverage",
        anchor="gene-coverage",
        description="Coverage metrics per gene; one tab per sample.",
        content="".join(nav_html) + "".join(pane_html) + script,
    )


# ---------------------------------------------------------------------------
# JSON sidecar
# ---------------------------------------------------------------------------

def write_per_sample_json(samples, json_dir, pedigree_pass):
    """
    Write one JSON file per sample matching GA4GH layout.
    """
    json_dir = Path(json_dir)
    json_dir.mkdir(parents=True, exist_ok=True)
    for sample, s in samples.items():
        doc = {
            "biosample": {"id": sample},
            "qc_metrics": {
                "aln_metrics":     {k: getattr(s, k) for k in GA4GH_ALIGNMENT_HEADERS},
                "variant_metrics": {k: getattr(s, k) for k in GA4GH_VARIANT_HEADERS},
            },
            "pedigree": {
                "family_id":              s.family_id,
                "pedigree_sex":           s.pedigree_sex,
                "inferred_sex": s.somalier_sex,
                "sex_check":              s.sex_check,
                "pedigree_validation":    pedigree_pass.get(sample, "na"),
            },
        }
        (json_dir / f"{sample}.metrics.json").write_text(json.dumps(doc, indent=2, default=str))


# ---------------------------------------------------------------------------
# Mode handling
# ---------------------------------------------------------------------------

# Sections of the built-in somalier MultiQC module that are useful in cohort
# mode but noisy in single-family mode.
DROP_SECTIONS_RMFAM = [
    "somalier-relatedness",
    "somalier-sexcheck",
]

def remove_mode_specific_sections(mode):
    """In family mode, drop the built-in somalier sections
    that are replaced by our custom Pedigree section
    """
    if mode != 'family':
        return
    drop = set(DROP_SECTIONS_RMFAM)
    for mod in multiqc.report.modules:
        mod.sections = [s for s in mod.sections if str(s.anchor) not in drop]

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Generate a GA4GH-compliant MultiQC report for the QC pipeline.")
    parser.add_argument("files", nargs="+", help="QC files to parse (staged by Nextflow).")
    parser.add_argument("--config", help="MultiQC config YAML.")
    parser.add_argument("--title", help="Report title.")
    parser.add_argument("--filename", help="Output report filename, e.g. Cohort_multiqc_report.html.")
    parser.add_argument("--mode", choices=["cohort", "family"], default="cohort",
                        help="Report mode. 'family' hides somalier built-in plots; 'cohort' keeps them. Default: cohort.")
    parser.add_argument("--ped", help="PED file driving the canonical sample list. Samples without QC data still appear (na cells).")
    parser.add_argument("--thresholds", help="YAML file overriding default QC thresholds.")
    parser.add_argument("--json-out", help="Directory to write per-sample JSON sidecars.")
    args = parser.parse_args()

    # 1) Load thresholds + optional PED
    thresholds = parse_thresholds(args.thresholds)
    ped_rows = parse_ped(args.ped) if args.ped else []

    # 2) Seed canonical sample list from PED (so every member appears in every table)
    samples: dict[str, SampleMetrics] = {}
    for row in ped_rows:
        s = samples.setdefault(row.sample_id, SampleMetrics())
        s.family_id = row.family_id
        s.pedigree_sex = row.sex

    # 3) Load the MultiQC config and parse all staged files.
    if args.config:
        multiqc.load_config(args.config)
    multiqc.parse_logs(args.files, preserve_module_raw_data=True)

    remove_mode_specific_sections(args.mode)

    # Suppress the built-in "General Statistics" table; our "General Status" replaces it.
    # parse_logs gives modules a chance to contribute via general_stats_addcols(); the
    # only reliable way to keep the resulting table out of the HTML is to clear the data.
    multiqc.report.general_stats_data = {}
    multiqc.report.general_stats_headers = {}

    # 4) Populate samples from each module
    parse_samtools(samples)
    parse_picard(samples)
    parse_verifybamid(samples)
    pairs = parse_somalier(samples, thresholds)

    # 5) Pipeline-specific custom files (VCF metrics, per-gene coverage, DRAGEN)
    vcf_files = [f for f in args.files if f.endswith(("_vcf_metrics.json", ".vcf_metrics.json"))]
    parse_vcf_metrics(vcf_files, samples)
    gene_files = [f for f in args.files if "coverage_by_gene" in f and f.endswith((".tsv", ".csv"))]
    gene_rows = parse_gene_coverage_files(gene_files)
    dragen_files = [f for f in args.files if f.endswith(
        (".mapping_metrics.csv", ".wgs_coverage_metrics.csv",
         ".vc_metrics.csv", ".ploidy_estimation_metrics.csv"))]
    parse_dragen_csv_files(dragen_files, samples)

    # 6) Per-sample pass/fail summary (sex_check is precomputed in parse_somalier)
    pedigree_pass = pedigree_validation_by_sample(pairs)

    # 7) Build the report sections under one "Quality Control" module, prepended
    #    to report.modules so it renders first.
    qc_module = multiqc.BaseMultiqcModule(name="Quality Control", anchor="custom_content")
    add_general_status_section(qc_module, samples, thresholds, pedigree_pass)
    add_metric_table(qc_module, samples, name="Alignment Metrics",
                     anchor="alignment-metrics", headers=GA4GH_ALIGNMENT_HEADERS)
    add_metric_table(qc_module, samples, name="Variant Calling Metrics",
                     anchor="variant-calling-metrics", headers=GA4GH_VARIANT_HEADERS)
    if args.mode == "family":
        # Per-gene coverage and the pedigree validation table are only meaningful per family.
        add_pedigree_section(qc_module, pairs, thresholds)
        add_gene_coverage_section(qc_module, gene_rows)
    multiqc.report.modules = [qc_module] + multiqc.report.modules

    # 8) JSON sidecar + write the HTML report
    if args.json_out:
        write_per_sample_json(samples, args.json_out, pedigree_pass)
    write_kwargs = {"force": True, "title": args.title or None, "zip_data_dir": True}
    if args.filename:
        write_kwargs["filename"] = args.filename
    multiqc.write_report(**write_kwargs)

    log.info("MultiQC report generated successfully.")


if __name__ == "__main__":
    main()
