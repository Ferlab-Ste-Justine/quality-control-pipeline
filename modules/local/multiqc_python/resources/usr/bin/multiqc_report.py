#!/usr/bin/env python3

"""
MultiQC report generation script for the quality-control-pipeline.
"""

import argparse
import logging
import multiqc
from multiqc import config
from multiqc.plots import table
import os
import glob
from collections import defaultdict
import json

# Initialise the logger
log = logging.getLogger('multiqc')

def get_alignment_and_pedigree_data():
    """Retrieve alignment and pedigree data from MultiQC's parsed data.

    Returns (alignment_stats_data, pedigree_data, somalier_sex_map)
    """
    alignment_stats_data = defaultdict(dict)
    pedigree_data = []
    somalier_sex_map = {}

    # Data from samtools stats (module may be missing)
    try:
        samtools_data = multiqc.get_module_data(module='samtools') or {}
    except ValueError:
        samtools_data = {}
    if isinstance(samtools_data, dict):
        for s_name, data in samtools_data.items():
            if data is None:
                continue
            if not any(k in data for k in ('raw_total_sequences', 'reads_mapped', 'reads_properly_paired_percent', 'reads_mapped_percent')):
                continue
            alignment_stats_data[s_name].update({
                'total_reads': data.get('raw_total_sequences'),
                'mapped_reads': data.get('reads_mapped'),
                'duplicate_reads': data.get('reads_duplicated'),
                'pct_reads_properly_paired': data.get('reads_properly_paired_percent'),
                'mean_insert_size': data.get('insert_size_average'),
                'insert_size_std_deviation': data.get('insert_size_standard_deviation'),
                'pct_reads_mapped': data.get('reads_mapped_percent'),
            })

    # Data from picard CollectWgsMetrics (module may be missing)
    try:
        picard_data = multiqc.get_module_data(module='picard') or {}
    except ValueError:
        picard_data = {}
    if isinstance(picard_data, dict):
        for key, pdata in picard_data.items():
            if pdata is None:
                continue
            for s_name, data in pdata.items():
                s_name = _clean_sample_name(s_name)
                if not isinstance(data, dict) or 'MEAN_COVERAGE' not in data:
                    continue
                mc = data.get('MEAN_COVERAGE')
                mc = float(mc) if mc is not None else None
                sd = data.get('SD_COVERAGE')
                sd = float(sd) if sd is not None else None
                mad = data.get('MAD_COVERAGE')
                mad = float(mad) if mad is not None else None
                p15 = data.get('PCT_15X') or data.get('PCT_15X')
                p15 = float(p15) if p15 is not None else None
                alignment_stats_data[s_name].update({
                    'mean_autosome_coverage': mc,
                    'mad_autosome_coverage': mad,
                    'pct_autosomes_15x': p15,
                })

    # Data from verifybamid (module may be missing)
    try:
        verify_data = multiqc.get_module_data(module='verifybamid') or {}
    except ValueError:
        verify_data = {}
    if isinstance(verify_data, dict):
        for s_name, data in verify_data.items():
            s_name = _clean_sample_name(s_name)
            if data is None:
                continue
            if 'FREEMIX' not in data:
                continue
            alignment_stats_data[s_name].update({
                'cross_contamination_rate': data.get('FREEMIX'),
            })

    # Data from somalier
    # somalier module returns a mapping of sample -> data or a list for relatedness
    try:
        raw_somalier = multiqc.get_module_data(module='somalier')
    except ValueError:
        raw_somalier = {}
    # normalize pedigree_data to list of dicts with sample_a/sample_b/relatedness
    pedigree_data = []
    if isinstance(raw_somalier, dict):
        # somalier module may have per-sample entries and/or relatedness entries
        for key, val in raw_somalier.items():
            if isinstance(val, dict):
                # collect sex predictions if present (per-sample keys)
                for sex_key in ('sex', 'sex_call', 'predicted_sex'):
                    if sex_key in val:
                        somalier_sex_map[key] = val.get(sex_key)
                        break
                # relatedness entries: often stored under a composite key like "SAMPLEA*SAMPLEB"
                if 'relatedness' in val:
                    # try to parse sample names from the parent key if present
                    sa = None
                    sb = None
                    if isinstance(key, str) and '*' in key:
                        parts = key.split('*', 1)
                        sa, sb = parts[0], parts[1]
                    merged = dict(val)
                    if sa:
                        merged['sample_a'] = sa
                    if sb:
                        merged['sample_b'] = sb
                    pedigree_data.append(merged)
                # sometimes relatedness list under 'relate'
                if 'relate' in val and isinstance(val.get('relate'), list):
                    for r in val.get('relate'):
                        if isinstance(r, dict):
                            pedigree_data.append(r)
    elif isinstance(raw_somalier, list):
        pedigree_data = raw_somalier
    else:
        pedigree_data = []

    return alignment_stats_data, pedigree_data, somalier_sex_map

def _clean_sample_name(path):
    bn = os.path.basename(path)
    if '.' in bn:
        bn = bn.split('.', 1)[0]
    return bn.strip('._-')

def parse_gene_coverage(f, data):
    """Parse gene coverage output."""
    s_name = _clean_sample_name(f)
    with open(f, 'r') as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            fields = line.split(',')
            try:
                data.append({
                    'sample': s_name,
                    'gene': fields[0],
                    'average_coverage': float(fields[2]),
                    'coverage15': float(fields[4]),
                    'coverage30': float(fields[5]),
                    'coverage100': float(fields[7]),
                })
            except (IndexError, ValueError) as e:
                log.warning(f"Could not parse line in {f}: {line} - {e}")

def add_general_status_section(module, alignment_stats_data, qc_thresholds, pedigree_data, somalier_sex_map, variant_metrics_map=None):
    """Add the General Status section to the report."""
    # Build a union of samples present in any data source so we show whichever info is available
    samples = set()
    samples.update(alignment_stats_data.keys())
    samples.update((rec.get('sample') for rec in (pedigree_data or []) if rec.get('sample')))
    samples.update(somalier_sex_map.keys())
    if variant_metrics_map:
        samples.update(variant_metrics_map.keys())

    # build pedigree pass set based solely on somalier's expected_relatedness field
    pedigree_pass_set = set()
    pedigree_compared_set = set()
    for rec in (pedigree_data or []):
        rel_raw = rec.get('relatedness')
        exp_raw = rec.get('expected_relatedness')
        try:
            if rel_raw is None or exp_raw is None:
                continue
            rel_f = float(rel_raw)
            exp_f = float(exp_raw)
        except Exception:
            # skip records that cannot be interpreted as numbers
            continue
        # mark that these samples had a comparison performed
        sa = rec.get('sample_a') or rec.get('s_name') or rec.get('sample') or rec.get('sample_a_name')
        sb = rec.get('sample_b') or rec.get('other_sample') or rec.get('sample_b_name')
        if sa:
            pedigree_compared_set.add(sa)
        if sb:
            pedigree_compared_set.add(sb)
        # pass if observed relatedness meets or exceeds expected relatedness
        if rel_f >= exp_f:
            if sa:
                pedigree_pass_set.add(sa)
            if sb:
                pedigree_pass_set.add(sb)

    status_data = {}
    for s_name in sorted(samples):
        data = alignment_stats_data.get(s_name, {})
        pct_mapped = data.get('pct_reads_mapped') if data else None
        pct_proper = data.get('pct_reads_properly_paired') if data else None
        contam = data.get('cross_contamination_rate') if data else None
        # Preserve None for missing mean coverage (don't default to 0.0X)
        mean_cov = data.get('mean_autosome_coverage') if data else None

        # aln quality: pass/fail/na
        if pct_mapped is None or pct_proper is None:
            aln_quality = 'na'
        else:
            aln_quality = 'pass' if (pct_mapped >= qc_thresholds['pct_reads_mapped'] and pct_proper >= qc_thresholds['pct_reads_properly_paired']) else 'fail'

        # contamination
        if contam is None:
            contamination = 'na'
        else:
            contamination = 'pass' if (contam <= qc_thresholds['cross_contamination_rate']) else 'fail'

        # vcf quality: present -> pass if total_snvs>0 else fail; absent -> na
        vcf_quality = None
        if variant_metrics_map is not None:
            vm = variant_metrics_map.get(s_name)
            if vm is None:
                vcf_quality = 'na'
            else:
                total_snvs = vm.get('total_snvs')
                vcf_quality = 'pass' if total_snvs > 0 else 'fail'

        # If there were pedigree comparisons available for this sample, use pass/fail.
        # Otherwise mark as 'na' (no expected_relatedness to compare).
        if s_name in pedigree_pass_set:
            pedigree_result = 'pass'
        elif s_name in pedigree_compared_set:
            pedigree_result = 'fail'
        else:
            pedigree_result = 'na'
        sex_result = 'pass' if s_name in (somalier_sex_map or {}) else ('na' if not somalier_sex_map else 'fail')

        status_data[s_name] = {
            'aln_quality': aln_quality,
            'contamination': contamination,
            'vcf_quality': vcf_quality,
            'sex_check': sex_result,
            'pedigree_validation': pedigree_result,
            'mean_coverage': mean_cov
        }

    # Build headers dynamically depending on which columns have any data
    status_headers = {
        'aln_quality': {'title': 'Aln Quality', 'description': 'Alignment quality (mapped/properly paired)', 'scale': 'pass-fail'},
        'contamination': {'title': 'Contamination', 'description': f"Contamination check (FREEMIX <= {qc_thresholds['cross_contamination_rate']})", 'scale': 'pass-fail'},
        'sex_check': {'title': 'Sex Check', 'description': 'Sex verification', 'scale': 'pass-fail'},
        'pedigree_validation': {'title': 'Pedigree Validation', 'description': 'Pedigree consistency', 'scale': 'pass-fail'},
        'mean_coverage': {'title': 'Mean Coverage', 'description': 'Mean autosome coverage', 'format': '{:,.2f}X'}
    }
    # Add vcf header only if variant metrics were provided
    if variant_metrics_map is not None:
        status_headers['vcf_quality'] = {'title': 'VCF Quality', 'description': 'Variant calling presence/quality', 'scale': 'pass-fail'}

    module.add_section(
        name='General Status',
        anchor='general-status',
        description='This section provides a summary of the quality control metrics.',
        plot=table.plot(status_data, headers=status_headers)
    )

def add_pedigree_section(module, pedigree_data, qc_thresholds):
    """Add the Pedigree section to the report.

    Expects `pedigree_data` as a list of dicts with keys like
    `sample_a`, `sample_b`, `relatedness`, `ibs0`, `ibs2`.
    Adds a final Pass/Fail column based on `qc_thresholds['relatedness_threshold']`.
    """
    if not pedigree_data:
        return

    # Build rows with sample columns and a pass/fail result based only on
    # somalier's `expected_relatedness` field (if present).
    rows = {}
    for rec in pedigree_data:
        sa = rec.get('sample_a') or rec.get('s_name') or rec.get('sample')
        sb = rec.get('sample_b') or rec.get('other_sample')
        # fallback: if neither present, try parsing a composite key in rec
        if not sa and not sb:
            pair = rec.get('pair') or rec.get('sample_pair')
            if pair and isinstance(pair, str) and '*' in pair:
                try:
                    sa, sb = pair.split('*', 1)
                except Exception:
                    sa, sb = (None, None)

        # parse observed and expected relatedness
        rel = None
        exp_rel = None
        try:
            rel_raw = rec.get('relatedness')
            if rel_raw is not None:
                rel = float(rel_raw)
        except Exception:
            rel = None
        try:
            exp_raw = rec.get('expected_relatedness')
            if exp_raw is not None:
                exp_rel = float(exp_raw)
        except Exception:
            exp_rel = None

        ibs0 = rec.get('ibs0')
        ibs2 = rec.get('ibs2')

        # Determine pass/fail: only if expected_relatedness present can we compare.
        if rel is None or exp_rel is None:
            passfail = 'na'
        else:
            passfail = 'pass' if rel >= exp_rel else 'fail'

        sa_disp = sa or 'unknown'
        sb_disp = sb or 'unknown'
        key = f"{sa_disp} | {sb_disp}"
        rows[key] = {
            'sample_a': sa,
            'sample_b': sb,
            'relatedness': rel,
            'expected_relatedness': exp_rel,
            'ibs0': ibs0,
            'ibs2': ibs2,
            'passfail': passfail
        }

    pedigree_headers = {
        'sample_a': {'title': 'Sample A', 'description': 'Sample A'},
        'sample_b': {'title': 'Sample B', 'description': 'Sample B'},
        'relatedness': {'title': 'Relatedness', 'description': 'Observed relatedness coefficient', 'format': '{:,.4f}'},
        'expected_relatedness': {'title': 'Expected Relatedness', 'description': 'Expected relatedness from somalier', 'format': '{:,.4f}'},
        'ibs0': {'title': 'IBS0', 'description': 'Number of sites with 0 shared alleles', 'format': '{:,.0f}'},
        'ibs2': {'title': 'IBS2', 'description': 'Number of sites with 2 shared alleles', 'format': '{:,.0f}'},
        'passfail': {'title': 'Pass/Fail', 'description': 'Pedigree validation result based on expected_relatedness'}
    }

    # module = multiqc.BaseMultiqcModule(
    #     name="Pedigree",
    #     anchor="custom_content",
    # )

    module.add_section(
        name="Pedigree",
        anchor="pedigree",
        description="This section shows the pedigree information from Somalier.",
        plot=table.plot(rows, headers=pedigree_headers)
    )

    # multiqc.report.modules = [module] + multiqc.report.modules

def add_gene_coverage_section(gene_coverage_data):
    """Add the Gene Coverage section to the report (searchable table)."""
    if not gene_coverage_data:
        return
    gene_coverage_headers = {
        'sample': {'title': 'Sample', 'description': 'Sample name'},
        'gene': {'title': 'Gene', 'description': 'Gene name'},
        'average_coverage': {'title': 'Avg. Coverage', 'description': 'Average coverage', 'format': '{:,.2f}X'},
        'coverage15': {'title': '% >= 15X', 'description': 'Percentage of bases with >= 15X coverage', 'format': '{:,.2f}%'},
        'coverage30': {'title': '% >= 30X', 'description': 'Percentage of bases with >= 30X coverage', 'format': '{:,.2f}%'},
        'coverage100': {'title': '% >= 100X', 'description': 'Percentage of bases with >= 100X coverage', 'format': '{:,.2f}%'},
    }

    module_genes = multiqc.BaseMultiqcModule(
        name="Per-Gene Coverage",
        anchor="custom_content",
    )

    # Convert list of dict rows into a dict keyed by sample_gene index so table.plot accepts it
    rows = {}
    for i, rec in enumerate(gene_coverage_data):
        sample = rec.get('sample') or f'sample_{i}'
        gene = rec.get('gene') or f'gene_{i}'
        key = f"{sample}__{gene}__{i}"
        rows[key] = {
            'sample': sample,
            'gene': gene,
            'average_coverage': rec.get('average_coverage'),
            'coverage15': rec.get('coverage15'),
            'coverage30': rec.get('coverage30'),
            'coverage100': rec.get('coverage100'),
        }

    # Build an HTML table so custom JS can search by gene (avoid automatic violin rendering)
    # Table columns are: Sample, Gene, Avg Coverage, %>=15X, %>=30X, %>=50X, %>=100X
    cols = ['sample', 'gene', 'average_coverage', 'coverage15', 'coverage30', 'coverage100']
    headers_html = ''.join(f"<th>{gene_coverage_headers[c]['title']}</th>" for c in cols)
    def _fmt(v):
        if v is None or v == '':
            return ''
        try:
            return "{:,.2f}".format(float(v))
        except Exception:
            return str(v)

    rows_html = []
    for key, rec in rows.items():
        cells = [rec.get('sample', ''), rec.get('gene', ''),
                 _fmt(rec.get('average_coverage', '')),
                 _fmt(rec.get('coverage15', '')),
                 _fmt(rec.get('coverage30', '')),
                 _fmt(rec.get('coverage100', ''))]
        row_html = '<tr>' + ''.join(f"<td>{cell}</td>" for cell in cells) + '</tr>'
        rows_html.append(row_html)

    table_html = (
        '<div class="mqc-table-wrapper">'
        + '<table id="gene_coverage_table" class="mqc-table table table-striped table-condensed">'
        + '<thead><tr>' + headers_html + '</tr></thead>'
        + '<tbody>' + '\n'.join(rows_html) + '</tbody>'
        + '</table></div>'
    )

    # Local JS/CSS for search, filter and pagination to avoid external CDNs
    local_table_script = '''
<style>
.mqc-gene-controls { margin-bottom: 0.5rem; }
.mqc-gene-controls .form-control { display: inline-block; vertical-align: middle; }
</style>
<script>
(function(){
    function $(sel, ctx){ return (ctx||document).querySelector(sel); }
    var table = document.getElementById('gene_coverage_table');
    if(!table) return;
    var tbody = table.tBodies[0];
    var rows = Array.prototype.slice.call(tbody.rows);
    // Pre-sort rows by gene name (column index 1) once to enable fast, single-key sorting
    try{
        rows.sort(function(a,b){
            var A = a.cells[1].textContent || '';
            var B = b.cells[1].textContent || '';
            return A.localeCompare(B);
        });
    }catch(e){/* ignore sort errors */}
    var pageSize = 15;
    var currentPage = 0;

    function render(){
        var queryEl = document.getElementById('gene-search');
        var query = queryEl ? queryEl.value.toLowerCase() : '';
        var filtered = rows.filter(function(r){
            return r.cells[1].textContent.toLowerCase().indexOf(query) > -1;
        });
        var total = filtered.length;
        var totalPages = Math.max(1, Math.ceil(total / pageSize));
        if(currentPage >= totalPages) currentPage = totalPages - 1;
        rows.forEach(function(r){ r.style.display = 'none'; });
        var start = currentPage * pageSize;
        var end = Math.min(start + pageSize, total);
        filtered.slice(start, end).forEach(function(r){ r.style.display = 'table-row'; });
        var info = document.getElementById('gc_pagination_info');
        if(info) info.textContent = (total === 0 ? '0' : (start+1) + '-' + end) + ' of ' + total;
    }

    // build controls
    var controls = document.createElement('div');
    controls.className = 'mqc-gene-controls';
    controls.innerHTML = '<input id="gene-search" placeholder="Search gene" class="form-control" style="width:220px;margin-right:8px;">'
        + '<select id="gene-page-size" class="form-control" style="width:90px;margin-right:8px;"><option>10</option><option>25</option><option>50</option><option>100</option></select>'
        + '<button id="gene-prev" class="btn btn-sm btn-outline-secondary" style="margin-right:4px;">Prev</button>'
        + '<button id="gene-next" class="btn btn-sm btn-outline-secondary" style="margin-right:8px;">Next</button>'
        + '<span id="gc_pagination_info"></span>';
    table.parentNode.insertBefore(controls, table);

    document.getElementById('gene-search').addEventListener('input', function(){ currentPage = 0; render(); });
    document.getElementById('gene-page-size').addEventListener('change', function(){ pageSize = parseInt(this.value, 10); currentPage = 0; render(); });
    document.getElementById('gene-prev').addEventListener('click', function(){ if(currentPage>0){ currentPage--; render(); }});
    document.getElementById('gene-next').addEventListener('click', function(){ currentPage++; render(); });
    document.getElementById('gene-page-size').value = pageSize;
    render();
})();
</script>
'''

    content = table_html + local_table_script

    module_genes.add_section(
        name="Per-Gene Coverage",
        anchor="gene-coverage",
        description="This section shows the coverage metrics for each gene.",
        content=content
    )

    multiqc.report.modules = [module_genes] + multiqc.report.modules

def parse_vcf_metrics(f, data):
    """Parse a per-sample VCF metrics JSON file into `data` list."""
    try:
        with open(f, 'r') as fh:
            j = json.load(fh)
            # Expect a top-level dict with a `sample_id` and metric keys
            if isinstance(j, dict) and 'sample_id' in j:
                j['_source_file'] = f
                data.append(j)
    except Exception as e:
        log.warning("Could not parse VCF metrics file %s: %s", f, e)

def add_qc_metrics_module(module, alignment_stats_data, variant_metrics_data):
    """Create a single module that contains Alignment Statistics and Variant Calling Metrics sections."""
    if not alignment_stats_data and not variant_metrics_data:
        return

    # Alignment section
    if alignment_stats_data:
        alignment_headers = {
            'total_reads': {'title': 'Total Reads', 'description': 'Total number of reads.', 'format': '{:,.2f}'},
            'mapped_reads': {'title': 'Mapped Reads', 'description': 'Number of mapped reads.', 'format': '{:,.2f}'},
            'pct_reads_mapped': {'title': '% Mapped', 'description': 'Percentage of mapped reads.', 'format': '{:,.2f}%'},
            'duplicate_reads': {'title': 'Duplicate Reads', 'description': 'Number of duplicate reads.', 'format': '{:,.2f}'},
            'pct_reads_properly_paired': {'title': '% Properly Paired', 'description': 'Percentage of properly paired reads.', 'format': '{:,.2f}%'},
            'mean_insert_size': {'title': 'Mean Insert Size', 'description': 'Mean insert size.', 'format': '{:,.2f}'},
            'insert_size_std_deviation': {'title': 'Insert Size SD', 'description': 'Insert size standard deviation.', 'format': '{:,.2f}'},
            'mean_autosome_coverage': {'title': 'Mean Autosome Coverage', 'description': 'Mean autosome coverage.', 'format': '{:,.2f}X'},
            'mad_autosome_coverage': {'title': 'MAD Autosome Coverage', 'description': 'Median Absolute Deviation of autosome coverage.', 'format': '{:,.2f}'},
            'pct_autosomes_15x': {'title': '% Autosomes >= 15X', 'description': 'Percentage of autosomes with coverage >= 15X.', 'format': '{:,.2f}%'},
            'cross_contamination_rate': {'title': 'Contamination Rate', 'description': 'Cross-contamination rate (FREEMIX).', 'format': '{:,.2f}'}
        }

        module.add_section(
            name='Alignment Metrics',
            anchor='alignment-metrics',
            description="This section summarizes the alignment metrics for each sample.",
            plot=table.plot(alignment_stats_data, headers=alignment_headers)
        )

    # Variant section
    if variant_metrics_data:
        # build rows similar to add_variant_metrics_section
        fields = [
            ('total_snvs', 'Total SNVs'),
            ('heterozygous_snvs', 'Het SNVs'),
            ('homozygous_snvs', 'Hom SNVs'),
            ('total_insertions', 'Insertions'),
            ('total_deletions', 'Deletions'),
            ('heterozygous_indels', 'Het indels'),
            ('homozygous_indels', 'Hom indels'),
            ('het_hom_ratio_snvs', 'Het/Hom SNV'),
            ('het_hom_ratio_indels', 'Het/Hom indel'),
            ('ins_dels_ratio', 'Ins/Dels'),
            ('ti_tv_ratio', 'Ti/Tv')
        ]

        rows = {}
        for rec in variant_metrics_data:
            sid = rec.get('sample_id') or rec.get('sample') or os.path.splitext(os.path.basename(rec.get('_source_file', '')))[0]
            key = sid
            row = {'sample': sid}
            for k, _ in fields:
                row[k] = rec.get(k)
            rows[key] = row

        headers = {'sample': {'title': 'Sample', 'description': 'Sample ID'}}
        for k, title in fields:
            headers[k] = {'title': title, 'description': title, 'format': '{:,.2f}'}

        module.add_section(
            name='Variant Calling Metrics',
            anchor='variant-calling-metrics',
            description='Per-sample variant calling summary metrics (VCF).',
            plot=table.plot(rows, headers=headers)
        )

def main():
    """ Main entry point """
    parser = argparse.ArgumentParser(description="Generate a custom MultiQC report.")
    parser.add_argument('files', nargs='+', help="List of QC files to parse")
    parser.add_argument('--config', help="MultiQC config file to use.")
    parser.add_argument('--outdir', default='.', help="Output directory for the report.")
    args = parser.parse_args()

    other_files = list(args.files) #[f for f in args.files if f != parquet_file]

    config_file = args.config
    multiqc.parse_logs(args.files, preserve_module_raw_data=True, config_files=[config_file] if config_file else None)

    # Get standard module data
    alignment_stats_data, pedigree_data, somalier_sex_map = get_alignment_and_pedigree_data()

    # Find and parse custom gene coverage files
    gene_coverage_data = []
    coverage_files = [f for f in other_files if 'coverage_by_gene' in f and f.endswith('.csv')]
    if not coverage_files:
        # auto-discover coverage CSVs in the repo under the data/ tree
        coverage_files = glob.glob(os.path.join('data', '**', '*coverage_by_gene*.csv'), recursive=True)
    for f in coverage_files:
        parse_gene_coverage(f, gene_coverage_data)

    # Find and parse VCF metrics files
    variant_metrics_data = []
    vcf_files = [f for f in other_files if f.endswith('.vcf_metrics.json') or f.endswith('_vcf_metrics.json')]
    if not vcf_files:
        vcf_files = glob.glob(os.path.join('data', '**', '*vcf_metrics.json'), recursive=True)
    for f in vcf_files:
        parse_vcf_metrics(f, variant_metrics_data)

    # Set QC thresholds
    qc_thresholds = {
        'pct_reads_mapped': 95,
        'pct_reads_properly_paired': 95,
        'mean_autosome_coverage': 25,
        'cross_contamination_rate': 0.02
    }
    # Allow overriding from MultiQC config (use robust fallback for different MultiQC versions)
    try:
        if hasattr(multiqc, 'get_config'):
            qc_thresholds = multiqc.get_config('qc_thresholds', qc_thresholds)
        elif hasattr(config, 'get_config'):
            qc_thresholds = config.get_config('qc_thresholds', qc_thresholds)
        elif hasattr(config, 'Set') and hasattr(config.Set, 'get'):
            try:
                qc_thresholds = config.Set.get('qc_thresholds', qc_thresholds)
            except Exception:
                pass
    except Exception:
        # If any unexpected error occurs, keep the defaults
        log.debug('Could not read qc_thresholds from MultiQC config, using defaults')

    module = multiqc.BaseMultiqcModule(
        name="Quality Control",
        anchor="custom_content",
    )
    # Add custom sections
    # build variant metrics map for general status (sample_id -> metrics)
    variant_metrics_map = {rec.get('sample_id'): rec for rec in variant_metrics_data if isinstance(rec, dict) and rec.get('sample_id')}
    add_general_status_section(module, alignment_stats_data, qc_thresholds, pedigree_data, somalier_sex_map, variant_metrics_map if variant_metrics_map else None)
    add_qc_metrics_module(module, alignment_stats_data, variant_metrics_data)
    add_pedigree_section(module, pedigree_data, qc_thresholds)
    add_gene_coverage_section(gene_coverage_data)

    multiqc.report.modules = [module] + multiqc.report.modules

    # Write the final report
    multiqc.write_report(force=True,
        exclude_modules=["general_stats"])

    log.info("Custom MultiQC report generated successfully!")

if __name__ == "__main__":
    main()
