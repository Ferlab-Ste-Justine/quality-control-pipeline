#!/usr/bin/env python3

"""
Process BED files to join regions and thresholds data, group by name and compute averages.

This script takes two compressed BED files:
1. regions.bed.gz with format: chr start end name mean
2. thresholds.bed.gz with header and threshold columns

It joins them on chr, start, end, name columns, groups by name,
and computes averages for all numeric columns.
"""
import argparse
import sys
import pandas as pd

__version__ = "1.0.0"

def parse_args():
    parser = argparse.ArgumentParser(description='Join and process BED files to compute averages by region name')

    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument( '--mean', '-m', required=True, help='Input regions BED file with columns: chr start end region mean')
    parser.add_argument('--thresholds', '-t', required=True, help='Input thresholds BED file with header')
    parser.add_argument('--output', '-o', required=True, help='Output file path')

    return parser.parse_args()

def read_bed_file(filepath, has_header=False):
    """Read a compressed BED file and return a pandas DataFrame."""
    try:
        if has_header:
            df = pd.read_csv(filepath, sep='\t', compression='gzip')
            df.columns = df.columns.str.lstrip('#')
        else:
            # For regions file without header, assume format: chr start end name mean
            df = pd.read_csv(
                filepath,
                sep='\t',
                compression='gzip',
                names=['chrom', 'start', 'end', 'region', 'mean'],
                header=None
            )
        return df
    except Exception as e:
        print(f"Error reading {filepath}: {e}", file=sys.stderr)
        sys.exit(1)

def join_and_aggregate(regions_df, thresholds_df):
    """Join the dataframes and compute averages grouped by name."""

    join_columns = ['chrom', 'start', 'end', 'region']

    # Perform inner join
    merged_df = pd.merge(regions_df, thresholds_df, on=join_columns)

    # Calculate region length
    merged_df['size'] = merged_df['end'] - merged_df['start']

    # totalcvg will not be precise because mean has only 2 decimal places. Converting to int at this point so the error does not propagate more.
    merged_df['totalcvg'] = (merged_df['mean'] * merged_df['size']).astype(int)

    # remove mean column
    merged_df = merged_df.drop(columns=['mean'])

    data_columns = [col for col in merged_df.columns if col not in join_columns]

    # Group by name and sum totals
    group_counts = merged_df.groupby('region').size().reset_index(name='count')
    grouped_regions = merged_df.groupby('region',group_keys=True)[data_columns].sum().reset_index()

    # Add count of elements in each group
    grouped = pd.merge(grouped_regions, group_counts, on='region')

    # Calculate mean coverage for regions (totalcvg / size)
    grouped['average_coverage'] = grouped['totalcvg'].div(grouped['size'])

    # Calculate percentages for threshold columns based on region length
    threshold_cols = [col for col in grouped.columns if col.endswith('X')]

    # Calculate proportion covered by at least X.
    for col in threshold_cols:
        pcov_col = f"coverage{col}".rstrip('X')  # e.g., coverage10 for 10X
        grouped[pcov_col] = grouped[col].div(grouped['size'])

    grouped.rename(columns={'region': '#gene'}, inplace=True)
    grouped = grouped.drop(columns=threshold_cols)

    return grouped

def main():
    """Main function."""
    args = parse_args()

    # Read input files
    regions_df = read_bed_file(args.mean, has_header=False)
    thresholds_df = read_bed_file(args.thresholds, has_header=True)

    # Join and aggregate
    result_df = join_and_aggregate(regions_df, thresholds_df)

    if result_df.empty:
        print("No results to output", file=sys.stderr)
        sys.exit(1)

    # Write output
    result_df.to_csv(args.output, sep='\t', index=False)

    print("Processing complete", file=sys.stderr)

if __name__ == "__main__":
    main()
