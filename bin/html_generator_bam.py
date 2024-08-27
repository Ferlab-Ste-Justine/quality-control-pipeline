#!/usr/bin/env python3

import pandas as pd
import os
import sys

def generate_html_report(mosdepth_txt, samtools_txt, plot_png, output_html, qc_file=None, median_coverage_csv=None, region_coverage_csv=None, other_error_csv=None):
    with open(output_html, 'w') as html_file:
        # Handle mosdepth analysis
        with open(mosdepth_txt, 'r') as file:
            mosdepth_content = file.read().split('\n\n')  # Splitting into paragraphs

        if mosdepth_content:
            html_file.write("<h2>Mosdepth Analysis Summary</h2>\n")
            
            # Check if QC file exists and write QC PASSED title if it does
            if qc_file and os.path.exists(qc_file):
                html_file.write('<h1 style="color: green;">QC PASSED</h1>\n')
            else:
                html_file.write('<h1 style="color: red;">QC FAILED</h1>\n')
            
            for i, paragraph in enumerate(mosdepth_content):
                html_file.write(f"<p>{paragraph}</p>\n")
                # Insert corresponding CSV content after each paragraph, if available
                if i == 1 and median_coverage_csv and os.path.exists(median_coverage_csv):
                    median_df = pd.read_csv(median_coverage_csv, sep='\t')
                    html_file.write(median_df.to_html(index=False))
                elif i == 0 and region_coverage_csv and os.path.exists(region_coverage_csv):
                    region_df = pd.read_csv(region_coverage_csv, sep='\t')
                    html_file.write(region_df.to_html(index=False))
                elif i == 2 and other_error_csv and os.path.exists(other_error_csv):
                    other_df = pd.read_csv(other_error_csv, sep='\t')
                    html_file.write(other_df.to_html(index=False))

        # Handle samtools analysis
        with open(samtools_txt, 'r') as file:
            samtools_content = file.read()
        html_file.write(f"<h2>Samtools Analysis</h2>\n<pre>{samtools_content}</pre>")

        # Add the PNG plot
        if os.path.exists(plot_png):
            html_file.write(f"<h2>Samtools Reads Plot</h2>\n")
            html_file.write(f"<img src='{plot_png}' alt='Samtools Reads Plot' style='max-width:100%; height:auto;'>")

    print(f"HTML file '{output_html}' has been generated successfully.")


if __name__ == "__main__":
    if len(sys.argv) < 8:
        print("Usage: script.py <mosdepth_txt> <samtools_txt> <plot_png> <qc_file> <output_html> <median_coverage_csv> <region_coverage_csv> <other_error_csv>")
        sys.exit(1)

    mosdepth_txt = sys.argv[1]
    samtools_txt = sys.argv[2]
    plot_png = sys.argv[3]
    output_html = sys.argv[4]
    qc_file = sys.argv[5] if len(sys.argv[5].strip()) > 0 else None

    # Optional CSV files
    median_coverage_csv = sys.argv[6] if len(sys.argv[5].strip()) > 0 else None
    region_coverage_csv = sys.argv[7] if len(sys.argv[6].strip()) > 0 else None
    other_error_csv = sys.argv[8] if len(sys.argv[7].strip()) > 0 else None

    generate_html_report(mosdepth_txt, samtools_txt, plot_png, output_html, qc_file, median_coverage_csv, region_coverage_csv, other_error_csv)
