#!/usr/bin/env python3

import pandas as pd
import os
import sys
import base64

def is_file_empty(file_path):
    return not os.path.exists(file_path) or os.path.getsize(file_path) == 0

def generate_html_report(mosdepth_txt, samtools_txt, plot_png, output_html, median_coverage_csv, region_coverage_csv, other_error_csv):
    if not os.path.exists(output_html):
        os.makedirs(output_html)

    # Determine QC status based on the contents of the CSV files
    if is_file_empty(median_coverage_csv) and is_file_empty(region_coverage_csv) and is_file_empty(other_error_csv):
        qc_status = "QC PASSED"
    else:
        qc_status = "QC FAILED"

    # Encode the PNG file in base64
    with open(plot_png, "rb") as image_file:
        encoded_string = base64.b64encode(image_file.read()).decode('utf-8')

    with open(f"{output_html}/qc_analysis_report.html", 'w') as html_file:
        # Handle mosdepth analysis
        with open(mosdepth_txt, 'r') as file:
            mosdepth_content = file.read().split('\n\n')  # Splitting into paragraphs

        if mosdepth_content:
            html_file.write("<h2>Mosdepth Analysis Summary</h2>\n")
            
            # Write QC status
            color = "green" if qc_status == "QC PASSED" else "red"
            html_file.write(f'<h1 style="color: {color};">{qc_status}</h1>\n')

            for i, paragraph in enumerate(mosdepth_content):
                html_file.write(f"<p>{paragraph}</p>\n")

                # Insert corresponding CSV content after each paragraph, if available
                if i == 1 and not is_file_empty(median_coverage_csv):
                    try:
                        median_df = pd.read_csv(median_coverage_csv, sep='\t')
                        html_file.write(median_df.to_html(index=False))
                    except Exception:
                        html_file.write("<p>No errors.</p>\n")
                elif i == 0 and not is_file_empty(region_coverage_csv):
                    try:
                        region_df = pd.read_csv(region_coverage_csv, sep='\t')
                        html_file.write(region_df.to_html(index=False))
                    except Exception:
                        html_file.write("<p>No errors.</p>\n")
                elif i == 2 and not is_file_empty(other_error_csv):
                    try:
                        other_df = pd.read_csv(other_error_csv, sep='\t')
                        html_file.write(other_df.to_html(index=False))
                    except Exception:
                        html_file.write("<p>No errors.</p>\n")

        # Handle samtools analysis
        with open(samtools_txt, 'r') as file:
            samtools_content = file.read()
        html_file.write(f"<h2>Samtools Analysis</h2>\n<pre>{samtools_content}</pre>")

        # Add the PNG plot as a base64-encoded image
        html_file.write(f"<h2>Samtools Reads Plot</h2>\n")
        html_file.write(f'<img src="data:image/png;base64,{encoded_string}" alt="Samtools Reads Plot" style="max-width:100%; height:auto;">')

    print(f"HTML file '{output_html}/qc_analysis_report.html' has been generated successfully.")


if __name__ == "__main__":
    if len(sys.argv) < 8:
        print("Usage: script.py <mosdepth_txt> <samtools_txt> <plot_png> <output_html> <median_coverage_csv> <region_coverage_csv> <other_error_csv>")
        sys.exit(1)

    mosdepth_txt = sys.argv[1]
    samtools_txt = sys.argv[2]
    plot_png = sys.argv[3]
    output_html = sys.argv[4]

    # Optional CSV files
    median_coverage_csv = sys.argv[5] if len(sys.argv[5].strip()) > 0 else None
    region_coverage_csv = sys.argv[6] if len(sys.argv[6].strip()) > 0 else None
    other_error_csv = sys.argv[7] if len(sys.argv[7].strip()) > 0 else None

    generate_html_report(mosdepth_txt, samtools_txt, plot_png, output_html, median_coverage_csv, region_coverage_csv, other_error_csv)
