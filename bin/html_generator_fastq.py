#!/usr/bin/env python3

import csv
import sys
import os


def csv_to_html(results_csv, errors_csv, output_html):
    # Check if the errors file contains a second line
    qc_status = "QC PASSED"
    if os.path.exists(errors_csv):
        with open(errors_csv, mode='r') as errors_file:
            lines = errors_file.readlines()
            if len(lines) > 1:
                qc_status = "QC FAILED"

    # Open the results CSV file for reading
    if os.path.exists(results_csv):
        with open(results_csv, mode='r') as csv_file:
            csv_reader = csv.reader(csv_file)
            headers = next(csv_reader)  # Get the headers from the first row

            if not os.path.exists(output_html):
                os.makedirs(output_html)

            # Start the HTML file
            with open(f"{output_html}/qc_analysis_report.html", mode='w') as html_file:
                html_file.write('<!DOCTYPE html>\n')
                html_file.write('<html>\n')
                html_file.write('<head>\n')
                html_file.write('<title>Analyse des resultats de fastqc</title>\n')
                html_file.write('<style>\n')
                html_file.write('table { border-collapse: collapse; width: 100%; }\n')
                html_file.write('th, td { border: 1px solid black; padding: 8px; text-align: left; }\n')
                html_file.write('th { background-color: #f2f2f2; }\n')
                html_file.write('</style>\n')
                html_file.write('</head>\n')
                html_file.write('<body>\n')

                # Write QC status
                color = "green" if qc_status == "QC PASSED" else "red"
                html_file.write(f'<h1 style="color: {color};">{qc_status}</h1>\n')

                html_file.write('<h1>Résultats de l\'analyse fastqc</h1>\n')
                html_file.write('<table>\n')

                # Write the headers in the table
                html_file.write('<tr>\n')
                for header in headers:
                    html_file.write(f'<th>{header}</th>\n')
                html_file.write('</tr>\n')

                # Write the data rows in the table
                for row in csv_reader:
                    html_file.write('<tr>\n')
                    for column in row:
                        html_file.write(f'<td>{column}</td>\n')
                    html_file.write('</tr>\n')

                # Close the HTML tags for the results section
                html_file.write('</table>\n')

    # Open the errors CSV file for reading
    if os.path.exists(errors_csv):
        with open(errors_csv, mode='r') as csv_file:
            csv_reader = csv.reader(csv_file)
            headers = next(csv_reader)  # Get the headers from the first row
            
            if qc_status == "QC FAILED":
                # Append to the HTML file
                with open(f"{output_html}/qc_analysis_report.html", mode='a') as html_file:
                    html_file.write('<h1>Erreurs trouvées lors de l\'analyse fastqc</h1>\n')
                    html_file.write('<table>\n')

                    # Write the headers in the table
                    html_file.write('<tr>\n')
                    for header in headers:
                        html_file.write(f'<th>{header}</th>\n')
                    html_file.write('</tr>\n')

                    # Write the data rows in the table
                    for row in csv_reader:
                        html_file.write('<tr>\n')
                        for column in row:
                            html_file.write(f'<td>{column}</td>\n')
                        html_file.write('</tr>\n')

                    # Close the HTML tags for the errors section
                    html_file.write('</table>\n')
                    html_file.write('</body>\n')
                    html_file.write('</html>\n')

    print(f"HTML file '{output_html}/qc_analysis_report.html' has been generated successfully.")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: csv_to_html.py <results_csv> <errors_csv> <output_html>")
        sys.exit(1)

    results_csv = sys.argv[1]
    errors_csv = sys.argv[2]
    output_html = sys.argv[3]
    csv_to_html(results_csv, errors_csv, output_html)