#!/usr/bin/env python3

import csv
import sys
import os


def csv_to_html(results_csv, errors_csv, output_html, qc_file=None):
    # Open the CSV file for reading
    if os.path.exists(results_csv):
        with open(results_csv, mode='r') as csv_file:
            csv_reader = csv.reader(csv_file)
            headers = next(csv_reader)  # Get the headers from the first row

            # Start the HTML file
            with open(output_html, mode='w') as html_file:
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

                # Check if QC file exists and write QC PASSED title if it does
                if qc_file and os.path.exists(qc_file):
                    html_file.write('<h1 style="color: green;">QC PASSED</h1>\n')
                else:
                    html_file.write('<h1 style="color: red;">QC FAILED</h1>\n')

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

                # Close the HTML tags
                html_file.write('</table>\n')

    # Open the CSV file for reading
    if os.path.exists(results_csv):
        with open(errors_csv, mode='r') as csv_file:
            csv_reader = csv.reader(csv_file)
            headers = next(csv_reader)  # Get the headers from the first row

            # Start the HTML file
            with open(output_html, mode='a') as html_file:

                html_file.write('<h1>Erreurs trouvees lors de l\'analyse fastqc</h1>\n')
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

                # Close the HTML tags
                html_file.write('</table>\n')
                html_file.write('</body>\n')
                html_file.write('</html>\n')

    print(f"HTML file '{output_html}' has been generated successfully.")


if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: csv_to_html.py <results_csv> <errors_csv> <output_html> <qc_file>")
        sys.exit(1)

    results_csv = sys.argv[1]
    errors_csv = sys.argv[2]
    output_html = sys.argv[3]
    qc_file = sys.argv[4] if len(sys.argv) == 4 else None
    csv_to_html(results_csv, errors_csv, output_html, qc_file)
