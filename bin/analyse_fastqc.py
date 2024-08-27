#!/usr/bin/env python3

import os
import sys
import pandas as pd

path_multiqc_general_stats = sys.argv[1]

multiqc_general_stats_df = pd.read_table(
    path_multiqc_general_stats
)

# Create an empty list to store results
comparison_results = []

multiqc_general_stats_df[
    'Sample_Prefix'] = multiqc_general_stats_df['Sample'].apply(
        lambda x: "_".join(x.split("_")[:-1]))

# Group by Sample_Prefix
grouped = multiqc_general_stats_df.groupby('Sample_Prefix')

# Iterate through each group and compare the total sequences
for name, group in grouped:
    if len(group) == 2:  # Ensure there are exactly 2 reads
        read1_total_sequences = group.iloc[0]['FastQC (raw)_mqc-generalstats-fastqc_raw-total_sequences']
        read2_total_sequences = group.iloc[1]['FastQC (raw)_mqc-generalstats-fastqc_raw-total_sequences']
        comparison_results.append({
            'Sample': name,
            'Read1_Total_Sequences': read1_total_sequences,
            'Read2_Total_Sequences': read2_total_sequences,
            'Are_Equal': read1_total_sequences == read2_total_sequences
        })
    else:
        # Handle cases where there's not exactly 2 reads
        comparison_results.append({
            'Sample_Prefix': name,
            'Error': 'Not exactly 2 reads were found'
        })

# Convert the results to a DataFrame for easy viewing
comparison_df = pd.DataFrame(comparison_results)

# Filter the DataFrame to only show rows where Are_Equal is False
unequal_df = comparison_df[comparison_df['Are_Equal'] == False]

if not unequal_df.empty:
    if not os.path.exists('errors'):
        os.makedirs('errors')

    # Define the output file path
    OUTPUT_FILE = "unequal_total_sequences_samples.csv"

    # Write the filtered rows to a CSV file
    unequal_df.to_csv(f"errors/{OUTPUT_FILE}", index=False)

# Save the entire DataFrame to a CSV file
comparison_df.to_csv('fastqc_analysis.csv', index=False)

if unequal_df.empty:
    with open('qc_pass', 'w') as file:
        file.write("All samples have equal total sequences for reads 1 and 2.")
else:
    with open('qc_fail', 'w') as file:
        file.write("Some samples have unequal total sequences for reads 1 and 2.")