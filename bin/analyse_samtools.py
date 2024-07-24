#!/usr/bin/env python3

import statistics
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

reads_mapped = []
sample_names = []
reads_MQ0 = []
reads_unmapped = []
reads_QC_failed_percent = []
flagged_low_reads_samples = []
flagged_high_MQ0_samples = []
flagged_high_unmapped_samples = []
flagged_high_reads_QC_failed_percent = []

samtools_df = pd.read_table("multiqc/multiqc_data/multiqc_samtools_stats.txt", 
                            usecols=['reads_mapped', 'Sample', 'reads_MQ0',
                                     'reads_unmapped',
                                     'reads_QC_failed_percent'],)
samtools_df.fillna(0, inplace=True)
samtools_df.sort_values(by=['reads_mapped'], inplace=True)

reads_mapped = samtools_df.loc[:, 'reads_mapped'].tolist()
sample_names = samtools_df.loc[:, 'Sample'].tolist()
reads_unmapped = samtools_df.loc[:, 'reads_unmapped'].tolist()
reads_MQ0 = samtools_df.loc[:, 'reads_MQ0'].tolist()
reads_QC_failed_percent = samtools_df.loc[:,
                                          'reads_QC_failed_percent'].tolist()

# print(reads_unmapped)
# print(sample_names)

average_reads_mapped = 0
average_reads_MQ0 = 0
average_reads_unmapped = 0
average_reads_QC_failed_percent = 0

if reads_mapped and reads_MQ0 and reads_unmapped and reads_QC_failed_percent:
    average_reads_mapped = sum(reads_mapped) / len(reads_mapped)
    average_reads_MQ0 = sum(reads_MQ0) / len(reads_MQ0)
    average_reads_unmapped = sum(reads_unmapped) / len(reads_unmapped)
    average_reads_QC_failed_percent = sum(reads_QC_failed_percent) / len(
          reads_QC_failed_percent)
    print("The average amount of mapped reads per sample is: " +
          f"{average_reads_mapped}")
    print("The average amount of reads with MQ of 0 per sample is: " +
          f"{average_reads_MQ0}")
    print("The average amount of unmapped reads per sample is: " +
          f"{average_reads_unmapped}")
    print("The average percentage of QC failed reads per sample is: " +
          f"{average_reads_QC_failed_percent}\n")
else:
    print("Missing some valid values to calculate the averages.\n")

min_reads_mapped = float("inf")
min_reads_mapped_index = 0
max_reads_unmapped = 0
max_reads_unmapped_index = 0
max_reads_MQ0 = 0
max_reads_MQ0_index = 0
max_reads_QC_failed_percent = 0
max_reads_QC_failed_percent_index = 0
stdev_reads_mapped = statistics.pstdev(reads_mapped)
stdev_reads_MQ0 = statistics.pstdev(reads_MQ0)
stdev_reads_unmapped = statistics.pstdev(reads_unmapped)
stdev_reads_QC_failed_percent = statistics.pstdev(reads_QC_failed_percent)

for i, read in enumerate(reads_mapped):
    if read < (average_reads_mapped - 2 * stdev_reads_mapped):
        print(f"The number of mapped reads for the sample {sample_names[i]}" +
              " is unusually low")
        flagged_low_reads_samples.append(i)
    if read < min_reads_mapped:
        min_reads_mapped = read
        min_reads_mapped_index = i

for i, mq0 in enumerate(reads_MQ0):
    if mq0 > (average_reads_MQ0 + 3 * stdev_reads_MQ0):
        print(f"The number of MQ0 reads for the sample {sample_names[i]}" +
              " is unusually high")
        flagged_high_MQ0_samples.append(i)
    if mq0 > max_reads_MQ0:
        max_reads_MQ0 = mq0
        max_reads_MQ0_index = i

for i, unmapped in enumerate(reads_unmapped):
    if unmapped > (average_reads_unmapped + (3 * stdev_reads_unmapped)):
        print(f"The number of unmapped reads for the sample {sample_names[i]}"
              + " is unusually high")
        flagged_high_unmapped_samples.append(i)
    if unmapped > max_reads_unmapped:
        max_reads_unmapped = unmapped
        max_reads_unmapped_index = i

for i, QC_failed_percent in enumerate(reads_QC_failed_percent):
    if QC_failed_percent > (average_reads_QC_failed_percent +
                            (3 * stdev_reads_QC_failed_percent)):
        print(f"The number of unmapped reads for the sample {sample_names[i]}"
              + " is unusually high")
        flagged_high_reads_QC_failed_percent.append(i)
    if QC_failed_percent > max_reads_QC_failed_percent:
        max_reads_QC_failed_percent = QC_failed_percent
        max_reads_QC_failed_percent_index = i

print("\nThe minimum amount of mapped reads for a sample is: " +
      str(min_reads_mapped) + " (" +
      str(sample_names[min_reads_mapped_index]) + ")")
print("The maximum amount of MQ0 reads for a sample is: " +
      str(max_reads_MQ0) + " (" +
      str(sample_names[max_reads_MQ0_index]) + ")")
print("The maximum amount of unmapped reads for a sample is: " +
      str(max_reads_unmapped) + " (" +
      str(sample_names[max_reads_unmapped_index]) + ")")
print("The maximum percentage of QC failed reads for a sample sample is: " +
      str(max_reads_QC_failed_percent) + " (" +
      str(sample_names[max_reads_QC_failed_percent_index]) + ")")

# Plotting mapping data
mapped = np.array(reads_mapped)
mq0 = np.array(reads_MQ0)

plt.barh(sample_names, reads_mapped, label='Mapped reads')
plt.barh(sample_names, reads_MQ0,
         left=reads_mapped, color='yellow', label='MQ0 reads')
plt.barh(sample_names, reads_unmapped, left=mapped+mq0, color='red',
         label='Unmapped reads')

min_reads_mapped_index = np.argmin(reads_mapped)
max_reads_mapped_index = np.argmax(reads_mapped)

min_reads_mapped_value = reads_mapped[min_reads_mapped_index]
max_reads_mapped_value = reads_mapped[max_reads_mapped_index]

plt.axvline(min_reads_mapped_value, color='black', linestyle='--',
            label=f'Min mapped reads: {min_reads_mapped_value} ' +
            f'({sample_names[min_reads_mapped_index]})')
plt.axvline(max_reads_mapped_value, color='green', linestyle='--',
            label=f'Max mapped reads: {max_reads_mapped_value} ' +
            f'({sample_names[max_reads_mapped_index]})')

max_unmapped_index = np.argmax(reads_unmapped)
highlight_sample = sample_names[max_unmapped_index]
highlight_band_top = max_unmapped_index + 0.5
highlight_band_bottom = max_unmapped_index - 0.5
plt.axhspan(highlight_band_bottom, highlight_band_top, color='red', alpha=0.3,
            label=f'Most unmapped reads: {highlight_sample}')

plt.title('Représentation graphique du nombre de mapped reads, du nombre\n' +
          'reads avec un mapping quality de 0 et du nombre de unmapped reads')
plt.legend()
plt.show()

# Plotting
