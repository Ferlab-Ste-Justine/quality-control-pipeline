import pandas as pd
import sys

couverture_region_minimale = float(sys.argv[1])
couverture_mediane_minimale = float(sys.argv[2])
nombre_std = float(sys.argv[3])
couverture_minimale_normale = float(sys.argv[4])
path_coverage_per_contig = sys.argv[5]
path_multiqc_general_stats = sys.argv[6]


mosdepth_per_region_df = pd.read_table(
    path_coverage_per_contig,
    usecols=['Sample', 'chr1', 'chr2',
             'chr3', 'chr4', 'chr5', 'chr6',
             'chr7', 'chr8', 'chr9',
             'chr10', 'chr11', 'chr12',
             'chr13', 'chr14', 'chr15',
             'chr16', 'chr17', 'chr18',
             'chr19', 'chr20', 'chr21',
             'chr22', 'chrX', 'chrY'],)
mosdepth_per_region_df.fillna(0, inplace=True)

problems = []

if mosdepth_per_region_df.isin([0]).any().any():
    print("ERREUR : probleme a cause de la courverture par region pour les" +
          " echantillons:")

    for index, value in mosdepth_per_region_df['Sample'].items():
        if 0 in mosdepth_per_region_df.loc[index].values:
            print(f"{mosdepth_per_region_df.loc[index, 'Sample']}")
            for column in mosdepth_per_region_df.columns:
                if column == 'Sample':
                    continue
                try:
                    cell_value = mosdepth_per_region_df.loc[index, column]
                    if cell_value is None or cell_value == "None" or cell_value == '':
                        cell_value = 0.0
                    if float(cell_value) <= couverture_region_minimale:
                        problems.append([mosdepth_per_region_df.loc[index,
                                                                    'Sample'],
                                         column, cell_value])
                except (ValueError, TypeError) as e:
                    problems.append([mosdepth_per_region_df.loc[
                        index, 'Sample'], column,
                                     f"Erreur: {e}"])

    problems_df = pd.DataFrame(problems, columns=['Sample', 'Chromosome',
                                                  'Coverage'])

    # Sort by Sample and Chromosome, handling chrX and chrY separately
    problems_df['Chromosome'] = problems_df['Chromosome'].apply(
        lambda x: (23 if x == 'chrX' else (24 if x == 'chrY' 
                                           else int(x.replace('chr', '')))))
    problems_df.sort_values(by=['Sample', 'Chromosome'], inplace=True)
    problems_df['Chromosome'] = problems_df['Chromosome'].apply(
        lambda x: ('chrX' if x == 23 else ('chrY' if x == 24 else f'chr{x}')))

    print("Pour plus de détails sur les régions manquantes," +
          " veuillez consulter le fichier de sortie " +
          "\'erreurs_de_couvertures_par_region.txt\'\n")

    try:
        # Ouvrir le fichier en mode écriture
        with open('erreurs_de_couvertures_par_region.txt', 'w') as f:
            f.write("Sample\tChromosome\tCoverage\n")
            for _, row in problems_df.iterrows():
                f.write(f"{row['Sample']}\t{row['Chromosome']}\t" +
                        f"{row['Coverage']}\n")
    except IOError as e:
        print(f"Erreur lors de l'ouverture ou de l'écriture du fichier: {e}")
    except Exception as e:
        print(f"Une erreur inattendue est survenue: {e}")
else:
    print("Toutes les régions ont ete correctement couvertes.")

mosdepth_general_df = pd.read_table(path_multiqc_general_stats,
                                    usecols=['Sample',
                                             'Mosdepth_mqc-generalstats-' +
                                             'mosdepth-median_coverage',
                                             'Mosdepth_mqc-generalstats-' +
                                             'mosdepth-1_x_pc', 'Mosdepth_mqc'
                                             + '-generalstats-mosdepth-5_x_pc',
                                             'Mosdepth_mqc-generalstats-' +
                                             'mosdepth-10_x_pc',
                                             'Mosdepth_mqc-generalstats-' +
                                             'mosdepth-30_x_pc',
                                             'Mosdepth_mqc-generalstats-' +
                                             'mosdepth-50_x_pc'],)
mosdepth_general_df.fillna(0, inplace=True)

if mosdepth_general_df['Mosdepth_mqc-generalstats-' +
                       'mosdepth-median_coverage'].isin([0]).any():
    print("ERREUR : probleme a cause de la courverture médianne pour les" +
          " echantillons:")

    problems = []

    for index, row in mosdepth_general_df.iterrows():
        try:
            cell_value = row['Mosdepth_mqc-generalstats-mosdepth-' +
                             'median_coverage']
            if cell_value is None or cell_value == "None" or cell_value == '':
                cell_value = 0.0
            if float(cell_value) <= couverture_mediane_minimale:
                sample_name = row['Sample']
                print(f"Sample: {sample_name}")
                problems.append([sample_name, cell_value])
        except (ValueError, TypeError) as e:
            problems.append([row['Sample'], f"Erreur: {e}"])

    # Create a DataFrame for problems
    problems_df = pd.DataFrame(problems, columns=['Sample', 'Median_coverage'])

    print("Pour plus de détails sur les régions manquantes, veuillez" +
          " consulter le fichier de sortie " +
          "'erreurs_de_couvertures_mediane.txt'")

    try:
        # Write the details to the output file
        with open('erreurs_de_couvertures_mediane.txt', 'w') as f:
            f.write("Sample\tMedian_coverage\n")
            for _, row in problems_df.iterrows():
                f.write(f"{row['Sample']}\t{row['Median_coverage']}\n")
    except IOError as e:
        print(f"Erreur lors de l'ouverture ou de l'écriture du fichier: {e}")
    except Exception as e:
        print(f"Une erreur inattendue est survenue: {e}")
else:
    print("Tous les echantillons ont une couverture medianne appropriee.")


# Fonction pour détecter les valeurs anormales
def detect_outliers(df, column, n_std=nombre_std):
    mean = df[column].mean()
    std = df[column].std()
    lower_bound = mean - n_std * std
    return df.loc[(df[column] < lower_bound) & (
        df[column] < couverture_minimale_normale)]


# Détecter les valeurs anormales pour chaque colonne de couverture
outliers = pd.DataFrame()
for column in [
    'Mosdepth_mqc-generalstats-mosdepth-1_x_pc',
    'Mosdepth_mqc-generalstats-mosdepth-5_x_pc',
    'Mosdepth_mqc-generalstats-mosdepth-10_x_pc',
    'Mosdepth_mqc-generalstats-mosdepth-30_x_pc',
    'Mosdepth_mqc-generalstats-mosdepth-50_x_pc'
]:
    column_outliers = detect_outliers(mosdepth_general_df, column)
    outliers = pd.concat([outliers, column_outliers[['Sample', column]]])

if not (outliers['Mosdepth_mqc-generalstats-mosdepth-1_x_pc'].isna().all() and outliers['Mosdepth_mqc-generalstats-mosdepth-5_x_pc'].isna().all()):
    print("Pour plus de détails sur les échantillons avec des valeurs de " +
          "couverture anormales, veuillez consulter le fichier de sortie " +
          "'erreurs_de_couvertures_anormales.txt'")
    try:
        # Écrire les détails dans le fichier de sortie
        with open('erreurs_de_couvertures_anormales.txt', 'w') as f:
            f.write("Sample\t1x\t5x\t10x\t30x\t50x\n")
            for _, row in outliers.iterrows():
                f.write(f"{row['Sample']}\t{row['Mosdepth_mqc-generalstats-mosdepth-1_x_pc']}\t{row['Mosdepth_mqc-generalstats-mosdepth-5_x_pc']}\t{row['Mosdepth_mqc-generalstats-mosdepth-10_x_pc']}\t{row['Mosdepth_mqc-generalstats-mosdepth-30_x_pc']}\t{row['Mosdepth_mqc-generalstats-mosdepth-50_x_pc']}\n")
    except IOError as e:
        print(f"Erreur lors de l'ouverture ou de l'écriture du fichier: {e}")
    except Exception as e:
        print(f"Une erreur inattendue est survenue: {e}")
else:
    print("Toutes les echantillons ont des couvertures normales.")
