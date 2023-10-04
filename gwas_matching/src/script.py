#!/usr/bin/env python

import pandas as pd
from glob import glob
import re

gwas_dict = {}
paths = set()

# Read file and create dict for O(1) lookup
gwas_raw = pd.read_csv('vars_overlapping_testable_gene_windows.txt', sep='\t')

for index, row in gwas_raw.iterrows():
    if row['Gene'] in gwas_dict:
        gwas_dict[row['Gene']].append([row['Chr'], int(row['Variant Position'])])
    else:
        gwas_dict[row['Gene']] = [[row['Chr'], int(row['Variant Position'])]]

def get_highest_peak_number(dir_path):
    """Returns the highest number of peaks based on filenames in the directory."""
    peak_numbers = []
    for file_name in glob(f"{dir_path}/*peaks_*.csv"):
        match = re.search(r'_(\d+)peaks_', file_name.split('/')[-1])
        if match:
            peak_numbers.append(int(match.group(1)))
    return max(peak_numbers) if peak_numbers else None

def get_file_path(dir_path, num_peaks):
    """Returns the file path based on number of peaks."""
    for file_name in glob(f"{dir_path}/*{num_peaks}peaks_*.csv"):
        match = re.search(r'_(\d+)peaks_', file_name.split('/')[-1])
        if match and int(match.group(1)) == num_peaks:
            return file_name
    return None

result = []

# Per Each Gene
for dir in glob("Results/*"):
    gene = dir.split("/")[1]
    try:
        selected_peaks = pd.read_csv(dir + "/selected_peaks.csv")
    except:
        print("Skipped", dir, "Because it doesn't contain selected peaks")
        continue

    mapping_dict = {
        'RFR_rf_ranker': 'rf_ranker',
        'RFR_perm_ranker': 'rf_permranker',
        'RFR_dropcol_ranker': 'rf_dropcolranker',
        'LinReg_perm_ranker': 'linreg_permranker',
        'LinReg_dropcol_ranker': 'linreg_dropcolranker'
    }

    selected_peaks['Ranker_Method'] = selected_peaks['Result'].str.extract(r'_(\w+_\w+_ranker)_')

    # Apply the mapping to the Ranker_Method column
    selected_peaks['Mapped_Ranker_Method'] = selected_peaks['Ranker_Method'].map(mapping_dict)

    # Per Each Model
    for index, model in selected_peaks.iterrows():
        model_specific_peaks_dir = dir + "/" + model['Mapped_Ranker_Method']
        peaks = model['nPeaks']

        # print(model_specific_peaks_dir)

        file_path_peaks = get_file_path(model_specific_peaks_dir, peaks)
        
        highest_peak_num = get_highest_peak_number(model_specific_peaks_dir)
        file_path_highest_peaks = get_file_path(model_specific_peaks_dir, highest_peak_num)

        # print("Model-specific peaks directory:", model_specific_peaks_dir)
        # print("File for number of peaks:", file_path_peaks)
        # print("File for highest number of peaks:", file_path_highest_peaks)
        # print("----------------------------------------------------")

        try:
            gwas_data = gwas_dict[gene]
        except:
            # print("Didn't process", gene, "because it isn't present in reference tsv")
            continue

        for snp in gwas_data:
            selected_df = pd.read_csv(file_path_peaks)
            total_df = pd.read_csv(file_path_highest_peaks)
            gwas_confirmed_selected = []
            total_confirmed = []

            for index, row in selected_df.iterrows():
                chromosome, var_range = row['Peak'].split(":")
                start = int(var_range.split("-")[0])
                end = int(var_range.split("-")[1])
                
                if chromosome != snp[0]:
                    continue

                if snp[1] >= start and snp[1] <= end:
                    gwas_confirmed_selected.append(snp[1])
                    result.append({
                        "GENE": gene,
                        "SNP_COORD": chromosome + ":" + str(snp[1]),
                        "ATAC_COORD": chromosome + ":" + str(start) + "-" + str(end),
                        "MODEL": model['Mapped_Ranker_Method'],
                        "SELECTED": "Yes"
                    })

            for index, row in total_df.iterrows():
                chromosome, var_range = row['Peak'].split(":")
                start = int(var_range.split("-")[0])
                end = int(var_range.split("-")[1])
                
                if chromosome != snp[0]:
                    continue

                if snp[1] >= start and snp[1] <= end:
                    if snp[1] not in gwas_confirmed_selected:
                        total_confirmed.append(snp[1])
                        result.append({
                            "GENE": gene,
                            "SNP_COORD": chromosome + ":" + str(snp[1]),
                            "ATAC_COORD": chromosome + ":" + str(start) + "-" + str(end),
                            "MODEL": model['Mapped_Ranker_Method'],
                            "SELECTED": "No"
                        })

result_df = df = pd.DataFrame(result)

print(result_df)