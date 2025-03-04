#!/path/to/python3
import pandas as pd
import numpy as np
import time
import re

# Start overall timing
start_time_total = time.time()

# 1. Read data.
df = pd.read_csv('hidden_file_1.txt', sep='\t')
# Mapping between sequence names and subtype-type
arglist = pd.read_csv('hidden_file_2.csv')

# 2. Filter duplicate annotated ORFs.
# For the file (df), if there are multiple rows with the same values in the 'sample_name_contig'
# and 'contig_orf' columns, then:
#   - Prefer to keep the row where the 'annotation_type' is 'arg'. If multiple rows have 'arg',
#     retain the one with the highest 'pident'; if 'pident' values are equal, retain the row with
#     the highest 'bitscore', and remove the others.
#   - If no row has 'annotation_type' equal to 'arg', then prefer rows with 'annotation_type' equal
#     to 'mge'. If multiple such rows exist, retain the one with the highest 'pident'; if 'pident'
#     values are equal, retain the one with the highest 'bitscore', and remove the others.
#   - If there is no row with 'annotation_type' equal to 'mge', then check whether all rows have
#     'annotation_type' equal to 'mrg': if so, retain the one with the highest 'pident'; if 'pident'
#     values are equal, retain the one with the highest 'bitscore', and remove the others; if not,
#     print the values in the 'sample_name_contig', 'contig_orf', and 'annotation_type' columns.
# Create a new column for ranking the annotation types.
df['type_rank'] = df['annotation_type'].map({'arg': 1, 'mge': 2, 'mrg': 3})
# Sort the data based on 'sample_name_contig', 'contig_orf', 'type_rank', 'pident', and 'bitscore'.
df.sort_values(['sample_name_contig', 'contig_orf', 'type_rank', 'pident', 'bitscore'], 
               ascending=[True, True, True, False, False], inplace=True)
# Remove duplicate rows, keeping only the first occurrence.
df.drop_duplicates(subset=['sample_name_contig', 'contig_orf'], keep='first', inplace=True)
# Remove the temporary 'type_rank' column.
df.drop(columns='type_rank', inplace=True)
# Save intermediate result.
all_annotation_deduplicate = df.copy()
all_annotation_deduplicate.to_csv('hidden_file_3.csv', index=False)

# 3. From all_annotation_deduplicate, select rows where 'annotation_type' is either 'arg' or 'mge'.
result = all_annotation_deduplicate[all_annotation_deduplicate['annotation_type'].isin(['arg', 'mge'])]
# For the rows with annotation_type 'arg', count the number of rows with annotation_type 'mge'
# whose absolute difference in 'orf_location' is within 5000.
# First, get the unique values from the 'sample_name_contig' column.
sample_name_contig_unique = result['sample_name_contig'].unique()
sample_name_contig_list = list(sample_name_contig_unique)
# Loop through each sample in the list.
for sample in sample_name_contig_list:
    # Select rows with annotation_type 'arg' (set A)
    rows_A = result[(result['sample_name_contig'] == sample) & (result['annotation_type'] == 'arg')]
    # Select rows with the same sample and annotation_type 'mge' (set B)
    rows_B = result[(result['sample_name_contig'] == sample) & (result['annotation_type'] == 'mge')]
    for index_1, row_A in rows_A.iterrows():
        # Count rows in B where the absolute difference in 'orf_location' is within 5000.
        num_matches = np.sum(np.abs(rows_B['orf_location'] - row_A['orf_location']) <= 5000)
        # Update the '5kbp_mge_num' column for the row in A with the count.
        result.at[index_1, '5kbp_mge_num'] = num_matches
        print(f"{result.at[index_1, '5kbp_mge_num']} MGEs around {sample}")
# Save intermediate result.
result.to_csv('hidden_file_4.csv', index=False)

# 4. Match arg-id to arg-type-subtype.
# Extract rows with 'annotation_type' equal to 'arg' into a new DataFrame.
all_arg_info = result[result['annotation_type'] == 'arg']
# Add two new columns: 'subtype' and 'type'
all_arg_info['subtype'] = None
all_arg_info['type'] = None
all_arg_info = all_arg_info.reset_index(drop=True)
all_arg_info_check = pd.DataFrame()
# For each row in all_arg_info, search for its 'sequence_id' within the 'gene' column of arglist (case-insensitive).
# If found, copy the corresponding 'subtype' and 'type' values from arglist to all_arg_info.
# Otherwise, append the row to all_arg_info_check.
for index, row in all_arg_info.iterrows():
    match = arglist[arglist['gene'].str.contains(row['sequence_id'], case=False)]
    if not match.empty:
        all_arg_info.loc[index, 'subtype'] = match['subtype'].values[0]
        all_arg_info.loc[index, 'type'] = match['type'].values[0]
    else:
        all_arg_info_check = all_arg_info_check.append(row)
    print(f"{index} has been matched")
# Save intermediate results.
all_arg_info.to_csv('hidden_file_5.csv', index=False)
if not all_arg_info_check.empty:
    all_arg_info_check.to_csv('hidden_file_6.txt', sep='\t', index=False)

# 5. Match antibiotics and consumption based on arg-type-subtype.
# Read data: correspondence between ARGs, antibiotic names, calculated clinical risk, and consumption.
antibiotic_list = pd.read_csv('hidden_file_7.tsv', sep='\t')
arg_clinical_availability = pd.read_csv('hidden_file_8.csv')
antibiotic_comsuption = pd.read_csv('hidden_file_9.csv')
# Add empty columns to all_arg_info.
all_arg_info['subtype_short'] = ''
all_arg_info['Classification'] = ''
all_arg_info['Clinical availability'] = 0.0
# Create an empty DataFrame for consumption check.
comsuption_check = pd.DataFrame(columns=all_arg_info.columns)
# For 'subtype_short', remove everything before and including '__' if present; otherwise, keep the original value.
all_arg_info['subtype_short'] = all_arg_info['subtype'].apply(lambda x: x.split('__')[-1] if '__' in x else x)
# For each row, if the 'subtype_short' is found in the 'ARO Term' column of arg_clinical_availability,
# then copy the 'Classification' and 'Clinical availability' values from arg_clinical_availability.
for index, row in all_arg_info.iterrows():
    if row['subtype_short'] in arg_clinical_availability['ARO Term'].values:
        match_row = arg_clinical_availability[arg_clinical_availability['ARO Term'] == row['subtype_short']].iloc[0]
        all_arg_info.at[index, 'Classification'] = match_row['Classification']
        all_arg_info.at[index, 'Clinical availability'] = match_row['Clinical availability']
        print(f"{index} 1_subtype_short")
# For rows where 'Classification' is still empty, check if 'subtype_short' is contained in antibiotic_list.
# If found, concatenate all matching 'Drug Class' values (separated by semicolons) and assign to 'Classification'.
for i, row in all_arg_info[all_arg_info['Classification'] == ''].iterrows():
    subtype_short = row['subtype_short'].lower()
    matches = antibiotic_list[antibiotic_list['Model Name'].str.contains(re.escape(subtype_short), case=False)]['Drug Class']
    if not matches.empty:
        all_arg_info.loc[i, 'Classification'] = ';'.join(matches.values)
        print(f"{i} 2_Clinical availability")
# For rows where 'Clinical availability' is still 0, further process the 'Classification' field.
for index, row in all_arg_info.iterrows():
    if row['Clinical availability'] == 0:
        classes = row['Classification'].split(';')
        for cls in classes:
            cls = cls.strip().lower()
            for value in antibiotic_comsuption['Classification'].values:
                if cls == value.lower():
                    comsuption_value = antibiotic_comsuption[antibiotic_comsuption['Classification'] == value]['Comsuption'].sum()
                    all_arg_info.at[index, 'Clinical availability'] += comsuption_value
                    print(cls)
    print(f"{index} 3_Clinical availability")
# Save intermediate result.
all_arg_info.to_csv('hidden_file_10.csv', index=False)

# 6. Calculate the average ARG abundance in mining areas based on arg-type-subtype.
# Read data: ARG abundance based on contig names.
arg_abundance = pd.read_csv('hidden_file_11.csv')
no_abundance_arg_cell = pd.read_csv('hidden_file_12.csv', header=None)
# Add an empty column for ARG abundance and create an empty DataFrame for abundance check.
all_arg_info['arg_abundance'] = ''
arg_abundance_check = pd.DataFrame(columns=all_arg_info.columns)

# For each row in all_arg_info, use 'sample_name' and 'subtype' to search in arg_abundance.
# If an exact match is found in the 'subtype' column of arg_abundance, then if the column name in arg_abundance
# (first row) matches exactly, copy the corresponding value to the 'arg_abundance' column.
# If no match is found in the header of arg_abundance, then check in no_abundance_arg_cell:
#   - if found, set the value to 0;
#   - otherwise, append the row to arg_abundance_check.
for index, row in all_arg_info.iterrows():
    sample_name = row['sample_name']
    subtype = row['subtype']
    match_in_arg_abundance = arg_abundance[arg_abundance['subtype'] == sample_name]
    if not match_in_arg_abundance.empty:
        if subtype in arg_abundance.columns:
            all_arg_info.loc[index, 'arg_abundance'] = match_in_arg_abundance[subtype].values[0]
        else:
            if subtype in no_abundance_arg_cell.values:
                all_arg_info.loc[index, 'arg_abundance'] = 0
            else:
                arg_abundance_check = arg_abundance_check.append(row)
    print(f"{sample_name}  {subtype}  has done")
# Save intermediate results.
all_arg_info.to_csv('hidden_file_13.csv', index=False)
arg_abundance_check.to_csv('hidden_file_14.csv', index=False)
all_arg_info.reset_index(inplace=True)

# 7. Calculate the pathogenic proportion of host contigs and compute risk based on arg-type-subtype.
# Extract required columns.
all_risk = all_arg_info[['sample_name', 'subtype', '5kbp_mge_num', 'Clinical availability', 'arg_abundance', 'bac_pathogen']]

# Define a function to compute a special mean for the 'bac_pathogen' column.
def special_mean(values):
    non_no_values = [value for value in values if value != 'NO']
    if len(non_no_values) == 0:
        return 0
    else:
        return len(non_no_values) / len(values)

# Separate samples: those with sample_name starting with 'A' (e.g., control samples) and others (e.g., polluted samples).
all_risk_ck = all_risk[all_risk['sample_name'].str.startswith('A')]
all_risk_pollution = all_risk[~all_risk['sample_name'].str.startswith('A')]

# Group by 'subtype' and aggregate mean values.
all_risk = all_risk.groupby('subtype').agg({
    '5kbp_mge_num': np.mean,
    'Clinical availability': np.mean,
    'arg_abundance': np.mean,
    'bac_pathogen': special_mean
}).reset_index()
# Add a 'risk' column as the product of the aggregated values.
all_risk['risk'] = all_risk['5kbp_mge_num'] * all_risk['Clinical availability'] * all_risk['arg_abundance'] * all_risk['bac_pathogen']

all_risk_ck = all_risk_ck.groupby('subtype').agg({
    '5kbp_mge_num': np.mean,
    'Clinical availability': np.mean,
    'arg_abundance': np.mean,
    'bac_pathogen': special_mean
}).reset_index()
all_risk_ck['risk'] = all_risk_ck['5kbp_mge_num'] * all_risk_ck['Clinical availability'] * all_risk_ck['arg_abundance'] * all_risk_ck['bac_pathogen']

all_risk_pollution = all_risk_pollution.groupby('subtype').agg({
    '5kbp_mge_num': np.mean,
    'Clinical availability': np.mean,
    'arg_abundance': np.mean,
    'bac_pathogen': special_mean
}).reset_index()
all_risk_pollution['risk'] = all_risk_pollution['5kbp_mge_num'] * all_risk_pollution['Clinical availability'] * all_risk_pollution['arg_abundance'] * all_risk_pollution['bac_pathogen']

# Save the risk calculation results.
all_risk.to_csv('hidden_file_15.csv', index=False)
all_risk_ck.to_csv('hidden_file_16.csv', index=False)
all_risk_pollution.to_csv('hidden_file_17.csv', index=False)

# End overall timing and print the total running time.
end_time_total = time.time()
print(f"Total Running Time : {end_time_total - start_time_total} seconds")
