#!/home/ubuntu/anaconda3/envs/args_oap/bin/python3
import pandas as pd
import numpy as np
import os
import concurrent.futures
import time
from Bio import SeqIO
import sys
import threading

# Path to blast-annotation result
folder_path = '/path/to/data/folder'
# Get .txt filenames
all_file_names = [f for f in os.listdir(folder_path) if f.endswith('.txt')]
# Deduplicate sample names
all_sample_name = set([f.rsplit('_', 1) for f in all_file_names])
# Create an empty list to store all data tables
all_data_list = []

# Start total time
start_time_total = time.time()

# Define a function to process each sample
for sample in all_sample_name:
    # Start time for sample
    start_time = time.time()
    # Files to search for
    search_files = [sample + '_arg.txt', sample + '_mrg.txt', sample + '_mge.txt']
    # Filter files containing the sample name
    files = {file1 for file1 in all_file_names if file1 in search_files}
    # Create a new DataFrame to process data
    all_data = pd.DataFrame()

    for file in files:
        data = pd.read_csv(os.path.join(folder_path, file), sep='\t', header=None)
        # Remove .txt from filename
        filename = file[:-4]
        # Split file name
        sample_name, annotation_type = filename.rsplit('_', 1)
        # Add new columns
        data.insert(0, 'sample_name', sample_name)
        data.insert(1, 'annotation_type', annotation_type)
        data.insert(2, 'contig', data.str.rsplit('_', n=1).str)
        # Append data to DataFrame
        all_data = pd.concat([all_data, data], ignore_index=True)
    
    # Add column names
    all_data.columns = ['sample_name', 'annotation_type', 'contig', 'contig_orf', 'sequence_id', 'pident', 'alignment_length', 'mismatch', 'gap_openings', 'query_start', 'query_end', 'subject_start', 'subject_end', 'evalue', 'bitscore']
    
    # Add new empty columns
    all_data['bac_name'] = ''
    all_data['bac_id'] = ''
    all_data['bac_pathogen'] = ''
    all_data['bac_level'] = ''
    all_data['orf_start'] = ''
    all_data['orf_end'] = ''
    all_data['orf_head'] = ''
    all_data['orf_location'] = ''
    all_data['5kbp_mge_num'] = 0
    
    # Read pathogen list file
    patonlist = pd.read_csv('/path/to/pathogen_list.txt', sep='\t', header=None)
    # Read kraken output and report files
    sample_r = sample.replace('-', '_')
    out_file = pd.read_csv(f'/path/to/kraken/output/{sample_r}.out', sep='\t', header=None)
    report_file = pd.read_csv(f'/path/to/kraken/report/{sample_r}.report', sep='\t', header=None)
    
    # Process each row in the DataFrame
    for index, row in all_data.iterrows():
        if row['sample_name'] == sample:
            # Match bacterial ID
            matching_rows1 = out_file[out_file == row['contig']]
            if not matching_rows1.empty:
                all_data.at[index, 'bac_id'] = matching_rows1.iloc[0, 2]
            index_value_bac_id = all_data.at[index, 'bac_id']
            if isinstance(index_value_bac_id, pd.Series):
                index_value_bac_id = index_value_bac_id.values
            # Match bacterial name
            matching_rows2 = report_file[report_file.iloc[:, 4] == index_value_bac_id]
            if not matching_rows2.empty:
                all_data.at[index, 'bac_name'] = matching_rows2.iloc[0, 5].strip()
                all_data.at[index, 'bac_level'] = matching_rows2.iloc[0, 3].strip()
            # Check if it's a pathogen
            index_value_bac_name = all_data.at[index, 'bac_name']
            if isinstance(index_value_bac_name, pd.Series):
                index_value_bac_name = index_value_bac_name.values
            matching_rows3 = patonlist[patonlist == index_value_bac_name]
            if not matching_rows3.empty:
                all_data.at[index, 'bac_pathogen'] = matching_rows3.values
            else:
                all_data.at[index, 'bac_pathogen'] = 'NO'
    
    # Process ORF information
    orf_records = SeqIO.to_dict(SeqIO.parse(f'/path/to/orf/fasta/{sample}_pro.fa', 'fasta'))
    sample_rows = all_data[all_data['sample_name'] == sample]
    for index, row in sample_rows.iterrows():
        orf_record = orf_records[row['contig_orf']]
        annotations = orf_record.description.split('#')
        all_data.at[index, 'orf_start'] = annotations.strip()
        all_data.at[index, 'orf_end'] = annotations.strip()
        all_data.at[index, 'orf_head'] = annotations.strip()
        all_data.at[index, 'orf_location'] = (int(annotations) + int(annotations)) / 2
    
    # Create a new column
    all_data['sample_name_contig'] = all_data['sample_name'] + "_" + all_data['contig']
    # Print running time for the sample
    end_time = time.time()
    print(f"{sample} Running time : {end_time - start_time} seconds")
    # Append the DataFrame to the list
    all_data_list.append(all_data)

# Print total running time
end_time_total = time.time()
print(f"Total Running Time : {end_time_total - start_time_total} seconds")

# Combine all DataFrames
all_data_combined = pd.concat(all_data_list, ignore_index=True)
# Save the result
all_data_combined.to_csv('all_annotation_pathogen_mge_num.txt', sep='\t', index=False)

# Check the result
result = pd.DataFrame()
groups = all_data_combined.groupby('sample_name_contig')
for name, group in groups:
    if not (group['bac_id'].nunique() == 1 and group['bac_name'].nunique() == 1 and group['bac_pathogen'].nunique() == 1):
        result = pd.concat([result, group], ignore_index=True)
result.to_csv('check_result.txt', sep='\t', index=False)
