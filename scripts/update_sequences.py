import argparse
import pandas as pd
import datetime
import os
import re
import time
from datetime import datetime

# Function to check if an accession number is real: it uses the entrez functionality of ncbi
from check_accession import extract_accession

# Function to extract the accession number from a FASTA header
def get_accession(header):
    match = re.search(r'\b[A-Z]{1,2}\d{5,6}\b', header)
    return match.group(0) if match else header.strip()
  
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='update_sequences',
        description="""
        Update sequences and accession files
        """,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('--in_seq', required=True, help="Input sequences file")
    parser.add_argument('--out_seq', required=True, help="Output sequences file")
    parser.add_argument('--dates', help="Date last updated file")
    parser.add_argument('--local_accession', help="Local accession file")
    parser.add_argument('--add_metadata', help="Additional extended metadata file")
    parser.add_argument('--metadata', help="metadata file")
    args = parser.parse_args()

    input_sequences = args.in_seq
    output_sequences = args.out_seq
    date_last_updated_file = args.dates
    local_accn_file = args.local_accession
    extended_meta_file = args.add_metadata


    # Get current date
    current_date = datetime.now().strftime("%Y-%m-%d")

    # Get download date for sequences
    t=time.ctime(os.path.getctime("data/sequences.fasta"))
    seqs_download_date=datetime.strptime(t, "%a %b %d %H:%M:%S %Y").strftime("%d-%m-%Y")

    # read table with extended (published) metadata
    extended_meta = pd.read_csv(extended_meta_file, sep='\t')
    extended_meta = extended_meta.loc[~extended_meta.accession.isna(), ['accession', 'strain']]

    # Read existing date_last_updated.txt file
    if os.path.exists(date_last_updated_file):
        df_dates = pd.read_csv(date_last_updated_file, sep='\t')
    else:
        df_dates = pd.DataFrame(columns=['accession', 'date_added'])

    # Read existing local_accn.txt file
    if os.path.exists(local_accn_file):
        df_local_accn = pd.read_csv(local_accn_file, sep='\t',)
    else:
        df_local_accn = pd.DataFrame(columns=['internal_accession', 'sample_name', 'gb_accession', 'date_added'])



    # Add entries for accessions or strain names present in the extended_meta
    for index, row in extended_meta.iterrows():
        accession, strain = row['accession'], row['strain']
        
        # Check if the accession or strain is already present
        if strain in df_local_accn['sample_name'].values and accession not in df_local_accn['gb_accession'].values:
            old_internal_accession = df_local_accn.loc[df_local_accn['sample_name'] == strain, 'internal_accession'].values[0] # get interal_accession from strain
            df_local_accn = pd.concat([df_local_accn, pd.DataFrame([[old_internal_accession, strain, accession, current_date]], columns=['internal_accession', 'sample_name', 'gb_accession', 'date_added'])])
            df_dates = pd.concat([df_dates, pd.DataFrame([[accession, current_date]], columns=['accession', 'date_added'])])

    # remove duplicated accessions or duplicated strain names
    df_local_accn = df_local_accn.drop_duplicates(subset=['internal_accession', 'sample_name'], keep='last')
    df_dates = df_dates.drop_duplicates(subset='accession', keep='last')


    # Step 1: Update date_last_updated.txt
    # Extract sequence names from the input FASTA file
    sequence_names = set()
    not_accession = set()
    local_accn_dic = {}
    with open(input_sequences, 'r') as seq_file:
        for line in seq_file:
            if line.startswith('>'):
                sequence_name = line[1:].strip()
                # Check if the sequence is new
                if sequence_name in df_dates['accession'].values or sequence_name in df_local_accn['sample_name'].values:
                    continue  # Skip this sequence if it's already in either the date or local accession files

                # Proceed with processing for new sequences
                if '|' in sequence_name:
                    local_accn_dic[sequence_name] = get_accession(sequence_name)
                    sequence_name = get_accession(sequence_name)
                    
                # Directly consider it a real accession if it matches the pattern
                elif re.match(r'^[A-Z]{2}\d{6}$', sequence_name): 
                    pass

                else:
                    is_real, accession = extract_accession(sequence_name)
                    
                    if is_real and not re.match(r'^[A-Z]{2}\d{6}$', accession):
                        is_real = False
                    if not is_real:
                        not_accession.add(sequence_name)
                        sequence_name = None
                    elif accession:
                        local_accn_dic[sequence_name] = accession
                        sequence_name = accession
                        
                if sequence_name:
                    sequence_names.add(sequence_name)


    # Update df_dates DataFrame
    new_entries = pd.DataFrame({'accession': list(sequence_names), 'date_added': current_date})
    df_dates = pd.concat([df_dates, new_entries]).drop_duplicates(subset='accession', keep='last')
    df_dates=df_dates.reset_index(drop=True)


    # Save updated date_last_updated.txt file; further down now
    # df_dates.to_csv(date_last_updated_file, sep='\t', index=False, header=True)

    # Step 2: Update local_accn.txt
    # Determine the maximum existing internal accession number
    if not df_local_accn.empty:
        # Replace NaN values with a placeholder before converting to integers
        existing_ids = df_local_accn['internal_accession'].str.extract(r'(\d+)', expand=False).fillna(0).astype(int)
        max_id = existing_ids.max() if not existing_ids.empty else 0
    else:
        max_id = 0


    # Create DataFrame for new entries
    new_entries = []
    if not not_accession:
        for name in not_accession:
            if name not in df_local_accn['sample_name'].values:
                # Assuming typical GenBank accession
                if len(name) == 8 and name[:2].isalpha() and name[2:].isdigit():  
                    gb_accession = name
                else:
                    max_id += 1
                    new_internal_accession = f"X{max_id:07d}"
                    gb_accession = 'NA'
                    new_entries.append([new_internal_accession, name, gb_accession, current_date])

    for name in local_accn_dic.keys():
        if re.match(r'^[A-Z]{2}\d{6}$', name):
            local_accn_dic[name] = None


    df = pd.DataFrame(local_accn_dic.items(), columns=['sample_name', 'gb_accession'])
    df=df.loc[df['gb_accession'].apply(lambda x: len(x) == 8 and x[:2].isalpha() and x[2:].isdigit()), :]
    df['date_added'] = current_date

    df_new_local_accn = pd.DataFrame(new_entries, columns=['internal_accession', 'sample_name', 'gb_accession', 'date_added'])

    # Concatenate dataframes
    df_local_accn = pd.concat([df_local_accn, df_new_local_accn, df])

    # Drop rows where all columns except 'sample_name' are NA
    df_local_accn = df_local_accn.dropna(subset='sample_name')

    # Then drop duplicates based on 'sample_name', keeping the last occurrence
    df_local_accn = df_local_accn.drop_duplicates(subset='sample_name', keep='last')

    # Define a lambda function to select 'gb_accession' if it is not NA, otherwise 'internal_accession'
    df_local_accn['seq_accession'] = df_local_accn.apply(
        lambda row: row['gb_accession'] if pd.notna(row['gb_accession']) else row['internal_accession'],
        axis=1)


    # Save updated local_accn.txt file
    df_local_accn.to_csv(local_accn_file, sep='\t', index=False, header=True)


    # Step 3: Replace headers in the sequences and save to output file
    # Create a mapping from sample_name to internal_accession

    # open date_last_updated_file and add new names
    df=df_local_accn.loc[:, ["seq_accession","date_added"]].reset_index(drop=True)
    df=df.rename(columns={"seq_accession":"accession"})
    df_dates=df_dates.reset_index(drop=True)

    # merge them
    df_dates_c = pd.concat([df_dates, df]).drop_duplicates(subset='accession', keep='last')

    # Read metadata
    meta = pd.read_csv(args.metadata, sep='\t')

    # Drop duplicates and convert to DataFrame
    accn_list = meta[['accession']].drop_duplicates()
    accn_list['date_added'] = seqs_download_date

    # Concatenate with the existing DataFrame
    df_dates_c = pd.concat([df_dates_c, accn_list]).sort_values(by="date_added", ascending=True).drop_duplicates(subset='accession', keep='last')
    df_dates_c.to_csv(date_last_updated_file, sep='\t', index=None)

    # Create the dictionary
    name_to_internal_accn = df_local_accn.set_index('sample_name')['seq_accession'].to_dict()

    # Replace the headers in fasta
    with open(input_sequences, 'r') as infile, open(output_sequences, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                sequence_name = line[1:].strip()
                new_header = name_to_internal_accn.get(sequence_name, sequence_name)
                outfile.write(f">{new_header}\n")
            else:
                outfile.write(line)