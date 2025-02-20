import argparse
import pandas as pd
import datetime
import os
import re
import time
from datetime import datetime
from check_accession import extract_accession # Function to check if an accession number is real: it uses the entrez functionality of ncbi

def load_dataframe(file_path, columns):
    """Load a DataFrame if the file exists, otherwise create an empty one."""
    if os.path.exists(file_path):
        return pd.read_csv(file_path, sep='\t')
    else:
        return pd.DataFrame(columns=columns)

# Function to check accession number from a FASTA header
def get_accession(header):
    """Extract accession number from a FASTA header."""
    match = re.search(r'\b[A-Z]{1,2}\d{5,6}\b', header)
    return match.group(0) if match else  header.strip()


def update_date(file_path, new_accessions, current_date):
    """Update date_last_updated file with new accessions."""
    df_dates = load_dataframe(file_path, ['accession', 'date_added'])

    # Filter out accessions already present
    new_entries = [
        {'accession': acc, 'date_added': current_date}
        for acc in new_accessions if acc not in df_dates['accession'].values
    ]

    # Append new entries and remove duplicates
    if new_entries:
        df_dates = pd.concat([df_dates, pd.DataFrame(new_entries)]).drop_duplicates(subset='accession', keep='first')

    # Save the updated file
    df_dates.to_csv(file_path, sep='\t', index=False, header=True)
    print(f"Updated {file_path} with {len(new_entries)} new entries.")


def internal_accession(local_df, must_have_local, current_date):
    """Add new local accessions and check if old ids now on Genbank"""

    # Determine the maximum existing internal accession number
    if not local_df.empty:
        existing_ids = local_df['internal_accession'].str.extract(r'(\d+)', expand=False).fillna(0).astype(int)
        max_id = existing_ids.max() if not existing_ids.empty else 0
    else:
        max_id = 0

    # Add a row with the new internal accession
    for accession in must_have_local:
        if not re.match(r'^[A-Z]{1,2}\d{5,6}$', accession):
            max_id += 1
            new_internal_accession = f"X{max_id:07d}"
            local_df = pd.concat([local_df, pd.DataFrame([[new_internal_accession, accession, 'NA', current_date, new_internal_accession]], 
                                                                         columns=['internal_accession', 'sample_name', 'gb_accession', 'date_added', 'seq_accession'])])
            print(f"Created new internal accession: {new_internal_accession} for {accession}")

    # Replace all invalid accessions in gb_accession column to NA
    local_df["gb_accession"] = local_df["gb_accession"].apply(lambda x: 'NA' if not re.match(r'^[A-Z]{1,2}\d{5,6}$', x) else x)
    
    # Check if any of the samples with internal accessions now have a GenBank accession
    for index, row in local_df.iterrows():
        if row['gb_accession'] == 'NA':
            is_real, gb_accession = extract_accession(row['sample_name'])
            if is_real:
                local_df.at[index, 'gb_accession'] = gb_accession
                local_df.at[index, 'seq_accession'] = gb_accession     # Update seq_accession based on gb_accession

                print(f"Updated {row['sample_name']} with GenBank accession: {gb_accession}")    

    # Remove NA sample names and keep only last occurence
    local_df = local_df.dropna(subset='sample_name')
    local_df = local_df.drop_duplicates(subset='sample_name', keep='last')

    return local_df



def process_fasta_headers(fasta_file, existing_accessions):
    """Extract new sequence names from a FASTA file."""
    new_accessions = set()

    with open(fasta_file, 'r') as infile:
        for line in infile:
            if line.startswith('>'):
                header = line[1:].strip()
                accession = get_accession(header)
                if accession and accession not in existing_accessions:
                    new_accessions.add(accession)
        
    return new_accessions

 
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
    parser.add_argument('--ingest_seqs', help = "original sequences from ingest")
    args = parser.parse_args()

    input_sequences = args.in_seq
    output_sequences = args.out_seq
    date_last_updated_file = args.dates
    local_accn_file = args.local_accession
    extended_meta_file = args.add_metadata
    sequences_ingest = args.ingest_seqs

    # Get current date
    current_date = datetime.now().strftime("%Y-%m-%d")

    # Get download date for sequences
    t=time.ctime(os.path.getctime(sequences_ingest))
    seqs_download_date=datetime.strptime(t, "%a %b %d %H:%M:%S %Y").strftime("%Y-%m-%d")

    # Load existing accessions from date_last_updated_file
    existing_accessions = set(
        load_dataframe(date_last_updated_file, ['accession', 'date_added'])['accession'].values
    )
    print("number of accessions:", len(existing_accessions))

    # Load Fasta file headers and get new sequences
    new_accessions = set(process_fasta_headers(input_sequences,existing_accessions))
    print("new accessions:", len(new_accessions))

    # Load local_accession file
    local_df = load_dataframe(local_accn_file, ["internal_accession","sample_name","gb_accession","date_added","seq_accession"])
    
    # Check if any of new_accessions is in local_df; if yes remove from new_accessions
    new_accessions = [accession for accession in new_accessions if accession not in local_df['sample_name'].values and accession not in local_df["seq_accession"].values]
    print("new accessions after removing local ones:", len(new_accessions))

    if len(new_accessions) > 0:
        # Check if real accession number, if not, add to must_have_local set
        must_have_local = set()
        for accession in new_accessions:
            if not re.match(r'^[A-Z]{1,2}\d{5,6}$', accession):
                must_have_local.add(accession)

        print(f"Accessions needing local IDs: {len(must_have_local)}")

        # update local_accn_file with new local ids and check old ids if now on Genbank
        local_df = internal_accession(local_df, must_have_local, seqs_download_date)
        
        local_df.to_csv(local_accn_file, sep='\t', index=False, header=True)

    # Create a dictionary from local_df with sample_name as key and seq_accession as value
    name_to_internal_accn = local_df.set_index('sample_name')['seq_accession'].to_dict()

    # Replace the sample_name in new_accessions with the seq_accession
    updated_new_accessions = [name_to_internal_accn.get(accession, accession) for accession in new_accessions]

    # Update date_last_updated_file with the updated new_accessions list
    update_date(date_last_updated_file, updated_new_accessions, seqs_download_date)

    print("number of accessions after replacing local ones:", len(updated_new_accessions))

    # Replace the headers in fasta with the new accession (in replacement_df)
    with open(input_sequences, 'r') as infile, open(output_sequences, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                sequence_name = line[1:].strip()
                new_header = name_to_internal_accn.get(sequence_name, sequence_name)
                outfile.write(f">{new_header}\n")
            else:
                outfile.write(line)