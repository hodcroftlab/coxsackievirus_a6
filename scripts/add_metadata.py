import pandas as pd
import numpy as np
import re
import sys
import string
# import ipdb
from fuzzywuzzy import process
import argparse
import yaml

#This file adds new metadata to the NCBI Virus metadata
#The files get checked and curated if necessary

# Function to check if an accession number is real: it uses the entrez functionality of ncbi
from check_accession import extract_accession

# Load configuration data from YAML file
with open('config/config.yaml', 'r') as file:
    config = yaml.safe_load(file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='add additional metadata',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input',  metavar=' ', help="input metadata")
    parser.add_argument('-o', '--output', metavar=' ', help="output metadata")
    parser.add_argument('--add', help="input additional data file")
    parser.add_argument('--local', help="input local accession file")
    parser.add_argument('--regions', help="file to specify regions: format = country region")
    parser.add_argument('--id', help="id: strain or accession", choices=["strain","accession"],default="accession")
    parser.add_argument('--update', help="date when sequences were added")
    args = parser.parse_args()

    id_field = args.id
    input_csv_meta = args.input
    output_csv_meta = args.output
    add_data = args.add # if several files, use more than one assignment (add_data_1, add_data_2,...)
    con_reg_table = args.regions
    local_accn = args.local
    last_updated_file = args.update

    # load data
    meta = pd.read_csv(input_csv_meta, keep_default_na=True, sep='\t', index_col=False)
    new_data = pd.read_csv(add_data, keep_default_na=True, sep='\t', index_col=False)
    local_accn_file= pd.read_csv(local_accn, keep_default_na=True, sep='\t', index_col=False)
    last_updated=pd.read_csv(last_updated_file, keep_default_na=True, sep='\t', index_col=False,names=["accession","date_added"])

    # Create a dictionary for quick lookup
    accession_dict = local_accn_file.set_index('sample_name')['seq_accession'].to_dict()

    # Replace missing accessions using the dictionary
    new_data['accession'] = new_data.apply(
        lambda row: accession_dict.get(row['strain'], row['accession']),
        axis=1
    )

    # replace None values with NaN
    new_data['accession'] = new_data['accession'].replace([None], np.nan)

                
    ## Remove duplicates based on id_field
    new_data = new_data.drop_duplicates(subset=id_field)
    meta = meta.drop_duplicates(subset=id_field)

    # If location is missing, replace it with division
    meta['location'] = meta['location'].mask(meta['location'].isna(), meta['division'])

    # step 1: merge both files with to accession number
    new_meta = pd.merge(meta, new_data, on=id_field, how='outer').dropna(subset="accession")

    # add date_added column
    new_meta= pd.merge(new_meta,last_updated, on=id_field,how='left')

    # Creating the new strain column based on the conditions
    new_meta['strain'] = new_meta['strain_y'].mask(new_meta['strain_y'] == new_meta['accession'], new_meta['strain_x'])  # Take strain_x if strain_y == accession
    new_meta['strain'] = new_meta['strain'].mask(new_meta['strain_x'] == new_meta['accession'], new_meta['strain_y'])  # Take strain_y if strain_x == accession
    new_meta['strain'] = new_meta['strain'].mask(new_meta['strain_x'].isna(), new_meta['strain_y'])  # Take strain_y if strain_x is NaN
    new_meta['strain'] = new_meta['strain'].fillna(new_meta['strain_x'])  # Take strain_x if strain_y is NaN

    new_meta.drop(columns="strain_x",inplace=True)

    # Keep only the dates from assign_publications.tsv table - except if they're NA
    new_meta['date'] = new_meta['date_y'].mask(sum([(new_meta['date_y'].isna()),(new_meta['date_y'] == "XXXX-XX-XX")])>=1, new_meta['date_x'])  # If date_y is unknown, keep date_x

    # Region: keep to most detailed one (longest string)
    new_meta['region'] = new_meta['region_x'].mask(new_meta['region_x'].isna(), new_meta['region_y'])
    new_meta['region'] = new_meta['region'].mask(new_meta['region_x'].str.len()<new_meta['region_y'].str.len(), new_meta['region_y'])

    # Country: keep the non-missing ones
    new_meta['country'] = new_meta['country_y'].mask(new_meta['country_y'].isna(), new_meta['country_x'])

    # Country: keep the non-missing ones
    new_meta['place'] = new_meta['location_y'].mask(new_meta['location_y'].isna(), new_meta['location_x'])

    # Clades: keep non-missing clades - subgenogroup
    new_meta['subgenogroup'] = new_meta['subgenogroup'].mask(new_meta['subgenogroup'].isna(), new_meta['clade'])

    # Isolation source: standardize
    # Define a mapping for full terms to their abbreviations and standardized names
    isolation_version = config['metadata']['isolation_source']
    isolation_forms = set(isolation_version.values())

    # Function to map non-standard terms to standard terms
    def standardize_isolation_source(isolation, threshold=75):
        if pd.isna(isolation):
            return np.nan

        # Normalize delimiters
        clean_isolation = (isolation.replace(',', ';')
                                    .replace(' or ', ';')
                                    .replace('/', ';')
                                    .replace('  ', ' ')
                                    .strip('; '))

        # Split on delimiters
        isolations = [iso.strip() for iso in clean_isolation.split(';') if iso.strip()]

        standardized_isolation = []
        for iso in isolations:
            iso_lower = iso.lower()

            # Exact match first
            if iso_lower in isolation_version:
                standardized_isolation.append(isolation_version[iso_lower])
            else:
                # Fuzzy match
                match = process.extractOne(iso_lower, isolation_version.keys(), score_cutoff=threshold)
                if match:
                    standardized_isolation.append(isolation_version[match[0]])
                else:
                    # Keep as title-case or classify as Other
                    standardized_isolation.append(iso.lower())

        # Deduplicate and sort
        standardized_isolation = sorted(set(standardized_isolation))

        return '; '.join([s.lower() for s in standardized_isolation])
    
    # Apply the standardization to both columns
    new_meta['sample_type'] = new_meta['sample_type'].apply(standardize_isolation_source)
    new_meta['isolation'] = new_meta['isolation'].apply(standardize_isolation_source)

    # Combine the two columns, prioritizing the standardized values
    new_meta['combined_isolation_source'] = new_meta['isolation'].mask(new_meta['isolation'].isna(), new_meta['sample_type'])

    # Replace the original isolation_source column
    new_meta['isolation_source'] = new_meta['combined_isolation_source']

    # Drop unnecessary columns
    new_meta = new_meta.drop(['isolation', 'combined_isolation_source'], axis=1)


    # Define a dictionary for old naming formats and spelling mistakes
    corrections = {
        'czech republic': 'czechia',
        'hongkong': 'hong_kong',
        'viet nam': 'vietnam',
        'uk': 'united kingdom',
        'ivory coast': 'côte d\'ivoire',
    }

    # Function to correct country names
    def correct_country_name(country):
        if pd.notna(country):
            country = country.strip().lower()
            corrected_country = corrections.get(country, country).title().replace('_', ' ')
            # print(f"Correcting '{country}' to '{corrected_country}'")  # Debugging statement
            return corrected_country
        return country

    # Apply the corrections to the 'country' column
    new_meta['country'] = new_meta['country'].apply(correct_country_name)

    # Check if the regions file is supplied
    if con_reg_table:
        # Read the regions file and create a dictionary for country-region mappings
        with open(con_reg_table) as f:
            regions = {line.split("\t")[0].strip().lower(): line.split("\t")[1].strip() for line in f.readlines()[1:]}

        # Function to get region from country
        def get_region(coun):
            if pd.notna(coun):
                coun = coun.strip().lower()
                return regions.get(coun, regions.get(coun.replace(' ', '_'), "NA")).replace('_', ' ').title()
            return "NA"

        # Update the 'region' column in the new_meta DataFrame with the new region values
        new_meta['region'] = new_meta['country'].apply(get_region)

    # Debugging statement to check the unique values in the 'country' column after correction
    # print("Unique countries after correction:", new_meta['country'].unique())

    new_meta['has_diagnosis'] =~new_meta['diagnosis'].isna()

    # Define a mapping for full terms to their abbreviations and standardized names
    short_versions = config['metadata']['symptom_list']
    major_versions = config['metadata']['major_symptoms']

    short_forms = set(short_versions.values())
    major_forms = set(major_versions.values())

    def clean_diagnosis(diagnosis, threshold=75):
        if pd.isna(diagnosis):
            return np.nan
        
        # Check if the diagnosis is already a short form
        if diagnosis in short_forms:
            return diagnosis
        
        # Remove punctuation and split multiple diagnoses
        clean_diag = diagnosis.replace(',', ';').replace(' or ', ';').replace('/', ';').replace('  ', ' ').strip()
        diagnoses = [diag.strip() for diag in clean_diag.split(';')]

        # Standardize diagnoses and replace full terms with abbreviations
        standardized_diagnoses = []
        for diag in diagnoses:
            diag_lower = diag.lower()
            if diag_lower in short_versions:
                standardized_diagnoses.append(short_versions[diag_lower])
            else:
                # Use fuzzy matching to handle typos
                match = process.extractOne(diag_lower, short_versions.keys(), score_cutoff=threshold)
                if match:
                    standardized_diagnoses.append(short_versions[match[0]])
                else:
                    # Check if the original diagnosis is in short_forms
                    standardized_diagnoses.append(diag if diag in short_forms else diag.title())
        # Reorder the symptoms so that Fatality is first, then HFMD, then CNS
        standardized_diagnoses = sorted(set(standardized_diagnoses), key=lambda x: (x != 'Fatality', x != 'HFMD', x != 'CNS', x))
        
        # Join the cleaned and standardized diagnoses back into a string
        return '; '.join(sorted(set(standardized_diagnoses)))

    def extract_major_diagnosis(cleaned_diagnosis, threshold=80):
        if pd.isna(cleaned_diagnosis) or cleaned_diagnosis == "":
            return np.nan
        
        # Check if the diagnosis is already a major category
        if cleaned_diagnosis in major_forms:
            return cleaned_diagnosis
        
        # Remove punctuation and split multiple diagnoses
        diagnoses = cleaned_diagnosis.split('; ')

        # Standardize diagnoses and replace full terms with major categories
        major_diagnoses = []
        for diag in diagnoses:
            if diag in major_versions:
                major_diagnoses.append(major_versions[diag])
        
        # Join the cleaned and standardized diagnoses back into a string
        return '; '.join(sorted(set(major_diagnoses)))

    # All the diagnosis
    new_meta['med_diagnosis_all'] = new_meta['diagnosis'].apply(lambda x: clean_diagnosis(x))

    # only the major diagnosis
    new_meta['med_diagnosis_major'] = new_meta['med_diagnosis_all'].apply(lambda x: extract_major_diagnosis(x))

    # Add filter for age and add age ranges
    new_meta['age_yrs'] = pd.to_numeric(new_meta['age_yrs'], errors='coerce')
    new_meta["has_age"] = ~new_meta["age_yrs"].isna()

    # parse gender
    new_meta['gender'] = new_meta['sex']
    new_meta['gender'] = new_meta['gender'].mask(new_meta['gender'].str.contains('female', case=False, na=False), 'F')
    new_meta['gender'] = new_meta['gender'].mask(new_meta['gender'].str.contains('male', case=False, na=False), 'M')

    #Define age bins and labels for years
    bins_years = [-np.inf, 1,4, 6, 11, 18, np.inf]
    labels_years = ['<=1 y/o', '1-3 y/o','4-5 y/o', '6-10 y/o', '11-17 y/o', '18+ y/o']

    # Define age bins and labels for months
    bins_months = [-np.inf, 0.25, 0.5, 1]
    labels_months = ['0-3 m/o', '4-6 m/o', '7-12 m/o']

    # Create age_range column using pd.cut for years
    new_meta['age_range'] = pd.cut(new_meta['age_yrs'], bins=bins_years, labels=labels_years, right=False).astype(str)

    # Handle ages less than 1 year old separately using pd.cut for months
    mask_months = new_meta['age_yrs'] < 1
    new_meta.loc[mask_months, 'age_range'] = pd.cut(new_meta.loc[mask_months, 'age_yrs'], bins=bins_months, labels=labels_months, right=False).astype(str)

    # rename length to NCBI_length_genome
    new_meta.rename(columns={"length": "NCBI_length_genome"}, inplace=True)

    # if ENPEN in origin, set ENPEN to True
    new_meta = new_meta.assign(ENPEN=new_meta['origin'].str.contains('ENPEN', case=False, na=False))

    # if ENPEN=TRUE; Authors in doi should be moved to 'authors'
    # e.g. Private: .... remove Private from authors, keep private in doi
    new_meta['authors'] = new_meta.apply(
        lambda row: row['doi'].replace("Private: ", "") if row['ENPEN'] else row['authors'],
        axis=1
    )
    new_meta['doi'] = new_meta.apply(
        lambda row: "Private" if row['ENPEN'] else row['doi'],
        axis=1
    )

    # write new metadata file to output
    new_meta2= new_meta.loc[:,['accession', 'accession_version', 'strain', 'date', 'region', 'place',
        'country', 'host', 'gender', 'age_yrs','age_range',"has_age", 'has_diagnosis','med_diagnosis_all','med_diagnosis_major',
        'isolation_source', 'NCBI_length_genome',
        'subgenogroup','date_released',
         'abbr_authors', 'authors', 'institution','ENPEN','doi',
        'qc.overallScore', 'qc.overallStatus',
        'alignmentScore', 'alignmentStart', 'alignmentEnd', 'genome_coverage','date_added']]

    new_meta2 = new_meta2.drop_duplicates(subset="accession",keep="first")
    new_meta2.to_csv(output_csv_meta, sep='\t', index=False)