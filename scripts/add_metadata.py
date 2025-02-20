import pandas as pd
import numpy as np
import re
import sys
import string
import ipdb
from fuzzywuzzy import process
import argparse

#This file adds new metadata to the NCBI Virus metadata
#The files get checked and curated if necessary

# Function to check if an accession number is real: it uses the entrez functionality of ncbi
from check_accession import extract_accession

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='add additional metadata',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-i', '--input',  metavar=' ', help="input metadata")
    parser.add_argument('-o', '--output', metavar=' ', help="output metadata")
    parser.add_argument('--add', help="input additional data file")
    parser.add_argument('--local', help="input local accession file")
    parser.add_argument('--regions', help="file to specify regions: format = country region")
    parser.add_argument('--id', help="id: strain or accession", choices=["strain","accession"],default="accession")
    parser.add_argument('--rename', help="output rule update_strain_names")
    parser.add_argument('--update', help="date when sequences were added")
    args = parser.parse_args()

    id_field = args.id
    input_csv_meta = args.input
    output_csv_meta = args.output
    add_data = args.add # if several files, use more than one assignment (add_data_1, add_data_2,...)
    con_reg_table = args.regions
    local_accn = args.local
    renamed_strains = args.rename
    last_updated_file = args.update
    

    # load data
    meta = pd.read_csv(input_csv_meta, keep_default_na=True, sep='\t', index_col=False)
    new_data = pd.read_csv(add_data, keep_default_na=True, sep='\t', index_col=False)
    local_accn_file= pd.read_csv(local_accn, keep_default_na=True, sep='\t', index_col=False)
    renamed_strains_df = pd.read_csv(renamed_strains, keep_default_na=True, sep='\t', index_col=False, on_bad_lines='skip')
    last_updated=pd.read_csv(last_updated_file, keep_default_na=True, sep='\t', index_col=False,names=["accession","date_added"])

    # Identify records that need updating
    needs_update = meta[~meta['accession'].isin(renamed_strains_df['accession'])]
    needs_update.to_csv("data/no_strain_correction.tsv", sep='\t', index=False)

    # Create a lookup dictionary for strain updates
    lookup_strain = renamed_strains_df.set_index('accession')['strain'].to_dict()

    # Update strains in metadata according to lookup
    meta['strain'] = meta.apply(
        lambda row: lookup_strain.get(row['accession'], row['strain']),
        axis=1
    )

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
    # Function to map non-standard terms to standard terms
    def standardize_isolation_source(value):
        # Add mappings for non-standard terms to standard terms
        mapping = {
            'stool': 'feces',
            'Stool': 'feces',
            'Faeces': 'feces',
            'CSF': 'cerebrospinal fluid',
            'Cerebrospinal Fluid': 'cerebrospinal fluid',
            'Cerebrospinal Fluid ': 'cerebrospinal fluid',
            'Throat Swab': 'oronasopharynx',
            'Throat swab': 'oronasopharynx',
            'Oral Swab': 'oronasopharynx',
            'Mouth Swab': 'oronasopharynx',
            'Nasopharyngeal Aspirate': 'oronasopharynx',
            'Brain': 'brain',
            'Cerebellum': 'brain',
            'Left brain': 'brain',
            'Brainstem': 'brain',
            'Vesicle Swab': 'swab',
            'Ulcer Swab': 'lesion',
            'Rectal swab': 'swab',
            'Rectal Swab': 'swab',
            'Mouth': 'oronasopharynx',
            'Throat': 'oronasopharynx',
            'Tonsil': 'oronasopharynx',
            'Vesicle': 'lesion',
            'Spinal cord ': 'spinal cord'
        }
        
        val=mapping.get(value, value)
        if pd.isna(val):
            return val
        val=val.title()

        # Return the mapped value if it exists, otherwise return the original value
        return val

    # Apply the standardization to both columns
    new_meta['sample_type'] = new_meta['sample_type'].apply(standardize_isolation_source)
    new_meta['isolation'] = new_meta['isolation'].apply(standardize_isolation_source)

    # Combine the two columns, prioritizing the standardized values
    new_meta['combined_isolation_source'] = new_meta['isolation'].mask(new_meta['isolation'].isna(), new_meta['sample_type'])

    # Replace the original isolation_source column
    new_meta['isolation_source'] = new_meta['combined_isolation_source']

    # Drop unnecessary columns
    new_meta = new_meta.drop(['isolation', 'combined_isolation_source'], axis=1)


    # Check if the regions file is supplied
    if con_reg_table:
        # Dictionary to store country-region mappings
        regions = {}
        
        # Read the regions file
        with open(con_reg_table) as f:
            regs = f.readlines()
        
        # Remove the first line (header)
        regs = regs[1:]
        
        # Populate the regions dictionary
        for x in regs:
            x = x.strip()
            pair = x.split("\t")
            if len(pair) > 1:
                regions[pair[0].strip().lower()] = pair[1].strip()
        
        # List to store the new region values
        newregion = []
        
        # Iterate over each country in the 'country' column of new_meta
        for coun in new_meta['country']:
            reg = "NA"  # Default value if no region is found
            if pd.notna(coun):
                coun = coun.strip().lower()
                coun_with_underscores = coun.replace(' ', '_')
                
                if coun not in regions:
                    if coun_with_underscores not in regions:
                        print(f"No region found for {coun}! Setting to NA")
                    else:
                        reg = regions[coun_with_underscores]
                else:
                    reg = regions[coun]
            
            reg = reg.replace('_', ' ')
            reg = reg.title()
            newregion.append(reg)
        
        # Update the 'region' column in the new_meta DataFrame with the new region values
        new_meta['region'] = newregion


    ## Diagnosis
    # ipdb.set_trace()

    new_meta['has_diagnosis'] =~new_meta['diagnosis'].isna()

    # Define a mapping for full terms to their abbreviations and standardized names
    short_versions = {
    'Poliomyelitis-like disease': 'AFP', 'acute flaccid paralysis': 'AFP',
    'Paralysis': 'AFP', 'Poliomyelitis-like paralysis': 'AFP',
    "AFP patient": 'AFP',
    'patient with AFP':'AFP',
    'Upper Respiratory Tract Infection': 'URTI',
    'asthmatic bronchitis': "Respiratory Symptoms",
    'Pharyngitis':'Pharyngitis',
    'Guillain-Barré syndrome': 'GBS', 'Guillain-Barrésyndrome': 'GBS',
    'Febrile illness': 'Fever', 'Pyrexia': 'Fever',
    'V&D': 'Vomiting; Diarrhea', 'Diarrhoea': 'Diarrhea', '(diarrhea)': 'Diarrhea',
    'Severe': "Severe HFMD", 'Light': "Light HFMD", 'mild': "Mild HFMD",
    'Atypical hand foot and mouth disease': "Atypical HFMD",
    'hand-foot and mouth disease': 'HFMD',
    'hand foot and mouth disease': 'HFMD',
    'hand, foot and mouth disease': 'HFMD',
    'herpangina': 'Herpangina',
    'Oral ulcer': 'Ulcer', 
    'vesicular': 'vesicular papules', 'Rash with vesicles': 'vesicular papules',
    '(death)': 'Fatality', 'fatal': 'Fatality', 'death': 'Fatality',
    "Myoclonic jerk": 'Myoclonus',
    'CNS symptoms': 'CNS', 'CNS involvement': 'CNS', 'Cns Disorder': 'CNS',
    'Neurological': 'CNS', 'severe neurological involvement': 'CNS',
    'polio': 'AFP', 'Difficulty Walking': 'AFP',
    'Neck Stiffness; Vomiting': 'Meningitis',
    'aseptic meningitis': 'Aseptic Meningitis', 
    'Meningitis': 'Meningitis',
    'Meningoencephalitis': 'Encephalitis', 'Encephalytis': 'Encephalitis', 'encephalitis':'Encephalitis',
    'Skin Rash': 'Rash', 'exanthema': 'Rash'
    }

    short_forms = set(short_versions.values())

    major_versions = {
        'Encephalitis': 'CNS',
        'Aseptic Meningitis': 'CNS',
        'Meningitis': 'CNS',
        'Opsomyoclonus Syndrome': 'CNS',
        'Acute Cardiogenic Shock': 'Fatality',
        'AFP':'AFP',
        'HFMD':'HFMD',
        'Atypical HFMD':"Atypical HFMD",
        'Fatality':'Fatality',
        'CNS':'CNS',
    }

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

    # add VP1 length as continuous variable
    ## read in blast output length
    ##TODO: change to counting the length of the sequence after alignment and without gaps - done after tree building
    # blast_results=pd.read_csv("vp1/results/blast_vp1_length.csv",names=["accession","l_vp1"])
    # new_meta2=pd.merge(new_meta,blast_results, on="accession", how='left')

    # bins_length = [-np.inf, 599, 699,799,899, np.inf]
    # labels_length = ['<600nt', '600-700nt', '700-800nt','800-900nt','>900nt']

    # # Create length_range column using pd.cut for length
    # new_meta2['length_VP1'] = pd.cut(new_meta2['l_vp1'], bins=bins_length, labels=labels_length, right=False).astype(str)

    # ipdb.set_trace()
    # write new metadata file to output
    new_meta2= new_meta.loc[:,['accession', 'accession_version', 'strain', 'date', 'region', 'place',
        'country', 'host', 'gender', 'age_yrs','age_range',"has_age", 'has_diagnosis','med_diagnosis_all','med_diagnosis_major',
        'isolation_source', 'length',#'length_VP1',
        'subgenogroup','lineage','date_released',
        'abbr_authors', 'authors', 'institution','doi',
        'qc.overallScore', 'qc.overallStatus',
        'alignmentScore', 'alignmentStart', 'alignmentEnd', 'genome_coverage','date_added']]


    new_meta2.to_csv(output_csv_meta, sep='\t', index=False)