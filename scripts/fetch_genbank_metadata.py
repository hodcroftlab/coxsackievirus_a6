#!/usr/bin/env python3
"""
Fetch GenBank metadata for viral accessions with config-driven symptom and isolation mappings.

This script can be used standalone or integrated into Snakemake workflows.

python fetch_genbank_metadata.py \
            --virus "<virus>" \
            --accession_file "<accession list file>"  \
            --output "<output-tsv file>"  \
            --genbank "<summary genbank metafile>" \
            --config config/config.yaml
"""

import argparse
import pandas as pd
import yaml
import re
import os
import sys
from Bio import Entrez, SeqIO
from dotenv import load_dotenv, find_dotenv
from typing import Dict, List, Optional, Tuple
# import ipdb
import word2number as w2n

class MetadataFetcher:
    """Handles fetching and parsing GenBank metadata for viral sequences."""
    
    def __init__(self, virus_name: str, config: Dict, email: str,location_list):
        """
        Initialize the metadata fetcher.
        
        Args:
            virus_name: Name of the virus (e.g., "Coxsackievirus A10")
            config: Configuration dictionary with symptom_list and isolation_source mappings
            email: Email for NCBI Entrez
        """
        self.virus_name = virus_name
        self.symptom_list = self._normalize_dict_keys(config['metadata']['symptom_list'])
        self.isolation_sources = self._normalize_dict_keys(config['metadata']['isolation_source'])
        self.location_list = location_list
        # Set up Entrez
        Entrez.email = email
        
    @staticmethod
    def _normalize_dict_keys(d: Dict) -> Dict:
        """Convert all dictionary keys to lowercase for case-insensitive matching."""
        return {k.lower(): v for k, v in d.items()}
    
    def parse_host(self, field: str) -> Tuple[Optional[str], Optional[str], Optional[str], 
                                                Optional[float], Optional[float], Optional[str]]:
        """
        Parse host/isolation field to extract metadata.
        
        Args:
            field: The field string to parse (from host or isolation_source)
            
        Returns:
            Tuple of (host, isolation, sex, age_yrs, age_mo, diagnosis)
        """
        try:
            sex, age_yrs, age_mo, diagnosis = None, None, None, None
            host, isolation = None, None

            # if not field:
            #     return host, isolation, sex, age_yrs, age_mo, diagnosis
            
            field_lower = field.lower()

            # Regular expressions for parsing
            host_regex = r"(homo sapiens|human|other hosts|child)"
            sex_regex = r"\b(male|female|m|f|man|woman|boy|girl|h)\b"
            age_regex = r"(\d+\.?\d*|one|two|three|four|five|six|seven|eight|nine|ten|eleven|twelve)\s*(years?|months?|days?|y|m|d)"
            symptom_pattern = "|".join(map(re.escape, self.symptom_list.keys()))
            diagnosis_regex = r"\b(?:" + symptom_pattern + r")\b"
            isolation_regex = r"\b(feces|throat|plasma|csf|cerebro spinal fluid|mouth swab|vesicle|stool|blood|serum)\b"
            
            # Match host
            host_match = re.search(host_regex, field_lower)
            if host_match:
                host = host_match.group(0).capitalize()
            
            # Match sex
            sex_match = re.search(sex_regex, field_lower)
            if sex_match:
                sex_value = sex_match.group(0).lower()
                if sex_value in ["f", "girl", "woman", "female"]:
                    sex = "F"
                elif sex_value in ["m", "h", "boy", "man", "male"]:
                    sex = "M"
            
            # Match age
            age_matches = re.findall(age_regex, field_lower)
            
            if age_matches:
                for match in age_matches:
                    value, unit = match
                    try:
                        value = float(value)
                    except ValueError:
                        value = w2n.word_to_num(value)
                    if "day" in unit or "d" in unit:
                        age_yrs = round(value / 365, 3)
                    elif "year" in unit or "y" in unit:
                        age_yrs = value
                    elif "month" in unit or "m" in unit:
                        age_mo = value
                    
                if len(age_matches) >1:
                    age_yrs = age_yrs + round(age_mo / 12, 3)
                    age_mo = age_yrs *12     

                elif age_yrs is None and age_mo is not None:
                    age_yrs = round(age_mo / 12, 3)
                elif age_yrs is not None:
                    age_mo = round(age_yrs * 12, 3)

            # Match diagnosis from symptom list
            diagnosis_matches = re.findall(diagnosis_regex, field_lower, re.IGNORECASE)
            diagnosis_list = []
            
            if diagnosis_matches:
                for match in diagnosis_matches:
                    standardized = self.symptom_list.get(match.lower())
                    
                    if standardized:
                        diagnosis_list.append(standardized)
            
            # Remove duplicates and join
            if diagnosis_list:
                diagnosis_list = list(set(diagnosis_list))
                diagnosis = "; ".join(diagnosis_list)
            
            # Match isolation source
            isolation = field_lower.split()[0].strip().lower()
            
            # Validate isolation (not too short or just digits)
            if len(isolation) < 2 or isolation.isdigit():
                isolation = None
            
            isolation_match = re.search(isolation_regex, field_lower)
            if isolation_match:
                isolation = isolation_match.group(0).strip().capitalize()
            
            return host, isolation, sex, age_yrs, age_mo, diagnosis
        
        except Exception as e:
            print(f"Error parsing field '{field}': {e}", file=sys.stderr)
            return host, isolation, sex, age_yrs, age_mo, diagnosis
    
    def get_metadata(self, accession: str) -> List[Optional[str]]:
        """
        Fetch metadata for a single accession from GenBank.
        
        Args:
            accession: GenBank accession number
            
        Returns:
            List of metadata fields: [accession, strain, country, location, region, 
                                     subgenogroup, lineage, date, collection_yr, sex, 
                                     age_yrs, age_mo, diagnosis, isolation, origin, doi]
        """
        try:
            # Fetch GenBank record
            handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            handle.close()
            
            # Initialize variables
            country = location = region = subgenogroup = lineage = date = collection_yr = sex = age_yrs = age_mo = diagnosis = isolation = origin = strain = doi = None
            
            # Extract reference information
            if record.annotations.get("references"):
                reference = record.annotations["references"][0]
                description = reference.title
                
                # Check for PubMed ID
                doi = f"PMID: {reference.pubmed_id}" if reference.pubmed_id else "genbank corrections"
                
                # Check for symptoms in description
                for keyword, standardized in self.symptom_list.items():
                    if keyword in description.lower():
                        diagnosis = standardized
                        break
                
                # Extract subgenogroup from description
                genotype_match = re.search(r'\b([A-Za-z0-9\-\.]+)\s*(?:subgenogroup|genogroup|type)\b|\b(?:subgenogroup|genogroup|type)\s*([A-Za-z0-9\-\.]+)\b', description, re.IGNORECASE)
                if genotype_match:
                    subgenogroup = genotype_match.group(1) or genotype_match.group(2)
            
            # Parse source feature
            for feature in record.features:
                if feature.type == "source":
                    organism = feature.qualifiers.get("organism", [None])[0]
                    
                    if organism != self.virus_name:
                        print(f"Accession {accession} is not {self.virus_name} (found: {organism})")
                        return [None] * 16
                    
                    # Extract strain
                    isolate = feature.qualifiers.get("isolate", [None])[0] or None
                    strain = feature.qualifiers.get("strain", [None])[0] or None

                    # Prefer isolate when strain is missing; normalize to safe strings for checks
                    if not strain:
                        strain = isolate
                    isolate_str = isolate or ""
                    strain_str = strain or ""

                    # If both present, prefer the longer identifier (avoid len(None))
                    if isolate_str and strain_str and len(strain_str) < len(isolate_str):
                        strain = isolate

                    # Extract location
                    location_str = feature.qualifiers.get("geo_loc_name", [None])[0]
                    if location_str and ":" in location_str:
                        country, location = location_str.split(":", 1)
                        country = country.strip()
                        location = location.strip()
                    else:
                        country = location_str
                    
                    # Extract host, isolation, and clinical data
                    isolation_field = feature.qualifiers.get("isolation_source", [None])[0]
                    host_field = feature.qualifiers.get("host", [None])[0]
                    note_field = feature.qualifiers.get("note", [None])[0]
                    
                    # Process isolation and host fields
                    if isolation_field or host_field:
                        isolation = isolation_field
                        origin = host_field
                        

                        # Check if isolation_field contains clinical data
                        if isolation_field and ('year' in isolation_field or 'month' in isolation_field or "patient" in isolation_field or "male" in isolation_field or
                            any(key in isolation_field.lower() for key in self.symptom_list)):
                            origin, isolation, sex, age_yrs, age_mo, diagnosis = self.parse_host(isolation_field)
                            origin = host_field
                        
                        
                        # Check host_field if it's not just "Homo sapiens"
                        if host_field and host_field != 'Homo sapiens':
                            origin,isolation, sex, age_yrs, age_mo, diagnosis = self.parse_host(host_field)
                            isolation= isolation_field

                    # get diagnosis from strain name (safe concatenation)
                    for keyword, standardized in self.symptom_list.items():
                        if (strain_str and keyword in strain_str.lower()) or (isolate_str and keyword in isolate_str.lower()):
                            diagnosis = (diagnosis + f"; {standardized}") if diagnosis else standardized
                            break     
                    
                    # get location from strain name (safely use strings)
                    if location is None:
                        for keyword in self.location_list:
                            if (strain_str and keyword.lower() in strain_str.lower()) or (isolate_str and keyword.lower() in isolate_str.lower()):
                                location = keyword
                                break  
                    
                    # Process note field
                    if note_field:
                        if "year" in note_field or "patient" in note_field:
                            origin,isolation, sex, age_yrs, age_mo, diagnosis = self.parse_host(note_field)
                            isolation= isolation_field
                            origin = host_field
                        
                        # Extract severity score
                        if "score" in note_field:
                            severity_matches = re.findall(r"severity\s*score:\s*(\w+).*outcome:\s*(\w+)", note_field)
                            if severity_matches:
                                severity, outcome = severity_matches[0]
                                if diagnosis:
                                    diagnosis += f"; severity {severity}; {outcome}"
                                else:
                                    diagnosis = f"severity {severity}; {outcome}"
                        
                        # Extract genotype from note
                        genotype_match = re.search(r'(type|genogroup)[\s:]*([A-Za-z0-9\-\.]+)', note_field, re.IGNORECASE)
                        if genotype_match:
                            subgenogroup = genotype_match.group(2)
                    
                    diagnosis = diagnosis.replace("HFMD; HFMD", "HFMD")
                    diagnosis = diagnosis.replace("AFP; AFP", "AFP")
                    
                    # Extract collection date
                    date = feature.qualifiers.get("collection_date", [None])[0]

                    # Safely handle missing dates and extract a 4-digit year if present
                    if not date:
                        collection_yr = None
                    else:
                        date_str = str(date).strip()
                        # look for a 4-digit year anywhere in the date string
                        m = re.search(r'(\d{4})', date_str)
                        if m:
                            try:
                                collection_yr = int(m.group(1))
                            except ValueError:
                                collection_yr = m.group(1)
                        else:
                            # fallback: keep original string if no 4-digit year found
                            collection_yr = date_str
                        
            return [accession, strain, country, location, region, subgenogroup, 
                   lineage, date, collection_yr, sex, age_yrs, age_mo, 
                   diagnosis, isolation, origin, doi]
        
        except Exception as e:
            print(f"Error fetching GenBank data for accession {accession}: {e}", file=sys.stderr)
            return [accession] + [None] * 15


def load_config(config_path: str) -> Dict:
    """Load YAML configuration file."""
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)


def load_accessions(input_file: str) -> pd.Series:
    """Load accession numbers from a file (one per line or tab-delimited)."""
    df = pd.read_table(input_file, header=None)
    return df[0].unique()

def normalize_columns_with_mapping(requested, canonical_columns, aliases):
    requested_to_canonical = {}
    canonical_order = []

    for col in requested:
        if col in canonical_columns:
            requested_to_canonical[col] = col
            canonical_order.append(col)
        elif col in aliases:
            requested_to_canonical[col] = aliases[col]
            canonical_order.append(aliases[col])
        else:
            print(
                f"Warning: Requested column '{col}' is not a valid column or alias",
                file=sys.stderr
            )

    return requested_to_canonical, canonical_order


def main():
    """Main function to orchestrate metadata fetching."""
    parser = argparse.ArgumentParser(
        description="Fetch GenBank metadata for viral accessions"
    )
    parser.add_argument(
        '--virus', '-v',
        required=True,
        help='Virus name (e.g., "Coxsackievirus A10")'
    )
    parser.add_argument(
        '--accession_file', '-a',
        required=True,
        help='Path to file containing accession numbers (one per line)'
    )
    parser.add_argument(
        '--output', '-o',
        required=True,
        help='Output TSV file path'
    )
    parser.add_argument(
        '--genbank',
        required=True,
        help='Summary GenBank Metadata TSV file path, includes previous runs'
    )
    parser.add_argument(
        '--latlongs', default='config/lat_longs.tsv',
        help='Path to lat/long TSV file (default: config/lat_longs.tsv)'
    )
    parser.add_argument(
        '--config', '-c',
        default='config/config.yaml',
        help='Path to config YAML file (default: config/config.yaml)'
    )
    parser.add_argument(
        '--email', '-e',
        help='Email for NCBI Entrez (overrides .env file)'
    )
    parser.add_argument(
        '--columns',
        nargs='+',
        default=["accession", "strain", "subgenogroup", "country", "location", 
                "date", "age_yrs", "sex", "diagnosis", "doi"],
        help='Columns to include in output (default: accession strain subgenogroup country location date age_yrs sex diagnosis doi)'
    )
    
    args = parser.parse_args()
    
    # Load email from environment or argument
    load_dotenv(find_dotenv())
    email = args.email or os.environ.get("EMAIL")
    
    if not email:
        print("Error: Email is required. Set EMAIL in .env file or use --email argument", 
              file=sys.stderr)
        sys.exit(1)
    
    # Load configuration
    # print(f"Loading configuration from {args.config}...")
    config = load_config(args.config)
    
    #locations list
    latlongs = pd.read_csv(args.latlongs, sep="\t")
    location_list = latlongs.iloc[:,1].tolist()

    # Load accessions
    # print(f"Loading accessions from {args.accession_file}...")
    accessions = load_accessions(args.accession_file)
    # print(f"Found {len(accessions)} unique accessions")
    
    # Initialize fetcher
    fetcher = MetadataFetcher(args.virus, config, email,location_list)
    
    # Fetch metadata
    # print(f"Fetching metadata for {args.virus}...")
    data = []
    total = len(accessions)
    
    for i, acc in enumerate(accessions, 1):
        print(f"Processing {i}/{total}: {acc}", end='\r')
        metadata = fetcher.get_metadata(acc)
        data.append(metadata)
    
    # print()  # New line after progress
    
    # Create DataFrame
    columns = ["accession", "strain", "country", "location", "region", 
               "subgenogroup", "lineage", "date", "collection_yr", "sex", 
               "age_yrs", "age_mo", "diagnosis", "isolation", "origin", "doi"]

    # Accepted user-facing aliases → canonical name
    COLUMN_ALIASES = {
        "gender": "sex",

        "place": "location",

        "year": "collection_yr",
        
        "collection_date": "date",

        "age": "age_yrs",
        "age_years": "age_yrs",
        "age_month": "age_mo",

        "symptoms": "diagnosis",

        "clade": "subgenogroup",
    }

    requested_to_canonical, canonical_cols = normalize_columns_with_mapping(
        args.columns,
        columns,
        COLUMN_ALIASES
    )

    df = pd.DataFrame(data, columns=columns)

    # Drop rows without accession
    df.dropna(subset=["accession"], inplace=True)

    # Select canonical columns
    df = df.loc[:, canonical_cols]

    # Rename back to requested column names
    df = df.rename(columns={v: k for k, v in requested_to_canonical.items()})
    
    # Save to file
    # print(f"Saving results to {args.output}...")
    df.to_csv(args.output, index=False, sep="\t")

    if os.path.exists(args.genbank):
        gb = pd.read_csv(args.genbank, sep="\t")
    else:
        gb = pd.DataFrame(columns=df.columns)
    
    gb = pd.concat([gb, df], ignore_index=True)

    gb.drop_duplicates(subset=["accession"], keep="last", inplace=True)
    # drop rows if they only have accession, and nothing else
    if gb.shape[1] > 1:
        # Drop rows that have only one non-NA value across all columns
        non_na_counts = gb.notna().sum(axis=1)
        gb = gb.loc[non_na_counts > 1].reset_index(drop=True)
    

    gb.to_csv(args.genbank, index=False, sep="\t")
    
    print(f"Done!  Processed {len(df)} records successfully.")


if __name__ == "__main__":
    main()