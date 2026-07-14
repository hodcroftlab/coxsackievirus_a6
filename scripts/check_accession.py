import requests
import re
from Bio import Entrez, SeqIO
import pandas as pd
from dotenv import load_dotenv, find_dotenv
import os

import time

# Function to check if an accession number is real: it uses the entrez functionality of ncbi
def extract_accession(name, extract="accession"):
    # Try to load from .env file (for local development)
    load_dotenv(find_dotenv())
    
    # Get email from environment variable (works for both local .env and GitHub Actions)
    email = os.environ.get("EMAIL") or os.environ.get("ENTREZ_EMAIL")
    
    if not email:
        raise ValueError("EMAIL or ENTREZ_EMAIL environment variable not set")
    
    Entrez.email = email
    # email address should be stored in an .env file in the base directory, in the format EMAIL=user@email.com
    # OR as a GitHub Secret named ENTREZ_EMAIL
    
    # Define the API endpoint and parameters
    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'
    params = {'db': 'nuccore', 'term': name, 'retmode': 'xml'}
    
    max_retries = 5
    backoff_factor = 5 #seconds

    time.sleep(1) # for every entry sleep first

    for attempt in range(max_retries):
        try:
            # Send request to the API
            response = requests.get(url, params=params)
            response.raise_for_status()  # Raise HTTPError for bad responses

            # Check if accession exists
            if '<Count>0</Count>' in response.text:
                return False, None

            # Extract all IDs from the response
            ids = re.findall(r'<Id>(\d+)</Id>', response.text)
            if len(ids) > 5:
                return False, None
            accessions = []

            for id_name in ids:
                # Fetch the sequence record using Entrez
                handle = Entrez.efetch(db="nucleotide", id=id_name, rettype="gb", retmode="text")
                seq_record = SeqIO.read(handle, "genbank")
                handle.close()  # Ensure we close the handle after reading

                # Check if the description contains "Coxsackievirus A16" or "CVA16"
                if re.search(r'\b(Coxsackievirus A16|CVA16)\b', seq_record.description, re.IGNORECASE):
                    # Check if the description also contains "VP1" or "complete genome"
                    if re.search(r'\bVP1\b', seq_record.description, re.IGNORECASE) or re.search(r'\bcomplete \b', seq_record.description, re.IGNORECASE) or re.search(r'\bpartial cds\b', seq_record.description, re.IGNORECASE):
                        # what to extract
                        if extract == "accession":
                            accession = seq_record.name
                        elif extract == "strain":
                            match = re.search(r'(strain|isolate)\s+([\w\-\.\/]+)', seq_record.description)
                            if match:
                                strain_name = match.group(2)
                                if re.search(r'\bcomplete \b', seq_record.description, re.IGNORECASE) and len(ids) > 1:
                                    accession = strain_name + "c"
                                else:
                                    accession = strain_name
                        else:
                            accession = seq_record.description
                        accessions.append(accession)

            if accessions:
                # Safely unpack accessions if it has at least one item
                accession = accessions[0]
                return True, accession
            else:
                return False, None

        except requests.exceptions.HTTPError as e:
            if response.status_code == 429:
                # Too many requests, wait and retry
                wait_time = backoff_factor * (2 ** attempt)
                print(f"Rate limit exceeded. Retrying in {wait_time} seconds...")
                time.sleep(wait_time)
            else:
                print(f"HTTP error occurred: {e}")
                break
        except Exception as e:
            # Handle other exceptions gracefully
            print(f"Error checking accession {name}: {e}")
            break

    return False, None


def extract_digits(string, n=4, position='first'):
    # Extract all digits from the string
    digits = ''.join(filter(str.isdigit, string))
    
    if position == 'first':
        # Return the first 'n' digits, or fewer if not enough digits are present
        return digits[:n]
    elif position == 'last':
        # Return the last 'n' digits, or fewer if not enough digits are present
        return digits[-n:] if len(digits) >= n else digits
    else:
        raise ValueError("Position should be 'first' or 'last'")
