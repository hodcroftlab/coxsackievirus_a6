#!/bin/bash

# Input and output file paths from arguments
input_file="$1"
output_file="$2"

# Define the columns to keep in the output file
required_columns=("accession" "strain")

# Read the header line from the input file
header_line=$(head -n 1 "$input_file")

# Extract header names into an array
IFS=$'\t' read -r -a headers <<< "$header_line"

# Create a mapping of header names to their indices
declare -A header_indices
for i in "${!headers[@]}"; do
    header_indices[${headers[$i]}]=$i
done

# Identify the indices of the required columns
for col in "${required_columns[@]}"; do
    if [[ -z "${header_indices[$col]}" ]]; then
        echo "Error: Required column '$col' not found in the header."
        exit 1
    fi
done

# Create a temporary output file for new entries
temp_output_file="$output_file.tmp"

# Create or clear the output file if it doesn't exist and write the header
if [[ ! -f "$output_file" ]]; then
    required_indices=$(for col in "${required_columns[@]}"; do echo -n "${header_indices[$col]} "; done)
    echo -e "$(for idx in $required_indices; do echo -n "${headers[$idx]}\t"; done | sed 's/\t$//')" > "$output_file"
fi

# Load existing accessions from the output file into a set
declare -A existing_accessions
while IFS=$'\t' read -r -a fields; do
    accession="${fields[0]}"  # Assumes accession is the first column
    existing_accessions["$accession"]=1
done < <(tail -n +2 "$output_file")

# Process the input file and append new entries to the temp file
while IFS=$'\t' read -r -a fields; do
    # Skip empty lines
    if [[ -z "${fields[${header_indices[accession]}]}" ]]; then
        continue
    fi

    # Extract values based on header indices
    accession="${fields[${header_indices[accession]}]}"
    strain="${fields[${header_indices[strain]}]}"

    # Check if accession is already in the output file
    if [[ -n "${existing_accessions[$accession]}" ]]; then
        echo "Skipping $accession, already present in the output file."
        continue
    fi

    # Fetch and update the strain if needed
    if [ "$accession" == "$strain" ]; then
        # Fetch the GenBank record
        gb_entry=$(efetch -db nucleotide -id "$accession" -format gb 2>/dev/null)
        
        if [ $? -eq 0 ]; then
            # Extract the strain from the FEATURES section
            strain_from_gb=$(echo "$gb_entry" | grep -oP '/strain="\K[^"]+')
            
            if [ -n "$strain_from_gb" ]; then
                strain=$strain_from_gb
            else
                strain="UNKNOWN"
            fi
        else
            strain="UNKNOWN"
        fi
    fi

    # Replace the strain value in the fields array
    fields[${header_indices[strain]}]=$strain

    # Write the selected fields to the temp output file
    updated_row=$(for idx in $required_indices; do echo -n "${fields[$idx]}\t"; done | sed 's/\t$//')
    echo -e "$updated_row" >> "$temp_output_file"
done < <(tail -n +2 "$input_file")

# Append the temp output to the main output file and remove duplicates
if [[ -f "$temp_output_file" ]]; then
    cat "$temp_output_file" >> "$output_file"
    sort -u "$output_file" -o "$output_file"
    rm "$temp_output_file"
fi

# Report completion
echo "Processing completed. Updated metadata saved to $output_file"