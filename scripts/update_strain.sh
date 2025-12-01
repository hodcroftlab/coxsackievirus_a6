#!/bin/bash

# Input and output file paths from arguments
input_file="$1"
backup_file="$2"
output_file="$3"

# Define column indices (1-based)
accession_col=1
strain_col=5

# Exit if required files are missing
if [[ ! -f "$input_file" ]]; then
    echo "Error: Input file not found: $input_file"
    exit 1
fi

if [[ ! -f "$backup_file" ]]; then
    echo "Error: Backup file not found: $backup_file"
    exit 1
fi

# Read existing accessions into associative array
declare -A existing_accessions
while IFS=$'\t' read -r accession strain; do
    existing_accessions["$accession"]=1
done < <(tail -n +2 "$backup_file")

# Create list of accessions where accession == strain
accessions_to_process=()
while IFS=$'\t' read -r line; do
    accession=$(echo "$line" | awk -F'\t' -v col=$accession_col '{print $col}')
    strain=$(echo "$line" | awk -F'\t' -v col=$strain_col '{print $col}')
    if [[ "$accession" == "$strain" ]]; then
        accessions_to_process+=("$accession")
    fi
done < <(tail -n +2 "$input_file")

# Remove accessions already in backup
accessions_to_process=($(comm -23 <(printf "%s\n" "${accessions_to_process[@]}" | sort) <(printf "%s\n" "${!existing_accessions[@]}" | sort)))

num_to_process=${#accessions_to_process[@]}
echo "Found $num_to_process new accessions with identical strain names."

# Copy backup to output
cp "$backup_file" "$output_file"

# Ensure backup file ends with a newline
tail -c1 "$output_file" | read -r _ || echo >> "$output_file"

# Process and fetch missing strain names
new_accessions=0
for accession in "${accessions_to_process[@]}"; do
    gb_entry=$(efetch -db nucleotide -id "$accession" -format gb 2>/dev/null)
    if [[ $? -eq 0 ]]; then
        strain_from_gb=$(echo "$gb_entry" | grep -oP '/strain="\K[^"]+')
        if [[ -n "$strain_from_gb" ]]; then
            strain="$strain_from_gb"
            ((new_accessions++))
        else
            strain="$accession"
        fi
    else
        strain="$accession"
    fi
    echo -e "$accession\t$strain" >> "$output_file"
done

# Sort and remove duplicates
sort -u "$output_file" -o "$output_file"

# Ensure header is present
if ! head -n 1 "$output_file" | grep -qP '^accession\s+strain'; then
    echo -e "accession\tstrain" | cat - "$output_file" > temp && mv temp "$output_file"
fi

echo -e "\nDone. $new_accessions updated strains written to: $output_file"