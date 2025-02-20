#!/bin/bash

# Input and output file paths from arguments
input_file="$1"
backup_file="$2"
output_file="$3"

temp_output_file="data/strain_names_temp.tsv"

# Define the positions of accession and strain columns based on the header (1-based index)
accession_col=1
strain_col=5

# Read existing accessions in the backup file into an array for faster lookup
declare -A existing_accessions
if [[ -f "$backup_file" ]]; then
    while IFS=$'\t' read -r accession strain; do
        existing_accessions["$accession"]=1
    done < <(tail -n +2 "$backup_file")  # Skip the header
else
    echo "Error: Backup file not found."
    exit 1
fi

# # Debug: Print the contents of the backup file
# echo "Contents of the backup file:"
# cat "$backup_file"

# # Debug: Print the list of existing accessions
# echo "Existing accessions:"
# for accession in "${!existing_accessions[@]}"; do
#     echo "$accession"
# done

# Initialize skipped and new_accessions counters
skipped=0
new_accessions=0

# Create a list of accessions to process
accessions_to_process=()

# Read the input file and check if accession is equal to strain name
while IFS=$'\t' read -r line; do
    accession=$(echo "$line" | cut -f$accession_col)
    strain=$(echo "$line" | cut -f$strain_col)
    # Debug: Print the current accession and strain
    # echo "Read from input file: accession=$accession, strain=$strain"
    if [[ "$accession" == "$strain" ]]; then
        accessions_to_process+=("$accession")
    fi
done < <(tail -n +2 "$input_file")  # Skip the header

# Debug: Print the list of accessions to process
# echo "Accessions to process: ${accessions_to_process[@]}"

# Remove accessions that are already in the backup file
accessions_to_process=($(comm -23 <(printf "%s\n" "${accessions_to_process[@]}" | sort) <(printf "%s\n" "${!existing_accessions[@]}" | sort)))

# Debug: Print the list of accessions after removing existing ones
echo "Accessions after removing existing ones: ${accessions_to_process[@]}"

# Copy all entries from the backup file to the output file
cp "$backup_file" "$output_file"

# Process the remaining accessions
for accession in "${accessions_to_process[@]}"; do
    # Fetch the GenBank record
    gb_entry=$(efetch -db nucleotide -id "$accession" -format gb 2>/dev/null)

    if [ $? -eq 0 ]; then
        # Extract the strain from the FEATURES section
        strain_from_gb=$(echo "$gb_entry" | grep -oP '/strain="\K[^"]+')

        if [ -n "$strain_from_gb" ]; then
            strain=$strain_from_gb
            ((new_accessions += 1))
            echo "Updated strain for $accession: $strain"
        else
            strain="$accession"
        fi
    else
        strain="$accession"
    fi

    # Append the updated entry to the output file
    echo -e "$accession\t$strain" >> "$output_file"
done

# Sort and remove duplicates in the output file
sort -u "$output_file" -o "$output_file"

# Add the header if not already present
if ! head -n 1 "$output_file" | grep -qP '^accession\s+strain'; then
    # Prepend the header to the output file
    echo -e "accession\tstrain" | cat - "$output_file" > temp && mv temp "$output_file"
fi

echo -e "\nProcessing completed. Updated metadata saved to $output_file.\n$skipped accessions were skipped.\n$new_accessions new strains were added."