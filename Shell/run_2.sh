#!/bin/bash

# Base directory containing the folders
base_dir="/nfshome/deshmukh/vaibhav/NaSICON_dft/database_matropolis_extra_nasi"

echo "Starting script in base directory: $base_dir"

# Check if base_dir exists
if [[ ! -d "$base_dir" ]]; then
    echo "Base directory does not exist: $base_dir"
    exit 1
fi

# Loop through each subdirectory in the base directory
for dir in "$base_dir"/*/; do
    # Get the directory name without the path
    dir_name=$(basename "${dir%/}")
    
    echo "Processing directory: $dir"

    # Construct the tar.bz2 filename based on the folder name
    tar_file="${dir}${dir_name}.tar.bz2"

    # Check if the tar.bz2 file exists in the directory
    if [[ -f "$tar_file" ]]; then
        echo "Found archive: $tar_file"

        # Change into the directory where the tar file is located
        cd "$dir" || { echo "Failed to change directory to $dir"; continue; }

        # Extract OUTCAR from the tar.bz2 file
        tar -xvjf "$tar_file" OUTCAR
        echo "Extracted OUTCAR into $dir"

        # Change back to the base directory
        cd "$base_dir" || { echo "Failed to change back to base directory"; exit 1; }
    else
        echo "Archive not found: $tar_file"
    fi
done

echo "Script completed."

