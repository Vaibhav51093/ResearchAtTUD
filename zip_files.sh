#!/bin/bash

# Function to gzip files with specified extensions
gzip_files() {
    local dir="$1"
    # Find .out and .lammps files that are not already gzipped and gzip them
    find "$dir" -type f \( -name "*.out" -o -name "*.lammps" \) ! -name "*.gz" -exec gzip {} \;
}

# Function to recursively process directories
process_directory() {
    local dir="$1"
    
    # Process files in the current directory
    gzip_files "$dir"
    
    # Recursively process each subdirectory
    for sub_dir in "$dir"/*; do
        if [ -d "$sub_dir" ]; then
            process_directory "$sub_dir"
        fi
    done
}

# Main script
main() {
    # Starting directory (can be changed if needed)
    local start_dir="."

    # Check if a directory argument is given
    if [ $# -gt 0 ]; then
        start_dir="$1"
    fi

    # Ensure the starting directory exists
    if [ ! -d "$start_dir" ]; then
        echo "The directory $start_dir does not exist."
        exit 1
    fi

    # Process the starting directory
    process_directory "$start_dir"
}

# Execute main function with all script arguments
main "$@"

