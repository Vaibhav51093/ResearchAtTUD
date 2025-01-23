#!/bin/bash

# Find all dump_nvt_prod.out.tar.bz files in the specified directory
find /nfshome/deshmukh/pyiron/projects/NASICON/project/ -type f -name "dump_nvt_prod.out.tar.bz" 2>/dev/null | while read -r file; do
    dir=$(dirname "$file")
    echo "Processing directory '${dir}'"

    # Check if the file is in the tar.bz2 format
    if [[ -e "dump_nvt_prod.out.tar.bz2" ]]; then
        echo "'dump_nvt_prod.out.tar.bz2' already exists in directory '${dir}'"
    else
        # Extract the file from the tar.bz and compress it to tar.bz2
        tar -xf "$file" -C "$dir"
        tar -cjf "dump_nvt_prod.out.tar.bz2" -C "$dir" "dump_nvt_prod.out"
        echo "Compressed 'dump_nvt_prod.out' to 'dump_nvt_prod.out.tar.bz2' in directory '${dir}'"

        # Remove the original files
        rm "dump_nvt_prod.out"
        rm "$file"
        echo "Deleted original files in directory '${dir}'"
    fi
done

