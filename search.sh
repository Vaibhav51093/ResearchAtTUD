#!/bin/bash


find /nfshome/deshmukh/pyiron/projects/NASICON/project/ \( -type f -name "dump_nvt_prod.out" -o -name "dump_nvt_prod.out.tar.gz" \) 2>/dev/null | while read -r file; do
    directory=$(dirname "$file")
    echo "file: $file"
    echo "directory: $directory"

    pushd "$directory" >/dev/null

    # Check if both files exist
    if [[ -f "dump_nvt_prod.out" && -f "dump_nvt_prod.out.tar.gz" ]]; then
        tar -xzf "dump_nvt_prod.out.tar.gz" "dump_nvt_prod.out"
        rm "dump_nvt_prod.out.tar.gz"
        gzip -f "dump_nvt_prod.out"
        rm "dump_nvt_prod.out"
        echo "Extracted 'dump_nvt_prod.out' from 'dump_nvt_prod.out.tar.gz', compressed to 'dump_nvt_prod.out.gz', and deleted both files"
    elif [[ -f "dump_nvt_prod.out.tar.gz" ]]; then
        # Extract dump_nvt_prod.out from dump_nvt_prod.out.tar.gz
        tar -xzf "dump_nvt_prod.out.tar.gz" "dump_nvt_prod.out"
        gzip -f "dump_nvt_prod.out"
        rm "dump_nvt_prod.out.tar.gz"
        echo "Extracted 'dump_nvt_prod.out' from 'dump_nvt_prod.out.tar.gz' and compressed to 'dump_nvt_prod.out.gz'"
    elif [[ -f "dump_nvt_prod.out" ]]; then
        # Compress dump_nvt_prod.out to dump_nvt_prod.out.gz
        gzip -f "dump_nvt_prod.out"
        rm "dump_nvt_prod.out"
        echo "Compressed 'dump_nvt_prod.out' to 'dump_nvt_prod.out.gz' and deleted the file"
    fi

    popd >/dev/null
done





