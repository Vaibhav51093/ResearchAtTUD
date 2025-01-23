#!/bin/bash


find /nfshome/deshmukh/vaibhav/frank_jÜlich/lmp_3_1_3/large/simu/simu_2/773/new/ -type f \( -name "dump_nvt_prod.out" -o -name "dump_nvt_prod.out.tar.gz" \) 2>/dev/null | while read -r file; do
    directory=$(dirname "$file")
    echo "file: $file"
    echo "directory: $directory"

    pushd "$directory" >/dev/null

    # If both files exist, then extract dump_nvt_prod.out from dump_nvt_prod.out.tar.gz, compress it to dump_nvt_prod.out.gz,
    # and delete both files
    if [[ -f "dump_nvt_prod.out" && -f "dump_nvt_prod.out.tar.gz" ]]; then
        tar -xzf "dump_nvt_prod.out.tar.gz" "dump_nvt_prod.out"
        rm "dump_nvt_prod.out.tar.gz"
        gzip "dump_nvt_prod.out"
        rm "dump_nvt_prod.out"
        echo "Extracted 'dump_nvt_prod.out' from 'dump_nvt_prod.out.tar.gz', compressed to 'dump_nvt_prod.out.gz', and deleted both files"
    # If only dump_nvt_prod.out.tar.gz exists, extract dump_nvt_prod.out from it, compress it to dump_nvt_prod.out.gz,
    # and delete both files
    elif [[ -f "dump_nvt_prod.out.tar.gz" ]]; then
        tar -xzf "dump_nvt_prod.out.tar.gz" "dump_nvt_prod.out"
        gzip "dump_nvt_prod.out"
        rm "dump_nvt_prod.out.tar.gz"
        rm "dump_nvt_prod.out"
        echo "Extracted 'dump_nvt_prod.out' from 'dump_nvt_prod.out.tar.gz', compressed to 'dump_nvt_prod.out.gz', and deleted both files"
    # If only dump_nvt_prod.out exists, compress it to dump_nvt_prod.out.gz and delete the file
    elif [[ -f "dump_nvt_prod.out" ]]; then
        gzip "dump_nvt_prod.out"
        rm "dump_nvt_prod.out"
        echo "Compressed 'dump_nvt_prod.out' to 'dump_nvt_prod.out.gz' and deleted the file"
    fi

    popd >/dev/null
done



