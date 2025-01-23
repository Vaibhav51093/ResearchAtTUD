#!/bin/bash


find /nfshome/deshmukh/vaibhav/frank_jÜlich/lmp_3_1_3/large/simu/simu_2/773/new/ -type f -name "dump_nvt_prod.out" -o -name "dump_nvt_prod.out.tar.gz" 2>/dev/null | while read -r file; do
    directory=$(dirname "$file")
    echo "file: $file"
    echo "directory: $directory"

    pushd "$directory" >/dev/null

    if [[ -f "dump_nvt_prod.out.tar.gz" ]]; then
        tar -xzf "dump_nvt_prod.out.tar.gz" "dump_nvt_prod.out" && rm "dump_nvt_prod.out.tar.gz"
    fi

    if [[ -f "dump_nvt_prod.out" ]]; then
        gzip -f "dump_nvt_prod.out" && rm "dump_nvt_prod.out"
    fi

    popd >/dev/null
done






