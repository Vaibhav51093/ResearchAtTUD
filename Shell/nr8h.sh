#!/bin/bash 

# Here the aim is to check the directory using find command 

find  /nfshome/deshmukh/vaibhav/frank_jÜlich/lmp_3_1_3/large/simu/simu_2/773/new/ -type f \( -name "dump_nvt_prod.out" -o -name " dump_nvt_prod.out.tar.gz"\) 2>/dev/null | while read -r file;directory; do
    echo "file: $file"
    echo "directory: $directory"

    pushd "$directory" >/dev/null

    # if both files exist, then extract dump_nvt_prod.out from dump_nvt_prod.out.tar.gz, delete dump_nvt_prod.out.tar.gz and compress dump_nvt_prod.out to dump_nvt_prod.out.gz and delete dump_nvt_prod.out
    if [[ -f dump_nvt_prod.out && -f dump_nvt_prod.out.tar.gz ]]; then
        tar -xzf dump_nvt_prod.out.tar.gz dump_nvt_prod.out
        rm dump_nvt_prod.out.tar.gz
        gzip dump_nvt_prod.out
        rm dump_nvt_prod.out
    fi

done