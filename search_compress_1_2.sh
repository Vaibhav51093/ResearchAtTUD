#!/bin/bash

# Find all dump_nvt_prod_5.out files in the specified directory
find  /nfshome/deshmukh -type f -name "dump_nvt_prod_5.out" 2>/dev/null | while read -r file; do
    dir=$(dirname "$file")
    echo "Processing directory '${dir}'"	
    # Go into the directory
    pushd "$dir" > /dev/null

    # Check if the file is in the tar.gz
    if [[ -e "dump_nvt_prod_5.out.tar.gz" ]]; then
	    echo "'dump_nvt_prod_5.out.tar.gz' already exists in directory '${dir}'"
        if ! tar -ztf dump_nvt_prod_5.out.tar.gz | grep -q 'dump_nvt_prod_5.out'; then
            # If not, add it to a new tar file and compress it
            tar -cf dump_nvt_prod_5.out.tar dump_nvt_prod_5.out
            gzip -f dump_nvt_prod_5.out.tar
            echo "Compressed 'dump_nvt_prod_5.out' in directory '${dir}'"
        fi
    else
        # If the tar.gz file doesn't exist, create a new tar file and compress it
        tar -cf dump_nvt_prod_5.out.tar dump_nvt_prod_5.out
        gzip -f dump_nvt_prod_5.out.tar
        echo "Compressed 'dump_nvt_prod_5.out' in directory '${dir}'"
    fi

    # Always delete the original file if it exists
    if [[ -e "dump_nvt_prod_5.out" ]]; then
        rm dump_nvt_prod_5.out
        echo "Deleted original 'dump_nvt_prod_5.out' in directory '${dir}'"
    fi

    # Return to the original directory
    popd > /dev/null
done

