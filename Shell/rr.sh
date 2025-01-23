#!/bin/bash

# Loop through each directory in the current location
for dir in */ ; do
    # Check if 'run.sh' exists in the directory
    if [ -f "$dir"run.sh ]; then
        # Go into the directory
        cd "$dir"
        # Execute 'run.sh' with bash
        sbatch run.sh
        # Go back to the parent directory
        cd ..
    else
        echo "run.sh not found in $dir"
    fi
done

