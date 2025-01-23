#!/bin/bash

# Base directory to search from
BASE_DIR="/nfshome/deshmukh/vaibhav"

# Create directories for categorized files
mkdir -p PythonFiles NotebookFiles ShellScripts

# Iterate over all directories and files in the base directory
for file in $(find "$BASE_DIR" -type f); do
    # Check the file extension and copy to the respective folder
    case "$file" in
        *.py)
            cp "$file" PythonFiles/ || echo "Failed to copy $file"
            ;;
        *.ipynb)
            cp "$file" NotebookFiles/ || echo "Failed to copy $file"
            ;;
        *.sh)
            cp "$file" ShellScripts/ || echo "Failed to copy $file"
            ;;
    esac
done

# Print completion message
echo "Files have been copied into PythonFiles, NotebookFiles, and ShellScripts folders."

