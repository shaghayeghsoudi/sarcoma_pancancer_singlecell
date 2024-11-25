#!/bin/bash

# Loop through all directories in the current location
for dir in */; do
    # Extract the part of the directory name after the first "_"
    dir_name="${dir%/}"  # Remove trailing slash
    prefix="${dir_name#*_}"  # Remove everything up to and including the first "_"

    # Check if a prefix was extracted
    if [[ -n "$prefix" ]]; then
        echo "Processing directory: $dir_name with prefix: $prefix"

        # Loop through all files in the directory
        for file in "$dir_name"/*; do
            # Extract the base name of the file
            file_name=$(basename "$file")

            # Prepend the prefix to the file name, separated by "_"
            mv "$file" "$dir_name/${prefix}_$file_name"
        done
    else
        echo "Skipping directory: $dir_name (no prefix found)"
    fi
done
