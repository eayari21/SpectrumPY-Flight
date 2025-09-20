#!/bin/bash

# Define the base directories
# directories=("../Pre_Env" "../Post_Env")
directories=("Data/")
hdf5_output_dir="./HDF5"  # Directory where .h5 files are stored

# Iterate through each directory
for dir in "${directories[@]}"; do
    echo "Processing directory: $dir"
    # Find all files with no extension and iterate through them
    find "$dir" -type f ! -name "*.*" | while read -r file; do
        # Extract the base filename without the path or extension
        base_filename=$(basename "$file")
        # Define the expected HDF5 output file path
        output_hdf5_file="$hdf5_output_dir/$base_filename.h5"
        
        # Check if the .h5 file already exists
        if [[ -f "$output_hdf5_file" ]]; then
            echo "Output file $output_hdf5_file already exists. Skipping processing for $file."
            continue
        fi

        echo "Processing file: $file"
        # Call the Python script with the file and handle errors
        python lmfit_idex_packet.py -f "$file" || {
            echo "Error processing file: $file. Skipping."
            continue
        }
    done
done


