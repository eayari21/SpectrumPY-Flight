#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PYTHONPATH="${SCRIPT_DIR}/src${PYTHONPATH:+:${PYTHONPATH}}"
PYTHON_BIN="${PYTHON:-python}"

# Collect CLI-specified inputs (files or directories)
inputs=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        -f|--file|--filename)
            if [[ -z "$2" ]]; then
                echo "Error: Missing value for $1." >&2
                exit 1
            fi
            inputs+=("$2")
            shift 2
            ;;
        *)
            inputs+=("$1")
            shift
            ;;
    esac
done

if [[ ${#inputs[@]} -eq 0 ]]; then
    # Default to historical behaviour of scanning the Data directory
    inputs=("Data/")
fi
hdf5_output_dir="./HDF5"  # Directory where .h5 files are stored

process_file() {
    local file="$1"
    if [[ ! -f "$file" ]]; then
        echo "Skipping missing file: $file"
        return
    fi

    local base_filename
    base_filename=$(basename "$file")
    local output_hdf5_file="$hdf5_output_dir/$base_filename.h5"

    if [[ -f "$output_hdf5_file" ]]; then
        echo "Output file $output_hdf5_file already exists. Skipping processing for $file."
        return
    fi

    echo "Processing file: $file"
    "$PYTHON_BIN" -m spectrumpy_flight.idex_packet -f "$file" || {
        echo "Error processing file: $file. Skipping."
        return
    }
}

for target in "${inputs[@]}"; do
    if [[ -d "$target" ]]; then
        echo "Processing directory: $target"
        find "$target" -type f ! -name "*.*" | while read -r file; do
            process_file "$file"
        done
    elif [[ -f "$target" ]]; then
        process_file "$target"
    else
        echo "Skipping unknown path (not a file or directory): $target"
    fi
done


