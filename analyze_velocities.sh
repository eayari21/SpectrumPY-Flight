#!/bin/bash

# Base directory for Scope_Data
BASE_DIR="/Users/etay8828/Desktop/IDEX_Cal/Scope"

# Python script path
PYTHON_SCRIPT="qd_quicklook.py"

# Loop through all subdirectories in Scope_Data
find "$BASE_DIR" -type d | while read -r DIR; do
  # Check if the directory contains any .trc files
  if ls "$DIR"/*.trc >/dev/null 2>&1; then
    # Extract the last part of the directory path for the -t argument
    LAST_PART=$(basename "$DIR")

    # Call the Python script with the directory and last part of the path
    echo "Processing: $DIR"
    python "$PYTHON_SCRIPT" -s "$DIR" -t "$LAST_PART"
  fi
done
