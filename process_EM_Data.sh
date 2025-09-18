#!/bin/bash

# Set the base directories to search through
BASE_DIRS=(
    "/Users/etay8828/Desktop/IDEX_Cal/Pre_Launch/2025_02_06_run2_split"
    "/Users/etay8828/Desktop/IDEX_Cal/Pre_Launch/2025_02_06_run1_split"
)

# Log file for errors
LOG_FILE="error_log.txt"

# Clear the log file at the start
> "$LOG_FILE"

# Function to process folders
process_folder() {
    local folder="$1"
    # Check if the folder contains .trc files
    if ls "$folder"/*.trc 1> /dev/null 2>&1; then
        # Extract the folder name to use as the experiment name
        experiment_name=$(basename "$folder")
        echo "Processing folder: $folder with experiment name: $experiment_name"

        # Run the script with the folder as the --trcdir argument
        python ImpactBook.py --trcdir "$folder" --experimentname "$experiment_name"
        if [ $? -ne 0 ]; then
            echo "Error processing $folder" >> "$LOG_FILE"
        fi
    fi
}

# Export the function for subshells
export -f process_folder
export LOG_FILE

# Loop through both base directories
for BASE_DIR in "${BASE_DIRS[@]}"; do
    echo "Processing base directory: $BASE_DIR"
    find "$BASE_DIR" -type d -exec bash -c 'process_folder "$0"' {} \;
done

echo "Processing complete. Check $LOG_FILE for errors."
