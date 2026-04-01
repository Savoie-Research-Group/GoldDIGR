#!/bin/bash

# --- Environment Setup ---
# Activate the same environment as the master workflow to ensure
# 'calculate_mult.py' is available and runs correctly.
echo "Activating Conda environment..."
module load anaconda/2022.10-py39
#source activate /scratch/negishi/li1724/100624_copy_yarp/
conda activate copy-classy-yarp
#conda activate IOP
module load intel-mkl
module load openmpi

# --- Configuration ---
# File containing the absolute paths to the XYZ files, one per line.
# This should be the SAME file used by 5_master_workflow.sh
XYZ_LIST_FILE="AAAS_folders.txt"

# Output files for a summary of job states
FINISHED_FILE="finished.txt"
CRASHED_FILE="crashed.txt"
UNTOUCHED_FILE="untouched.txt"

# --- Script Body ---
if [ ! -f "$XYZ_LIST_FILE" ]; then
    echo "Error: XYZ list file not found at '$XYZ_LIST_FILE'"
    exit 1
fi

# Ensure calculate_mult.py is findable
if ! command -v python &> /dev/null || ! python -c "import sys; sys.path.append('.'); import calculate_mult" &> /dev/null; then
    if [ ! -f "calculate_mult.py" ]; then
        echo "Error: 'calculate_mult.py' not found in the current directory."
        echo "Please run this script from the same directory as 'calculate_mult.py'."
        exit 1
    fi
fi

echo "Starting job status check..."
echo "Input file: $XYZ_LIST_FILE"
echo "---"

# Clear/initialize output files
> "$FINISHED_FILE"
> "$CRASHED_FILE"
> "$UNTOUCHED_FILE"

# Initialize counters
count_finished=0
count_crashed=0
count_untouched=0

while IFS= read -r DIR_PATH || [[ -n "$DIR_PATH" ]]; do
    if [ -z "$DIR_PATH" ]; then
        continue
    fi

    # PATH_NO_SLASH is the base path, e.g., /path/to/block/16
    PATH_NO_SLASH=${DIR_PATH%/}
    # XYZ_FILE_PATH is the input geometry file, e.g., /path/to/block/16.xyz
    XYZ_FILE_PATH="${PATH_NO_SLASH}.xyz"

    if [ ! -f "$XYZ_FILE_PATH" ]; then
        echo "  -> WARNING: Input XYZ file not found: $XYZ_FILE_PATH. Skipping all charge states for this entry."
        continue
    fi
    
    echo "Checking base path: $PATH_NO_SLASH"

    # Loop over the three charge states
    for charge in -1 0 1; do
        # Calculate the multiplicity just as the master script does
        multiplicity=$(python calculate_mult.py "$XYZ_FILE_PATH" "$charge")

        if [ $? -ne 0 ]; then
            echo "    -> ERROR: Failed to calculate multiplicity for $XYZ_FILE_PATH (Charge: $charge). Skipping this job."
            continue
        fi

        # This is the full path to the expected job directory, e.g., /path/to/block/16_-1_2
        FULL_DIR_PATH="${PATH_NO_SLASH}_${charge}_${multiplicity}"
        
        # Define the paths for the zip file and the job directory
        ZIP_FILE="${FULL_DIR_PATH}.zip"
        JOB_DIR="${FULL_DIR_PATH}"

        # --- Status Check Logic ---

        # 1. Finished: The .zip file exists.
        if [ -f "$ZIP_FILE" ]; then
            echo "    -> FINISHED: ${ZIP_FILE}"
            echo "$FULL_DIR_PATH" >> "$FINISHED_FILE"
            ((count_finished++))
            
        # 2. Crashed: The .zip file does NOT exist, but the job directory DOES.
        #    (This implies the job started but did not finish and zip)
        elif [ -d "$JOB_DIR" ]; then
            echo "    -> CRASHED:  ${JOB_DIR}"
            echo "$FULL_DIR_PATH" >> "$CRASHED_FILE"
            ((count_crashed++))

        # 3. Untouched: Neither the .zip file nor the job directory exists.
        else
            echo "    -> UNTOUCHED: ${FULL_DIR_PATH}"
            echo "$FULL_DIR_PATH" >> "$UNTOUCHED_FILE"
            ((count_untouched++))
        fi
    done

done < "$XYZ_LIST_FILE"

echo "---"
echo "Job status check finished."
echo "Summary:"
echo "  Finished:  $count_finished (see $FINISHED_FILE)"
echo "  Crashed:   $count_crashed (see $CRASHED_FILE)"
echo "  Untouched: $count_untouched (see $UNTOUCHED_FILE)"

