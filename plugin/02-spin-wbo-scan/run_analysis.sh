#!/bin/bash

# --- Safety Options ---
# -e: Exit immediately if a command fails
# -u: Treat unset variables as an error
# -o pipefail: Exit if any command in a pipeline fails
set -euo pipefail

module load conda/2025.02
conda activate another-yarp
#conda activate IOP
module load intel-mkl
module load openmpi

# --- Main Configuration ---
# The file to check for inside the zip archives
CHECK_FILE="finished_irc.trj"

# The directory to search for zips.
# Uses the first script argument (e.g., ./run_analysis.sh /test)
# or defaults to the current directory ('.') if no argument is given.
SEARCH_DIR="${1:-.}"

# --- Python Analysis Configuration ---
# The name of your python script.
# IMPORTANT: This script MUST be in the same directory as run_analysis.sh
PYTHON_SCRIPT_NAME="spin_wbo_scan.py"

# Arguments for the python script
ARG_CHARGE=1
ARG_MULT=-1
ARG_UHF_MAX=6
ARG_NPROCS=1
# The python script will create this sub-directory inside the temp folder
ARG_WORKDIR="xTB-scan" 

# --- End Configuration ---

# Get the absolute path of the directory this script (run_analysis.sh) is in
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
PYTHON_SCRIPT_PATH="$SCRIPT_DIR/$PYTHON_SCRIPT_NAME"

# Check if the python script actually exists before we start
if [ ! -f "$PYTHON_SCRIPT_PATH" ]; then
    echo "🚨 Error: Python script not found!"
    echo "Please place '$PYTHON_SCRIPT_NAME' in the same directory as this script."
    echo "Expected location: $PYTHON_SCRIPT_PATH"
    exit 1
fi
echo "🐍 Python script found: $PYTHON_SCRIPT_PATH"

# Get the absolute path of the directory we're running from
START_DIR=$(pwd)
# Get the absolute path of the directory to search
ABS_SEARCH_DIR=$(readlink -f "$SEARCH_DIR")

echo "🔍 Starting processing in: $ABS_SEARCH_DIR"
echo "---"

# Change to the search directory
cd "$ABS_SEARCH_DIR"

# Loop through all files ending in .zip
# Loop through all files ending in .zip
for zip_file in ./*.zip; do
    
    # Make sure it's a file before processing
    [ -f "$zip_file" ] || continue

    # Skip zip files whose name indicates a crash
    if [[ "$zip_file" == *crash* ]]; then
        echo "Skipping (marked crash): $zip_file"
        continue
    fi

    # Derive charge from zip basename: <...>_<charge>_<mult>.zip
    # Examples:
    #   38_C-TSa1_-1_2.zip  -> charge = -1
    #   38_C-TSa1_0_1.zip    -> charge = 0
    zip_base="${zip_file##*/}"   # strip leading ./ or path → 38_C-TSa1_-1_2.zip
    zip_base="${zip_base%.zip}"  # drop extension         → 38_C-TSa1_-1_2
    tmp="${zip_base%_*}"         # drop last _mult        → 38_C-TSa1_-1
    ARG_CHARGE="${tmp##*_}"      # take last token        → -1 or 0
    echo "  Charge parsed from name: $ARG_CHARGE"

    # Get the absolute path to the zip file
    abs_zip_file=$(readlink -f "$zip_file")
    
    # Define the analysis directory name by removing the '.zip' extension
    # e.g., ./1-1-0.zip -> ./1-1-0
    analysis_dir="${zip_file%.zip}"
    abs_analysis_dir=$(readlink -f "$analysis_dir")

    echo "Processing: $abs_zip_file"

    # 0. Skip already-processed zip files (those that already contain xTB-scan/)
    if { unzip -l "$abs_zip_file" || true; } | grep -q "xTB-scan/"; then
        echo "  Skipping: xTB-scan folder already present in zip (already processed)."
        echo "---"
        continue
    fi

    # 1. Check if the target file exists in the zip (lists without extracting)
    if { unzip -l "$abs_zip_file" || true; } | grep -q "$CHECK_FILE"; then
        
        echo "  ✅ Found $CHECK_FILE. Extracting for analysis..."

        # 2. Create directory and extract ONLY the required file
	# 2.1. Clean any stale analysis dir and create it fresh
        if [ -d "$abs_analysis_dir" ]; then
            echo "  Found existing analysis directory (likely from an unfinished run); removing: $abs_analysis_dir"
            rm -rf "$abs_analysis_dir"
        fi

        # -q = Quiet
        # -j = Junk paths (extract file to root of destination, not its internal path)
        # -d = Destination directory
        #unzip -qj "$abs_zip_file" "$CHECK_FILE" -d "$abs_analysis_dir"
	echo "$abs_analysis_dir"/"$CHECK_FILE"
	unzip -oqj "$abs_zip_file" "*/$CHECK_FILE" -d $abs_analysis_dir/
	#unzip -j 38_C-TSa1_0_1.zip 38_C-TSa1_0_1/finished_irc.trj -d 38_C-TSa1_0_1/
        
        echo "  Running analysis in: $abs_analysis_dir"
        
        # Change into the analysis directory to run the script
        # This is CRITICAL because the python script expects
        # 'finished_irc.trj' to be in its current directory.
        cd "$abs_analysis_dir"

        # --- 🤖 EXECUTING YOUR ANALYSIS 🤖 ---
        python "$PYTHON_SCRIPT_PATH" \
          --trj "$CHECK_FILE" \
          --charge "$ARG_CHARGE" \
          --mult "$ARG_MULT" \
          --uhf-max "$ARG_UHF_MAX" \
          --workdir "$ARG_WORKDIR" \
          --nprocs "$ARG_NPROCS"
        
        echo "  Analysis complete."
        # --- END OF ANALYSIS ---
        

        # 3. Append new/modified files back to the zip and clean up
        echo "  Updating zip file with new contents (e.g., $ARG_WORKDIR/)..."
        
        # We are already in $abs_analysis_dir
        
        # Update the original zip file
        # -u = Update (adds new files)
        # -r = Recursive (to include the new $ARG_WORKDIR folder)
        # -q = Quiet
        # This will add the entire 'xTB-scan' directory to the zip.
	cd "$abs_analysis_dir"/../
        zip -urq "$abs_zip_file" $(basename "$abs_analysis_dir")/*
        
        # Go back to the search directory to continue the loop
        cd "$ABS_SEARCH_DIR"
        
        # Remove the temporary analysis folder
        echo "  Cleaning up $abs_analysis_dir..."
        rm -r "$abs_analysis_dir"
        
    else
        # This block runs if the check file was not found
        echo "  Skipping: $CHECK_FILE not found."
    fi
    
    echo "---" # Separator
done

# Go back to the directory the user started from
cd "$START_DIR"

echo "🎉 All processing complete."

