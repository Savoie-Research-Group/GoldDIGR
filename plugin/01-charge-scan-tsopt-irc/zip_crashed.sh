#!/bin/bash

# --- Configuration ---
CRASHED_LIST="crashed.txt"
echo "--- Crashed Job Archiver ---"

# --- 1. Validate Input ---
if [ ! -f "$CRASHED_LIST" ]; then
    echo "Error: Input file '$CRASHED_LIST' not found. Exiting."
    exit 1
fi

# --- 2. Initialize Counters ---
count_total=0
count_success=0
count_skipped_zip=0
count_skipped_dir=0
count_failed=0

echo "Scanning '$CRASHED_LIST' to find and zip crashed directories..."

# --- 3. Process File List ---
while IFS= read -r DIR_PATH || [[ -n "$DIR_PATH" ]]; do
    if [ -z "$DIR_PATH" ]; then
        continue
    fi
    
    ((count_total++))

    # Construct the full, absolute path to the .zip file we want to create
    ZIP_OUTPUT_PATH="${DIR_PATH}-crashed.zip"

    # --- Check 1: Does the target directory exist? ---
    if [ ! -d "$DIR_PATH" ]; then
        echo "  -> Skipping: Directory not found, may have been processed: $DIR_PATH"
        ((count_skipped_dir++))
        continue
    fi
    
    # --- Check 2: Does the output zip file *already* exist? ---
    if [ -f "$ZIP_OUTPUT_PATH" ]; then
        echo "  -> Skipping: Crashed zip file already exists: $ZIP_OUTPUT_PATH"
        ((count_skipped_zip++))
        continue
    fi
    
    echo "  -> Processing: $DIR_PATH"

    # Get the parent directory (e.g., /path/to/parent)
    PARENT_DIR=$(dirname "$DIR_PATH")
    # Get the name of the directory to zip (e.g., 28_-1_2)
    DIR_NAME=$(basename "$DIR_PATH")
    
    # Construct the relative name for the zip file (e.g., 28_-1_2-crashed.zip)
    ZIP_NAME="${DIR_NAME}-crashed.zip"
    
    # --- 4. Zip and Remove ---
    # We (cd) to the parent directory first. This ensures that:
    # 1. The zip file is created in the correct location (the parent dir).
    # 2. The paths *inside* the zip file are relative (e.g., "28_-1_2/output.log")
    #    and not absolute (e.g., "/scratch/bell/li1724/.../28_-1_2/output.log").
    # We use 'zip -q -r' for 'quiet' and 'recursive'.
    
    echo "     Zipping to: $ZIP_OUTPUT_PATH"
    (cd "$PARENT_DIR" && zip -q -r "$ZIP_NAME" "$DIR_NAME")
    
    # Check if the zip command was successful ($? -eq 0)
    if [ $? -eq 0 ]; then
        echo "     Success. Removing original directory: $DIR_PATH"
        rm -rf "$DIR_PATH"
        ((count_success++))
    else
        echo "     ❌ ERROR: Zip command failed for $DIR_PATH."
        echo "     Original directory has NOT been removed."
        ((count_failed++))
    fi

done < "$CRASHED_LIST"

# --- 5. Final Report ---
echo "---"
echo "✅ Crashed job archival finished."
echo "Summary:"
echo "  Processed paths: $count_total"
echo "  Successfully zipped & removed: $count_success"
echo "  Failed to zip: $count_failed"
echo "  Skipped (zip already exists): $count_skipped_zip"
echo "  Skipped (directory not found): $count_skipped_dir"
