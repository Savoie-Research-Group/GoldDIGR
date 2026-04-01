#!/bin/bash

# --- Configuration ---
FINISHED_LIST="finished.txt"
# The name of the final archive file (it will be created *inside* the base directory)
TAR_OUTPUT_FILE="finished_results.tar"

echo "--- Zip Archiver and Cleaner ---"

# --- 1. Get and Validate User Input ---

# Ask the user for the base directory
read -p "Enter the absolute base directory for 'tar' (e.g., /scratch/bell/li1724/SI-Downloads/SI_Agent/New-ACS/): " TAR_BASE_DIR

# Validate the input
if [ -z "$TAR_BASE_DIR" ]; then
    echo "Error: No directory provided. Exiting."
    exit 1
fi

# Remove trailing slash (if any) for path consistency
TAR_BASE_DIR=${TAR_BASE_DIR%/}

if [ ! -d "$TAR_BASE_DIR" ]; then
    echo "Error: Base directory not found at '$TAR_BASE_DIR'. Exiting."
    exit 1
fi

if [ ! -f "$FINISHED_LIST" ]; then
    echo "Error: Input file '$FINISHED_LIST' not found. Exiting."
    exit 1
fi

# --- 2. Prepare File List for Tar ---

# Create a secure temporary file to store the list of *relative* file paths
# We will pass this file to 'tar'
TEMP_FILE_LIST=$(mktemp)

# Set a trap: This ensures the temporary file is *always* deleted when the script exits,
# whether it finishes successfully, fails, or is cancelled (e.g., Ctrl+C).
trap 'rm -f "$TEMP_FILE_LIST"' EXIT

echo "Scanning '$FINISHED_LIST' and preparing file list..."

file_count=0
while IFS= read -r BASE_PATH || [[ -n "$BASE_PATH" ]]; do
    if [ -z "$BASE_PATH" ]; then
        continue
    fi

    # Construct the full, absolute path to the .zip file
    ZIP_FILE_ABS_PATH="${BASE_PATH}.zip"

    # Check if the zip file actually exists
    if [ -f "$ZIP_FILE_ABS_PATH" ]; then
        # Calculate the file's path *relative* to the user's base directory
        # 'realpath --relative-to=FROM TO' is perfect for this.
        # We redirect stderr (2>/dev/null) to suppress warnings if the file isn't in the dir.
        RELATIVE_PATH=$(realpath --relative-to="$TAR_BASE_DIR" "$ZIP_FILE_ABS_PATH" 2>/dev/null)

        # Check if realpath succeeded (exit code $? is 0)
        if [ $? -eq 0 ]; then
            # Add the relative path to our temporary list file
            echo "$RELATIVE_PATH" >> "$TEMP_FILE_LIST"
            ((file_count++))
        else
            # This happens if the zip file is *not* inside the TAR_BASE_DIR
            echo "Warning: File '$ZIP_FILE_ABS_PATH' is not inside the base directory '$TAR_BASE_DIR'. Skipping."
        fi
    else
        echo "Warning: Zip file not found, skipping: $ZIP_FILE_ABS_PATH"
    fi
done < "$FINISHED_LIST"

# --- 3. Check if any files were found ---
if [ $file_count -eq 0 ]; then
    echo "No valid .zip files were found to archive. Exiting."
    # The 'trap' will clean up $TEMP_FILE_LIST
    exit 0
fi

echo "Found $file_count .zip files to archive."
echo "---"

# --- 4. Create Tar and Remove Files ---

# Get the absolute path for the *output* file.
# This places it inside the user's provided base directory.
TAR_OUTPUT_ABS_PATH="${TAR_BASE_DIR}/${TAR_OUTPUT_FILE}"

echo "Creating archive at: $TAR_OUTPUT_ABS_PATH"
echo "This will also REMOVE the $file_count original .zip files after they are archived."

# Run the tar command
# 1. (cd "$TAR_BASE_DIR" && ...): Change to the base directory first.
# 2. tar -cvf "$TAR_OUTPUT_ABS_PATH": Create a new, verbose archive at the specified *absolute* path.
# 3. -T "$TEMP_FILE_LIST": Read the list of *relative* file paths from our temp file.
# 4. --remove-files: Securely removes files *after* they are successfully added to the archive.
(cd "$TAR_BASE_DIR" && tar -cvf "$TAR_OUTPUT_ABS_PATH" -T "$TEMP_FILE_LIST" --remove-files)

# --- 5. Final Report ---
if [ $? -eq 0 ]; then
    echo "---"
    echo "✅ Success! Archive created at $TAR_OUTPUT_ABS_PATH"
    echo "All $file_count .zip files have been safely archived and removed."
else
    echo "---"
    echo "❌ Error: 'tar' command failed."
    echo "The .zip files have NOT been removed."
    echo "Please check the output above for errors."
fi

# The 'trap' command will automatically clean up $TEMP_FILE_LIST upon exit
echo "Script finished."

