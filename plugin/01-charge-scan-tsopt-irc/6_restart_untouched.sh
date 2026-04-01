#!/bin/bash

# activate a correct anaconda environment
module load anaconda/2022.10-py39
#source activate /scratch/negishi/li1724/100624_copy_yarp/
conda activate copy-classy-yarp
#conda activate IOP
module load intel-mkl
module load openmpi
# --- Configuration ---
# File containing the absolute paths to the *full job directories*
# (e.g., .../19_-1_2), one per line.
UNTOUCHED_LIST_FILE="untouched.txt"
# SLURM account (replace if 'standby' is not correct for you)
SLURM_ACCOUNT="standby"
# A unique name for this batch of jobs to track them in the queue
BATCHNAME="AAAS"
# The maximum number of jobs from this batch to have in the queue at one time
JOB_LIMIT=1000
# How many seconds to wait before re-checking the queue when the limit is reached
CHECK_EVERY=30

# --- Script Body ---
if [ ! -f "$UNTOUCHED_LIST_FILE" ]; then
    echo "Error: Untouched jobs file not found at '$UNTOUCHED_LIST_FILE'"
    exit 1
fi

echo "Starting restart workflow from '$UNTOUCHED_LIST_FILE'..."
echo "Job limit: ${JOB_LIMIT} for batch '${BATCHNAME}'."

while IFS= read -r FULL_DIR_PATH || [[ -n "$FULL_DIR_PATH" ]]; do
    if [ -z "$FULL_DIR_PATH" ]; then
        continue
    fi

    echo "Processing job: $FULL_DIR_PATH"

    # --- NEW: Parse Charge, Multiplicity, and Base Path ---
    # Example FULL_DIR_PATH: /path/to/block/19_-1_2
    
    # Get multiplicity (e.g., "2")
    multiplicity=${FULL_DIR_PATH##*_}
    
    # Get path without multiplicity (e.g., "/path/to/block/19_-1")
    TEMP_PATH=${FULL_DIR_PATH%_${multiplicity}}
    
    # Get charge (e.g., "-1")
    charge=${TEMP_PATH##*_}
    
    # Get the original base path (e.g., "/path/to/block/19")
    PATH_NO_SLASH=${TEMP_PATH%_${charge}}
    
    # Get the path to the original XYZ file (e.g., "/path/to/block/19.xyz")
    XYZ_FILE_PATH="${PATH_NO_SLASH}.xyz"

    # --- Validation ---
    if [ -z "$charge" ] || [ -z "$multiplicity" ]; then
        echo "  -> ERROR: Could not parse charge/multiplicity from '$FULL_DIR_PATH'. Skipping."
        continue
    fi
    
    if [ ! -f "$XYZ_FILE_PATH" ]; then
        echo "  -> ERROR: Original input file not found: '$XYZ_FILE_PATH'. Skipping."
        continue
    fi
    # --- End New Parsing ---


    # --- Job Name Generation (from original script) ---
    # Extract the index from the path using a regex.
    index=$(echo "$PATH_NO_SLASH" | grep -oP '(?:/All-XYZ|/all-pdf-xyz)/\K[^/]+')

    if [ -z "$index" ]; then
        echo "  -> WARNING: Could not extract index for path: $PATH_NO_SLASH. Using 'NA' for index."
        index="NA"
    fi

    # Extract the block name (parent directory of the .xyz file)
    block=$(basename "$(dirname "$XYZ_FILE_PATH")")
    BASENAME=$(basename "$XYZ_FILE_PATH" .xyz)
    
    # Construct the new, more specific job name.
    JOB_NAME="${BATCHNAME}-${index}-${block}-${BASENAME}-${charge}"
    # --- End Job Name Generation ---


    # --- Job Status Checks (from original script) ---
    if [ -f "${FULL_DIR_PATH}.zip" ]; then
        echo "  -> Skipping ${FULL_DIR_PATH}: Result zip file already exists."
        continue
    fi
    
    if squeue -u "$USER" --noheader -o "%.200j" | grep -q "${JOB_NAME}$"; then
        echo "  -> Skipping ${FULL_DIR_PATH}: Job is already in the SLURM queue."
        continue
    fi

    if [ -d "${FULL_DIR_PATH}" ]; then
        echo "  -> Found stale directory for ${FULL_DIR_PATH}. Cleaning up and restarting."
        rm -rf "${FULL_DIR_PATH}"
    fi

    # --- Job Limit and Submission Logic (from original script) ---
    while true; do
        current_jobs=$(squeue -u "$USER" -h -o "%j" | grep -F "$BATCHNAME" | wc -l)

        if (( current_jobs < JOB_LIMIT )); then
            echo "  --> Preparing and submitting job for ${FULL_DIR_PATH}. Current jobs: ${current_jobs}/${JOB_LIMIT}."
            
            mkdir -p "$FULL_DIR_PATH"
            cp "$XYZ_FILE_PATH" "${FULL_DIR_PATH}/input.xyz"
            
            # Use the parsed charge and multiplicity to fill templates
            sed -e "s/CHARGE_PLACEHOLDER/${charge}/g" -e "s/MULT_PLACEHOLDER/${multiplicity}/g" TSOPT.yaml.template > "${FULL_DIR_PATH}/TSOPT.yaml"
            sed -e "s/CHARGE_PLACEHOLDER/${charge}/g" -e "s/MULT_PLACEHOLDER/${multiplicity}/g" TSOPT-IRC.yaml.template > "${FULL_DIR_PATH}/TSOPT-IRC.yaml"
            
            template_content=$(<job_tsopt.submit.template)
            modified_content="${template_content//XXX/${BASENAME}}"
            modified_content="${modified_content//CHARGE_PLACEHOLDER/${charge}}"
            modified_content="${modified_content//JOBNAME_PLACEHOLDER/${JOB_NAME}}"
            echo "$modified_content" > "${FULL_DIR_PATH}/job_tsopt.submit"

            (cd "$FULL_DIR_PATH" && sbatch job_tsopt.submit)
            
            break # Job submitted, exit the 'while true' loop
        else
            printf '[%s] Job limit reached. %d jobs in queue (limit %d). Waiting %ds...\n' \
                   "$(date '+%H:%M:%S')" "$current_jobs" "$JOB_LIMIT" "$CHECK_EVERY"
            sleep "$CHECK_EVERY"
        fi
    done
    
done < "$UNTOUCHED_LIST_FILE"

echo "Workflow finished. All untouched jobs have been submitted."

