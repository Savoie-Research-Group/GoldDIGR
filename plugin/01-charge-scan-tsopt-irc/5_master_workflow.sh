#!/bin/bash

# activate a correct anaconda environment
module load anaconda/2022.10-py39
#source activate /scratch/negishi/li1724/100624_copy_yarp/
conda activate copy-classy-yarp
#conda activate IOP
module load intel-mkl
module load openmpi

# --- Configuration ---
# File containing the absolute paths to the XYZ files, one per line.
XYZ_LIST_FILE="AAAS_folders.txt"
# SLURM account (replace if 'standby' is not correct for you)
SLURM_ACCOUNT="standby"
# A unique name for this batch of jobs to track them in the queue
BATCHNAME="PAAAS"
# The maximum number of jobs from this batch to have in the queue at one time
JOB_LIMIT=100
# How many seconds to wait before re-checking the queue when the limit is reached
CHECK_EVERY=30

# --- Script Body ---
if [ ! -f "$XYZ_LIST_FILE" ]; then
    echo "Error: XYZ list file not found at '$XYZ_LIST_FILE'"
    exit 1
fi

echo "Starting workflow with a job limit of ${JOB_LIMIT} for batch '${BATCHNAME}'."

while IFS= read -r DIR_PATH || [[ -n "$DIR_PATH" ]]; do
    if [ -z "$DIR_PATH" ]; then
        continue
    fi

    PATH_NO_SLASH=${DIR_PATH%/}
    XYZ_FILE_PATH="${PATH_NO_SLASH}.xyz"

    # --- MODIFICATION START ---
    # Extract the index from the path using a regex.
    # It looks for the string immediately following '/All-XYZ/' or '/all-pdf-xyz/'.
    index=$(echo "$DIR_PATH" | grep -oP '(?:/All-XYZ|/all-pdf-xyz)/\K[^/]+')

    # If the index isn't found, use a default placeholder to avoid an empty name.
    if [ -z "$index" ]; then
        echo "  -> WARNING: Could not extract index for path: $DIR_PATH. Using 'NA' for index."
        index="NA"
    fi

    # Extract the block name by getting the basename of the parent directory
    # e.g., /path/to/block_name/file.xyz -> block_name
    block=$(basename "$(dirname "$XYZ_PATH")")

    # --- MODIFICATION END ---


    echo "Processing molecule: $XYZ_FILE_PATH (Index: $index)"
    BASENAME=$(basename "$XYZ_FILE_PATH" .xyz)

    for charge in -1 0 1; do
        multiplicity=$(python calculate_mult.py "$XYZ_FILE_PATH" "$charge")

        if [ $? -ne 0 ]; then
            echo "Error calculating multiplicity for $BASENAME with charge $charge. Skipping."
            continue
        fi

        FULL_DIR_PATH="${PATH_NO_SLASH}_${charge}_${multiplicity}"
        
        # --- MODIFICATION ---
        # Construct the new, more specific job name.
        JOB_NAME="${BATCHNAME}-${index}-${block}-${BASENAME}-${charge}"
        # ------------------

        if [ -f "${FULL_DIR_PATH}.zip" ]; then
            echo "  -> Skipping ${FULL_DIR_PATH}: Result zip file already exists."
            continue
        fi
	if squeue -u "$USER" --noheader -o "%.200j" | grep -q "${JOB_NAME}$"; then
	#if squeue -u "$USER" --noheader -o "%.200j" | grep -q -F -x "$JOB_NAME"; then
        #if squeue -u "$USER" -n "$JOB_NAME" --noheader | grep -q "$JOB_NAME"; then
            echo "  -> Skipping ${FULL_DIR_PATH}: Job is already in the SLURM queue."
	    #echo "FOUND $JOB_NAME"; exit
            continue
        fi
	#echo "NOT FOUND $JOB_NAME"
	#exit

        if [ -d "${FULL_DIR_PATH}" ]; then
            echo "  -> Found stale directory for ${FULL_DIR_PATH}. Cleaning up and restarting."
            rm -rf "${FULL_DIR_PATH}"
        fi

        while true; do
            current_jobs=$(squeue -u "$USER" -h -o "%j" | grep -F "$BATCHNAME" | wc -l)

            if (( current_jobs < JOB_LIMIT )); then
                echo "  --> Preparing and submitting job for ${FULL_DIR_PATH}. Current jobs: ${current_jobs}/${JOB_LIMIT}."
                
                mkdir -p "$FULL_DIR_PATH"
                cp "$XYZ_FILE_PATH" "${FULL_DIR_PATH}/input.xyz"
                sed -e "s/CHARGE_PLACEHOLDER/${charge}/g" -e "s/MULT_PLACEHOLDER/${multiplicity}/g" TSOPT.yaml.template > "${FULL_DIR_PATH}/TSOPT.yaml"
                sed -e "s/CHARGE_PLACEHOLDER/${charge}/g" -e "s/MULT_PLACEHOLDER/${multiplicity}/g" TSOPT-IRC.yaml.template > "${FULL_DIR_PATH}/TSOPT-IRC.yaml"
                
                template_content=$(<job_tsopt.submit.template)
                modified_content="${template_content//XXX/${BASENAME}}"
                modified_content="${modified_content//CHARGE_PLACEHOLDER/${charge}}"
                # The placeholder in the template is now replaced with the new dynamic JOB_NAME
                modified_content="${modified_content//JOBNAME_PLACEHOLDER/${JOB_NAME}}"
                echo "$modified_content" > "${FULL_DIR_PATH}/job_tsopt.submit"

                (cd "$FULL_DIR_PATH" && sbatch job_tsopt.submit)
                
                break
            else
                printf '[%s] Job limit reached. %d jobs in queue (limit %d). Waiting %ds...\n' \
                       "$(date '+%H:%M:%S')" "$current_jobs" "$JOB_LIMIT" "$CHECK_EVERY"
                sleep "$CHECK_EVERY"
            fi
        done
    done
done < "$XYZ_LIST_FILE"

echo "Workflow finished. All jobs have been submitted."

