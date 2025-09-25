# LDSC

#!/bin/bash

# Path to the directory containing GWAS .txt files
## Expects GWAS .txt with columns: SNP, A1, A2, beta, pval, gwas_ntotal
GWAS_DIR="/path/to/directory/containing/GWAS/files"
# Path to the cloned LDSC GitHub repository
LDSC_DIR="/path/to/the/cloned/LDSC/github/repo"

# Step 1: Navigate to the cloned LDSC directory
cd "$LDSC_DIR" || { echo "Error: Cannot navigate to LDSC directory"; exit 1; }

# Step 2: Check if the LDSC environment exists
if ! conda env list | grep -q "^ldsc"; then
echo "Creating LDSC environment..."
conda env create --file environment.yml
else
 echo "LDSC environment already exists. Skipping creation."
fi

# Step 3: Activate the LDSC environment
conda activate ldsc || { echo "Error: Failed to activate LDSC environment"; exit 1; }

# Generate timestamp for folder names
TIMESTAMP=$(date +"%Y_%d_%m_%H_%M_%S")
MUNGED_DIR="${GWAS_DIR}/ldsc_files_${TIMESTAMP}"
RESULTS_DIR="${GWAS_DIR}/LDSC_results_${TIMESTAMP}"

# Step 4: Create directories for munged files and results
mkdir -p "$MUNGED_DIR"
mkdir -p "$RESULTS_DIR"

# Step 5: Iterate through all .txt files in the GWAS directory
for GWAS_FILE in "$GWAS_DIR"/*.txt; do
# Extract the base filename without the directory
BASENAME=$(basename "$GWAS_FILE" .txt)

# Extract `gwas_ntotal` value from the file (assuming it's in the first row of the `gwas_ntotal` column)
GWAS_NTOTAL=$(awk -F' ' '
NR == 1 { for (i=1; i<=NF; i++) if ($i == "gwas_ntotal") col=i } 
NR == 2 { if (col) print $col; exit }
' "$GWAS_FILE")


# Validate the extracted value to ensure it is numeric
if ! [[ "$GWAS_NTOTAL" =~ ^[0-9]+$ ]]; then
echo "Warning: Invalid or missing gwas_ntotal for $GWAS_FILE. Skipping file."
continue
fi

# Run the munge_sumstats.py script
./munge_sumstats.py \
--sumstats "$GWAS_FILE" \
--snp SNP \
--a1 A1 \
--a2 A2 \
--signed-sumstats beta,0 \
--p pval \
--N "$GWAS_NTOTAL" \
--chunksize 500000 \
--out "${MUNGED_DIR}/${BASENAME}_munged" \
--merge-alleles w_hm3.snplist

echo "Finished processing $GWAS_FILE"
done

echo "All files processed. Munged files are in $MUNGED_DIR."

# Step 6: Perform rg analysis for all combinations of .sumstats.gz files in MUNGED_DIR
MUNGED_FILES=(${MUNGED_DIR}/*_munged.sumstats.gz)

for ((i=0; i<${#MUNGED_FILES[@]}; i++)); do
 # Current file
 FILE_1=${MUNGED_FILES[$i]}
 
 # Create a comma-separated list of all other files
 FILE_OTHERS=$(printf ",%s" "${MUNGED_FILES[@]}")
 FILE_OTHERS=${FILE_OTHERS:1} # Remove leading comma
 
 # Run rg analysis
 ./ldsc.py \
 --rg "${FILE_1},${FILE_OTHERS}" \
 --ref-ld-chr eur_w_ld_chr/ \
 --w-ld-chr eur_w_ld_chr/ \
 --out "${RESULTS_DIR}/$(basename ${FILE_1%_munged.sumstats.gz})_vs_all"
 done
 
 # Step 7: Extract and merge Summary of Genetic Correlation Results tables
 SUMMARY_FILE="${RESULTS_DIR}/00_summary_of_all_ldsc_rg_analyses_${TIMESTAMP}.txt"
 
 # Initialize the summary file with headers
 echo -e "p1\tp1_userfacinglabel\tp2\tp2_userfacinglabel\trg\tse\tz\tp\th2_obs\th2_obs_se\th2_int\th2_int_se\tgcov_int\tgcov_int_se" > "$SUMMARY_FILE"
 
 # Iterate over all log files in RESULTS_DIR
 for LOG_FILE in "${RESULTS_DIR}"/*_vs_all.log; do
 # Extract the "Summary of Genetic Correlation Results" table
 awk '/^Summary of Genetic Correlation Results/,/^Analysis finished/' "$LOG_FILE" | \
 awk 'NR > 2 && !/^Analysis finished/ && NF > 0 { print $0 }' | \
 while read -r LINE; do
 # Extract p1 and p2 from the line
 P1=$(echo "$LINE" | awk '{print $1}')
 P2=$(echo "$LINE" | awk '{print $2}')
 
 # Generate user-facing labels for p1 and p2
 P1_USER_LABEL=$(echo "$P1" | sed -E 's|.*/([^/]+)_build.*|\1|')
 P2_USER_LABEL=$(echo "$P2" | sed -E 's|.*/([^/]+)_build.*|\1|')
 
 # Append the processed line to the summary file with user-facing labels
 echo -e "$P1\t$P1_USER_LABEL\t$P2\t$P2_USER_LABEL\t$(echo "$LINE" | awk '{$1=$2=""; print $0}' | sed 's/^ *//')" >> "$SUMMARY_FILE"
 done
 done
 
 # Remove any trailing blank lines in the final file
 sed -i '/^$/d' "$SUMMARY_FILE"
 
 echo "Summary of Genetic Correlation Results from all .log files, including user-facing labels, has been saved to $SUMMARY_FILE."
 