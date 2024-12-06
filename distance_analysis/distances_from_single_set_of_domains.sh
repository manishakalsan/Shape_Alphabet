#!/bin/bash

# Check if correct number of arguments are provided
if [ "$#" -lt 5 ]; then
    echo "Usage: $0 <tf_list_file> <bed_folder> <bed_file> <n_hits> <output_folder>"
    exit 1
fi

# Arguments
TF_LIST_FILE=$1    # Text file containing the list of TFs (one TF per line)
TF_BED_FOLDER=$2   # Folder containing BED files for each TF binding site
BED_FILE=$3        # BED file to compare against
N_HITS=$4	   # Number of closest hits to be reported for each region
OUTPUT_FOLDER=$5   # Output folder for intermediate results and final plots

# Create output folder if it doesn't exist
mkdir -p "$OUTPUT_FOLDER"
echo "Output folder: $OUTPUT_FOLDER"

# Sort the main BED file for comparison
BED_FILE_SORTED="$OUTPUT_FOLDER/$(basename "$BED_FILE")_sorted.bed"
sort -k1,1 -k2,2n "$BED_FILE" > "$BED_FILE_SORTED"

# Check if sorting succeeded
if [ ! -f "$BED_FILE_SORTED" ]; then
    echo "Error: sorting of $BED_FILE failed."
    exit 1
fi

# Initialize CSV list for later Python visualization
CSV_LIST=""

# Read the TF names from the TF_LIST_FILE and process each one
while IFS= read -r TF_NAME; do
    TF_BED="$TF_BED_FOLDER/${TF_NAME}_tfbs.bed"
    
    # Check if the BED file exists for the given TF
    if [ ! -f "$TF_BED" ]; then
        echo "Warning: BED file for $TF_NAME not found in $TF_BED_FOLDER, skipping $TF_NAME..."
        continue
    fi
    
    # Output files for distances
    DIST_FILE="$OUTPUT_FOLDER/${TF_NAME}_distances.txt"
    CSV_OUTPUT="$OUTPUT_FOLDER/${TF_NAME}_distances.csv"
        
    # Sort the BED file for the TF
    TF_BED_SORTED="$OUTPUT_FOLDER/${TF_NAME}_sorted.bed"
    sort -k1,1 -k2,2n "$TF_BED" > "$TF_BED_SORTED"

    # Check if sorting succeeded
    if [ ! -f "$TF_BED_SORTED" ]; then
        echo "Error: sorting of $TF_BED failed for $TF_NAME."
        continue
    fi
    
    # Calculate distances for the first region set from the BEDPE file
    echo "Calculating distances for $TF_NAME ..."
    bedtools closest -a "$TF_BED_SORTED" -b "$BED_FILE_SORTED" -k "$N_HITS" -d > "$DIST_FILE"

    # Check if the distance calculation succeeded
    if [ $? -ne 0 ]; then
        echo "Error: Distance calculation failed for $TF_NAME."
        continue
    fi
    
    # Extract distances and save them in CSV format for easier plotting
    echo "TF,RegionSet,Distance" > "$CSV_OUTPUT"
    cut -f7 "$DIST_FILE" | awk '{print "'$TF_NAME',domains,"$1}' >> "$CSV_OUTPUT"
    
    # Append CSV path to the list for Python visualization
    CSV_LIST="$CSV_LIST $CSV_OUTPUT"
    
done < "$TF_LIST_FILE"

echo "Processing complete. CSV files created: $CSV_LIST"



