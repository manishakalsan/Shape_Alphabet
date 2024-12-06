#!/bin/bash

# print each command being executed
set -x

# Check if correct number of arguments are provided
if [ "$#" -lt 3 ]; then
    echo "Usage: $0 <tf_list_file> <bed_folder> <bedpe_file> <n_hits> <output_folder>"
    exit 1
fi

# Arguments
TF_LIST_FILE=$1    # Text file containing the list of TFs (one TF per line)
TF_BED_FOLDER=$2   # Folder containing BED files for each TF binding site
BEDPE_FILE=$3      # BEDPE file
OUTPUT_FOLDER=$5   # Output folder for intermediate results and final plots
N_HITS=$4	   # Number of closest hits to be reported for each region

# Create output folder if it doesn't exist
mkdir -p "$OUTPUT_FOLDER"
echo "Output folder: $OUTPUT_FOLDER"

 # Change the format of the files to add chr prefix to chr columns
 awk 'BEGIN {OFS="\t"} NR==1 {print $0} NR>1 {$1="chr"$1; $4="chr"$4; print $0}' $BEDPE_FILE > "$BEDPE_FILE_formatted"

# Temporary file names
TEMP_FIRST_BEDPE="$OUTPUT_FOLDER/bedpe_first_region.bed"
TEMP_SECOND_BEDPE="$OUTPUT_FOLDER/bedpe_second_region.bed"
    
# Extract first and second regions from the BEDPE file
cut -f1-3 "$BEDPE_FILE_formatted" > "$TEMP_FIRST_BEDPE"
cut -f4-6 "$BEDPE_FILE_formatted" > "$TEMP_SECOND_BEDPE"
    
# Sort BEDPE file's first and second region sets before processing
sort -k1,1 -k2,2n "$TEMP_FIRST_BEDPE" > "$TEMP_FIRST_BEDPE_sorted.bed"
sort -k1,1 -k2,2n "$TEMP_SECOND_BEDPE" > "$TEMP_SECOND_BEDPE_sorted.bed"

# Read the TF names from the TF_LIST_FILE and process each one
while IFS= read -r TF_NAME; do
    TF_BED="$TF_BED_FOLDER/${TF_NAME}_tfbs.bed"

    # Check if the BED file exists for the given TF
    if [ ! -f "$TF_BED" ]; then
        echo "Warning: BED file for $TF_NAME not found in $TF_BED_FOLDER, skipping..."
        continue
    fi

    # Output files for distances
    DIST_FIRST="$OUTPUT_FOLDER/${TF_NAME}_first_distances.txt"
    DIST_SECOND="$OUTPUT_FOLDER/${TF_NAME}_second_distances.txt"
    CSV_OUTPUT="$OUTPUT_FOLDER/${TF_NAME}_distances.csv"
    
    # Sort the BED file for the TF
    TF_BED_SORTED="$OUTPUT_FOLDER/${TF_NAME}_sorted.bed"
    sort -k1,1 -k2,2n "$TF_BED" > "$TF_BED_SORTED"

    # Check if sorted files exist
    if [ ! -f "$TEMP_FIRST_BEDPE_sorted.bed" ] || [ ! -f "$TEMP_SECOND_BEDPE_sorted.bed" ]; then
        echo "Error: BEDPE sorted files not found, skipping $TF_NAME."
        continue
    fi

    # Calculate distances for the first region set from the BEDPE file
    echo "Calculating distances for $TF_NAME (first BEDPE region)..."
    bedtools closest -a "$TF_BED_SORTED" -b "$TEMP_FIRST_BEDPE_sorted.bed" -k "$N_HITS" -d > "$DIST_FIRST"

    # Check if the distance calculation succeeded
    if [ $? -ne 0 ]; then
        echo "Error: Distance calculation failed for $TF_NAME (first BEDPE region)."
        continue
    fi
    
    # Calculate distances for the second region set from the BEDPE file
    echo "Calculating distances for $TF_NAME (second BEDPE region)..."
    bedtools closest -a "$TF_BED_SORTED" -b "$TEMP_SECOND_BEDPE_sorted.bed" -k "$N_HITS" -d > "$DIST_SECOND"
    
    # Check if the distance calculation succeeded
    if [ $? -ne 0 ]; then
        echo "Error: Distance calculation failed for $TF_NAME (second BEDPE region)."
        continue
    fi

    # Extract distances and save them in CSV format for easier plotting
    echo "TF,RegionSet,Distance" > "$CSV_OUTPUT"
    cut -f7 "$DIST_FIRST" | awk '{print "'$TF_NAME',first,"$1}' >> "$CSV_OUTPUT"
    cut -f7 "$DIST_SECOND" | awk '{print "'$TF_NAME',second,"$1}' >> "$CSV_OUTPUT"
    
    # Append CSV path to the list for Python visualization
    CSV_LIST="$CSV_LIST $CSV_OUTPUT"
    
done 

