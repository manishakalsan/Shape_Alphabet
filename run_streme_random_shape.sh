#!/bin/bash

tfnames=("NR2C2" "CTCF" "YBX1" "ATF2" "CEBPB" "RAD21" "IKZF1" "ZNF207" "NR2F2" "NR2F6" "NRF1")

for TFname in "${tfnames[@]}"; do
    # Define the corresponding fasta file
    fasta_file="${TFname}_tfbs.fasta"
    shuffled_file="shuffled_${TFname}.fasta"

    # Check if the file exists
    if [[ -f "$fasta_file" ]]; then
        # Define the output directory
        output_dir="${TFname}_streme_random_shape"

        # Run meme command with the specified parameters
        streme -alph random_shape_alphabet -oc "$output_dir" -minw 5 -maxw 15 -nmotifs 3 -p "$fasta_file" -n "$shuffled_file"

        echo "Processed $fasta_file and output to $output_dir"
    else
        echo "File $fasta_file not found!"
    fi
done


for TFname in "${tfnames[@]}"; do
    # Define the corresponding fasta file
    fasta_file="${TFname}_tfbs.fasta"
    shuffled_file="shuffled_${TFname}.fasta"

    # Check if the file exists
    if [[ -f "$fasta_file" ]]; then
        # Define the output directory
        output_dir="${TFname}_streme_random_shape_cd"

        # Run meme command with the specified parameters
        streme -alph random_shape_alphabet -oc "$output_dir" -minw 5 -maxw 15 -nmotifs 3 -objfun cd -p "$fasta_file" -n "$shuffled_file"

        echo "Processed $fasta_file and output to $output_dir"
    else
        echo "File $fasta_file not found!"
    fi
done


