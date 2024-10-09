#!/bin/bash

tfnames=("NR2C2" "CTCF" "YBX1" "ATF2" "CEBPB" "RAD21" "IKZF1" "ZNF207" "NR2F2" "NR2F6" "NRF1")

for TFname in "${tfnames[@]}"; do
    # Define the corresponding fasta file
    fasta_file="${TFname}_tfbs.fasta"

    # Check if the file exists
    if [[ -f "$fasta_file" ]]; then
        # Define the output directory
        output_dir="${TFname}_shape_meme"

        # Run meme command with the specified parameters
        meme -alph shape_alphabet -nmotifs 3 -minw 5 -maxw 15 -mod zoops -maxsize 100000000 -oc "$output_dir" "$fasta_file"

        echo "Processed $fasta_file and output to $output_dir"
    else
        echo "File $fasta_file not found!"
    fi
done


for TFname in "${tfnames[@]}"; do
    # Define the corresponding fasta file
    fasta_file="shuffled_${TFname}.fasta"

    # Check if the file exists
    if [[ -f "$fasta_file" ]]; then
        # Define the output directory
        output_dir="shuffled_${TFname}_shape_meme"

        # Run meme command with the specified parameters
        meme -alph shape_alphabet -nmotifs 3 -minw 5 -maxw 15 -mod zoops -maxsize 100000000 -oc "$output_dir" "$fasta_file"

        echo "Processed $fasta_file and output to $output_dir"
    else
        echo "File $fasta_file not found!"
    fi
done
