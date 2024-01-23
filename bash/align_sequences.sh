#!/bin/bash
set -x  # This will print each command and its arguments to the terminal

serotypes=("Dengue_1" "Dengue_2" "Dengue_3" "Dengue_4")
reference_serotype=("1" "2" "3" "4")

# Loop over the serotypes
for i in "${!serotypes[@]}"; do
    # Define the input sequence file and reference sequence file based on the serotype
    sequences="results/Unaligned_${serotypes[i]}.fasta"
    reference="reference_genomes/reference_denv${reference_serotype[i]}.fasta"
    genemap="genemap/genemap_denv${reference_serotype[i]}.gff"
    output_dir="results/Aligned_${serotypes[i]}"

    # Check if the sequence and reference files exist
    if [[ -f "$sequences" && -f "$reference" && -f "$genemap" ]]; then
        # Run the Nextstrain command within Docker
        nextalign run \
          --reference="$reference" \
          --genemap="$genemap" \
          --output-all="$output_dir" \
          "$sequences"

        echo "Alignment for DENV serotype ${serotypes[i]} completed."
    else
        echo "Error: Missing data files for DENV serotype ${serotypes[i]}."
    fi
done

set +x  # Turn off debugging information
