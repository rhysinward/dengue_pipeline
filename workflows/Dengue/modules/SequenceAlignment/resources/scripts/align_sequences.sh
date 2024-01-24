#!/bin/bash
set -x  # This will print each command and its arguments to the terminal

serotypes=("Dengue_1" "Dengue_2" "Dengue_3" "Dengue_4")
reference_serotype=("1" "2" "3" "4")
sequence_dir=$1
output_dir=$2

# Default sequence folder
if [[ -z "$sequence_dir" ]]; then
    sequence_dir="sequences"
fi

# Default output folder
if [[ -z "$output_dir" ]]; then
    output_dir="results"
fi

# Loop over the serotypes
for i in "${!serotypes[@]}"; do
    # Define the input sequence file and reference sequence file based on the serotype
    sequences="${sequence_dir}/Unaligned_${serotypes[i]}.fasta"
    reference="resources/reference_genomes/reference_denv${reference_serotype[i]}.fasta"
    genemap="resources/genemap/genemap_denv${reference_serotype[i]}.gff"
    output_serotype="${output_dir}/Aligned_${serotypes[i]}"

    # Check if the sequence and reference files exist
    if [[ -f "$sequences" && -f "$reference" && -f "$genemap" ]]; then
        # Run the Nextstrain command within Docker
        nextalign run \
          --reference="$reference" \
          --genemap="$genemap" \
          --output-all="$output_serotype" \
          "$sequences"

        echo "Alignment for DENV serotype ${serotypes[i]} completed."
    else
        echo "Error: Missing data files for DENV serotype ${serotypes[i]}."
    fi
done

set +x  # Turn off debugging information
