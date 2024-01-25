#!/bin/bash
#
# This script aligns the sequences for each serotype using Nextalign.

usage() {
    echo "Usage: $0 [sequence_dir] [output_dir] [resources_dir]"
    echo "  sequence_dir: Directory containing the input sequences."
    echo "  output_dir: Directory to store the output files."
    echo "  resources_dir: Directory containing the reference genomes and genemap files."
    echo "  If no arguments are provided, the default values will be used."
    echo "  Default sequence_dir: sequences"
    echo "  Default output_dir: results"
    echo "  Default resources_dir: resources"
    exit 0
}

# Check if the help option is provided
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    usage
fi

set -x  # This will print each command and its arguments to the terminal

serotypes=("Dengue_1" "Dengue_2" "Dengue_3" "Dengue_4")
reference_serotype=("1" "2" "3" "4")
sequence_dir=$1
output_dir=$2
resources_dir=$3

# Default sequence folder
if [[ -z "$sequence_dir" ]]; then
    sequence_dir="sequences"
fi

# Default output folder
if [[ -z "$output_dir" ]]; then
    output_dir="results"
fi

# Default resources folder
if [[ -z "$resources_dir" ]]; then
    resources_dir="resources"
fi

# Loop over the serotypes
for i in "${!serotypes[@]}"; do
    # Define the input sequence file and reference sequence file based on the serotype
    sequences="${sequence_dir}/Unaligned_${serotypes[i]}.fasta"
    reference="${resources_dir}/reference_genomes/reference_denv${reference_serotype[i]}.fasta"
    genemap="${resources_dir}/genemap/genemap_denv${reference_serotype[i]}.gff"
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
