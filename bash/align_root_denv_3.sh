#!/bin/bash
set -x  # This will print each command and its arguments to the terminal

# Define the specific file names
sequence="reference_genomes/rooting_genotype_4_dengue3.fasta"
reference="reference_genomes/reference_denv3.fasta"  # Adjust if your reference genome has a different name
genemap="genemap/genemap_denv3.gff"  # Adjust if your genemap file has a different name
output_dir="reference_genomes/aligned_root_dengue3"

# Check if the sequence, reference, and genemap files exist
if [[ -f "$sequence" && -f "$reference" && -f "$genemap" ]]; then
    # Run the Nextstrain command within Docker
    docker run --rm -v "$(pwd)":/data nextstrain/base \
    nextalign run \
      --reference="/data/$reference" \
      --genemap="/data/$genemap" \
      --output-all="/data/$output_dir" \
      "/data/$sequence"

    echo "Alignment for $sequence completed."
else
    echo "Error: Missing one or more data files for $sequence."
fi

set +x  # Turn off debugging information
