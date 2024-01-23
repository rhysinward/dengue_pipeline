#!/bin/bash

# Define the base directory for the Dengue analysis
BASE_DIR="/Users/rhysinward/Documents/Dengue_anaysis"

# Loop through each Dengue dataset
for i in {1..4}
do
    # Set variables for each dataset
    TREE_FILE="${BASE_DIR}/Dengue_${i}_cleaned.fasta.treefile"
    ALIGNMENT="${BASE_DIR}/Dengue_${i}_cleaned.fasta"
    METADATA="${BASE_DIR}/Dengue_${i}_cleaned_metadata.csv"
    REFERENCE_SEQUENCE="${BASE_DIR}/reference_genomes/reference_dengue_denv${i}.gb"
    OUTPUT_DIR="${BASE_DIR}/dengue_${i}"

    # Create output directory if it doesn't exist
    mkdir -p "${OUTPUT_DIR}"

    # Refine
    augur refine --tree "$TREE_FILE" --alignment "$ALIGNMENT" --metadata "$METADATA" --output-tree "${OUTPUT_DIR}/timetree_denv${i}.nwk" --output-node-data "${OUTPUT_DIR}/branch_lengths.json" --timetree --coalescent opt --date-confidence --date-inference marginal --clock-filter-iqd 3 --root best

    # Ancestral
    augur ancestral --tree "${OUTPUT_DIR}/timetree_denv${i}.nwk" --alignment "$ALIGNMENT" --inference joint

    # Translate
    augur translate --tree "${OUTPUT_DIR}/timetree_denv${i}.nwk" --ancestral-sequences "${OUTPUT_DIR}/Dengue_${i}_cleaned_mutations.json" --reference-sequence "$REFERENCE_SEQUENCE" --output "${OUTPUT_DIR}/aa_muts.json"

    # Traits
    augur traits --tree "${OUTPUT_DIR}/timetree_denv${i}.nwk" --metadata "$METADATA" --output "${OUTPUT_DIR}/traits.json" --columns Country State --confidence

    # Export for visualisation in Auspice
    augur export v2 --tree "${OUTPUT_DIR}/timetree_denv${i}.nwk" --metadata "$METADATA" --node-data "${OUTPUT_DIR}/branch_lengths.json" "${OUTPUT_DIR}/traits.json" "${OUTPUT_DIR}/Dengue_${i}_cleaned_mutations.json" "${OUTPUT_DIR}/aa_muts.json" --auspice-config "${BASE_DIR}/auspice_config.json" --output "${BASE_DIR}/auspice/dengue_${i}.json"
done
