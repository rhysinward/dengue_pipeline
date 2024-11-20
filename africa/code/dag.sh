#!/bin/bash

# Activate the appropriate Conda environment
conda activate dag  # Replace with your environment name

# Generate the DAG
snakemake --dag | dot -Tpdf > workflow_dag.pdf

echo "DAG generated and saved to workflow_dag.pdf"
