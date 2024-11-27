#!/usr/bin/env bash

if [ ! -d "results" ]; then
    mkdir results
fi
if [ ! -d "results/metadata" ]; then
    ln -s ../../FilterGenotypesData/results/out results/metadata
fi
if [ ! -d "results/metadata" ]; then
    echo "No input data found. Please run the FilterGenotypesData pipeline first."
    exit 1
fi
if [ ! -d "results/fasta" ]; then
    ln -s ../../FilterGenotypesData/results/out results/fasta
fi
if [ ! -d "results/fasta" ]; then
    echo "No input data found. Please run the FilterGenotypesData pipeline first."
    exit 1
fi
rm -rf results/out
rm -rf workspace
snakemake --use-conda --cores 1 --configfile=config/dengue_target.yaml _test
