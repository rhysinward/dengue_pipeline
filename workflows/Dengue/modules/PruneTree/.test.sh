#!/usr/bin/env bash

if [ ! -d "results" ]; then
    mkdir results
fi
if [ ! -d "results/in" ]; then
    ln -s ../../TreeTimeBioconda/results/out results/in
    if [ ! -d "results/in" ]; then
        echo "No input data found. Please run the TreeTimeBioconda pipeline first."
        exit 1
    fi
fi
rm -rf results/out
snakemake --use-conda --cores 1 --configfile=config/.test.yaml _test
