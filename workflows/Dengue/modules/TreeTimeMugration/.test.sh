#!/usr/bin/env bash

if [ ! -d "results" ]; then
    mkdir results
fi
if [ ! -d "results/nexus" ]; then
    ln -s ../../PruneTree/results/out results/nexus
    if [ ! -d "results/nexus" ]; then
        echo "No input data found. Please run the PruneTree pipeline first."
        exit 1
    fi
fi
if [ ! -d "results/metadata" ]; then
    ln -s ../../Reformatting/results/out results/metadata
    if [ ! -d "results/metadata" ]; then
        echo "No input data found. Please run the Reformatting pipeline first."
        exit 1
    fi
fi
rm -rf results/out
snakemake --use-conda --cores 1 --configfile=config/.test.yaml _test
