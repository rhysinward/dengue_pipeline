#!/usr/bin/env bash

if [ ! -d "results" ]; then
    mkdir results
fi
if [ ! -d "results/generator" ]; then
    ln -s ../../PartitionIntoGenotypes/results/out results/generator
fi
if [ ! -d "results/generator" ]; then
    echo "No input data found. Please run the PartitionIntoGenotypes pipeline first."
    exit 1
fi
if [ ! -d "results/target" ]; then
    ln -s ../../SubsampleLineages/results/out results/target
fi
if [ ! -d "results/target" ]; then
    echo "No input data found. Please run the SubsampleLineages pipeline first."
    exit 1
fi
rm -rf results/out
snakemake --use-conda --cores 1 --configfile=config/.test.yaml _test
