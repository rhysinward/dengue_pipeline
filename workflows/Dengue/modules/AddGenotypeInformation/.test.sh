#!/usr/bin/env bash

if [ ! -d "results" ]; then
    mkdir results
fi
if [ ! -d "results/genotype" ]; then
    ln -s ../../AssignSerotypeAndGenotypeNextclade/results/out results/genotype
fi
if [ ! -d "results/genotype" ]; then
    echo "No input data found. Please run the AssignSerotypeAndGenotypeNextclade pipeline first."
    exit 1
fi
if [ ! -d "results/metadata" ]; then
    ln -s ../../ProcessDengueData/results/out results/metadata
fi
if [ ! -d "results/metadata" ]; then
    echo "No input data found. Please run the ProcessDengueData pipeline first."
    exit 1
fi
rm -rf results/out
snakemake --use-conda --cores 1 --configfile=config/.test.yaml _test
