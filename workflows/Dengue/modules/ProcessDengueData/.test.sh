#!/usr/bin/env bash

if [ ! -d "results" ]; then
    mkdir results
fi
unlink results/in || true
if [ ! -d "results/in" ]; then
    ln -s ../../ProcessGenbankData/results/out results/in
fi
if [ ! -d "results/in" ]; then
    echo "No input data found. Please run the ProcessGenbankData pipeline first."
    exit 1
fi
rm -rf results/out
snakemake --use-conda --cores 1 --configfile=config/.test.yaml _test
