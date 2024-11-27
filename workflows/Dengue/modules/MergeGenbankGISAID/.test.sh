#!/usr/bin/env bash

if [ ! -d "results" ]; then
    mkdir results
fi
if [ ! -d "results/genbank" ]; then
    ln -s ../../ProcessGenbankData/results/out results/genbank
fi
if [ ! -d "results/genbank" ]; then
    echo "No input data found. Please run the ProcessGenbankData pipeline first."
    exit 1
fi
if [ ! -d "results/gisaid" ]; then
    ln -s ../../ProcessGisaidData/results/out results/gisaid
fi
if [ ! -d "results/gisaid" ]; then
    echo "No input data found. Please run the ProcessGisaidData pipeline first."
    exit 1
fi
rm -rf results/out
snakemake --use-conda --cores 1 --configfile=config/.test.yaml _test
