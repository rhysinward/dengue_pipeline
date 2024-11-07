#!/usr/bin/env bash

if [ ! -d "results" ]; then
    mkdir results
fi
if [ ! -d "results/processdenguedata" ]; then
    ln -s ../../ProcessDengueData/results/out results/processdenguedata
    if [ ! -d "results/processdenguedata" ]; then
        echo "No input data found. Please run the ProcessDengueData pipeline first."
        exit 1
    fi
fi
if [ ! -d "results/splitgenomeandqc" ]; then
    ln -s ../../SplitGenomeAndQC/results/out results/splitgenomeandqc
    if [ ! -d "results/splitgenomeandqc" ]; then
        echo "No input data found. Please run the SplitGenomeAndQC pipeline first."
        exit 1
    fi
fi
rm -rf results/out
snakemake --use-conda --cores 1 --configfile=config/.test.yaml _test
