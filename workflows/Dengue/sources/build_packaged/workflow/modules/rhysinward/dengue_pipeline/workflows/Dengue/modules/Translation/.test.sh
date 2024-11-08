#!/usr/bin/env bash

if [ ! -d "results" ]; then
    mkdir results
fi
if [ ! -d "results/treetime" ]; then
    ln -s ../../TreeTime/results/out results/treetime
    if [ ! -d "results/treetime" ]; then
        echo "No input data found. Please run the TreeTime pipeline first."
        exit 1
    fi
fi
if [ ! -d "results/mutations" ]; then
    ln -s ../../Mutations/results/out results/mutations
    if [ ! -d "results/mutations" ]; then
        echo "No input data found. Please run the Mutations pipeline first."
        exit 1
    fi
fi
rm -rf results/out
snakemake --use-conda --cores 1 --configfile=config/.test.yaml _test
