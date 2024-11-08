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
if [ ! -d "results/reformatting" ]; then
    ln -s ../../Reformatting/results/out results/reformatting
    if [ ! -d "results/reformatting" ]; then
        echo "No input data found. Please run the Reformatting pipeline first."
        exit 1
    fi
fi
rm -rf results/out
snakemake --use-conda --cores 1 --configfile=config/.test.yaml _test
