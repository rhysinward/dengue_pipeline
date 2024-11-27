#!/usr/bin/env bash

if [ ! -d "results" ]; then
    mkdir results
fi
if [ ! -d "results/in" ]; then
    mkdir results/in
    touch "results/in/pre_A.fasta"
    touch "results/in/pre_B.fasta"
    touch "results/in/pre_C.fasta"
    touch "results/in/mark"
fi
if [ ! -d "results/in" ]; then
    echo "Error generating input data."
    exit 1
fi
rm -rf results/out
snakemake --use-conda --cores 1 --configfile=config/.test.yaml _test
