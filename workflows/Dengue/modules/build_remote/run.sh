#!/usr/bin/env bash

CONDA_SUBDIR=osx-64 snakemake --snakefile workflow/Snakefile --cores 4 --use-conda dengue_export_all
