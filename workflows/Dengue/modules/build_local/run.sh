#!/usr/bin/env bash

snakemake --snakefile workflow/Snakefile --cores 4 --use-conda dengue_export_all
