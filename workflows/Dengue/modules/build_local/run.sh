#!/usr/bin/env bash

snakemake --snakefile workflow/Snakefile --cores 1 --use-conda dengue_export_all
