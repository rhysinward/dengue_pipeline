#!/usr/bin/env bash

snakemake --snakefile workflow/Snakefile --cores 1 --use-conda extractphylotree_target plotexportsandimports_target
