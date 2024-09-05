#!/bin/bash

sbatch /project2/gilad/umans/oxygen_eqtl/data/snakemake_cellranger.batch \
"/scratch/midway2/umans/miniconda3/envs/chromium" \
"-s /project2/gilad/umans/oxygen_eqtl/data/Snakefile_cellranger2" \
"--configfile /project2/gilad/umans/oxygen_eqtl/data/config.yaml" \
"--config proj_dir=/project2/gilad/umans/oxygen_eqtl/data/" \
"--cluster-config /project2/gilad/umans/oxygen_eqtl/data/cluster.json"
