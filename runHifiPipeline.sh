#!/bin/bash
echo "loading mamba"
module load mamba
echo "mamba loaded"

echo "loading environment"
mamba activate /home/tmathieu/.conda/envs/enad5
echo "environment loaded"

echo "running snakemake"
snakemake --snakefile phased_analysis_vcf.smk -j 100 -p --resources mem_mb=100000 --rerun-incomplete