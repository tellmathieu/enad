#!/bin/bash
echo "loading mamba"
module load mamba
echo "mamba loaded"

echo "loading environment"
mamba activate /home/tmathieu/.conda/envs/enad5
echo "environment loaded"

fastqc -o qc -f fastq fastq/horse_01_200.fastq
