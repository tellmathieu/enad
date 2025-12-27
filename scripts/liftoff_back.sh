#!/bin/bash
echo "loading mamba"
module load mamba
echo "mamba loaded"

echo "loading environment"
mamba activate /home/tmathieu/.conda/envs/enad5
echo "environment loaded"

liftoff -g resources/horse_05_608_hap1.gff -o backtoEquCab3.0_hap1.gff -u backtoEquCab3.0_hap1_unmapped.txt resources/reference_genome_chr_only.fa horse_05_608.hap2/horse_05_608.hap2_scaffold_chr_only.fa
