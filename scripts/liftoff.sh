#!/bin/bash
echo "loading mamba"
module load mamba
echo "mamba loaded"

echo "loading environment"
mamba activate /home/tmathieu/.conda/envs/enad5
echo "environment loaded"

gff="resources/EquCab3.0_annotation_chr_only.gff"
hap="horse_05_608_hap2.gff"
ref="resources/reference_genome_chr_only.fa"
target="horse_05_608.hap2/horse_05_608.hap2_scaffold_chr_only.fa"

liftoff -g $gff -o $hap -u ${hap}_unmapped.txt $target $ref
