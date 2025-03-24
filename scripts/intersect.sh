#!/bin/bash

echo "loading mamba"
module load mamba
echo "mamba loaded"

echo "loading environment"
mamba activate /home/tmathieu/.conda/envs/enad5
echo "environment loaded"

gff="horse_05_608_hap2.gff"
vcf="sniffles/phasedSVs.compto5-2_filtered.vcf"
outgff="horse_05_608_hap2_overlap_vcf.txt"

bedtools intersect -a $gff -b $vcf -wo > $outgff