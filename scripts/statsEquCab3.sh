#!/bin/bash

echo "loading mamba"
module load mamba
echo "mamba loaded"

echo "loading environment"
mamba activate /home/tmathieu/.conda/envs/enad5
echo "environment loaded"

dir="EquCab3"
ref="/home/tmathieu/eNAD/resources/GCF_002863925.1_EquCab3.0_genomic.fna"
samples=(horse_01_200 horse_02_205 horse_03_238 horse_04_56 horse_05_608 horse_06_1212)
stats="EquCab3stats.txt"
mkdir $dir

for sample in ${samples[@]}
do
	minimap2 \
      -a \
      -x map-hifi \
      --eqx \
      -t 16 \
      $ref \
      $sample/${sample}_scaffold_chr_only.fa \
      > $dir/${sample}_compref.bam

	echo $sample >> $stats
	bamtools stats -in $dir/${sample}_compref.bam >> $stats
	echo " " >> $stats
done
