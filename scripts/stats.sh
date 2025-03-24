#!/bin/bash

echo "loading mamba"
module load mamba
echo "mamba loaded"

echo "loading environment"
mamba activate /home/tmathieu/.conda/envs/enad5
echo "environment loaded"

echo horse_01_200 > stats.txt
bamtools stats -in snifflesPrimary/horse_01_200.primary_compto5_sorted.bam >> stats.txt
echo " " >> stats.txt

echo horse_02_205 >> stats.txt
bamtools stats -in snifflesPrimary/horse_02_205.primary_compto5_sorted.bam >> stats.txt
echo " " >> stats.txt

echo horse_03_238 >> stats.txt
bamtools stats -in snifflesPrimary/horse_03_238.primary_compto5_sorted.bam >> stats.txt
echo " " >> stats.txt

echo horse_04_56 >> stats.txt
bamtools stats -in snifflesPrimary/horse_04_56.primary_compto5_sorted.bam >> stats.txt
echo " " >> stats.txt

echo horse_05_608 >> stats.txt
bamtools stats -in snifflesPrimary/horse_05_608.primary_compto5_sorted.bam >> stats.txt
echo " " >> stats.txt

echo horse_06_1212 >> stats.txt
bamtools stats -in snifflesPrimary/horse_06_1212.primary_compto5_sorted.bam >> stats.txt
echo " " >> stats.txt