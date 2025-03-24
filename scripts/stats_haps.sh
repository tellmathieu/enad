#!/bin/bash
echo "loading mamba"
module load mamba
echo "mamba loaded"

echo "loading environment"
mamba activate /home/tmathieu/.conda/envs/enad5
echo "environment loaded"

list=("horse_01_200.compto5-1" "horse_01_200.compto5-2" "horse_02_205.compto5-1" "horse_02_205.compto5-2" "horse_03_238.compto5-1" "horse_03_238.compto5-2" "horse_04_56.compto5-1" "horse_04_56.compto5-2" "horse_05_608.compto5-1" "horse_05_608.compto5-2" "horse_06_1212.compto5-1" "horse_06_1212.compto5-2")

for wild in ${list[@]};
do
	echo ${wild}
	echo sniffles/${wild}_sorted.bam
	cov=$(samtools depth -a sniffles/${wild}_sorted.bam | awk '{sum+=$3} END {print sum/NR}') 
	echo ${wild}	coverage	${cov} >> sniffles/refStats.txt
	bamtools stats -in sniffles/${wild}_sorted.bam >> sniffles/refStats.txt
done
