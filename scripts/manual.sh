#!/bin/bash
echo "loading mamba"
module load mamba
echo "mamba loaded"

echo "loading environment"
mamba activate /home/tmathieu/.conda/envs/enad5
echo "environment loaded"

wild='horse_06_1212'

minimap2 \
	-a \
	-x map-hifi \
	--eqx \
	-t 20 \
	resources/reference_genome_chr_only.fa \
	fastq/${wild}.fastq \
	> refAlign/${wild}.sam

samtools view -@ 20 -S -b refAlign/${wild}.sam > refAlign/temp_${wild}.bam

rm refAlign/${wild}.sam

samtools sort -@ 20 -o refAlign/${wild}.sorted.bam refAlign/temp_${wild}.bam
wait
samtools index -@ 20 -b refAlign/${wild}.sorted.bam

#rm refAlign/${wild}.sam
#rm refAlign/temp_${wild}.bam
