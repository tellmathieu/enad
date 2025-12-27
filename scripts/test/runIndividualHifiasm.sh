#!/bin/bash

module load mamba

mamba activate enad

mkdir outputGenomes/horse_02_205

hifiasm -o outputGenomes/horse_02_205/horse_02_205.asm -t 10 /home/tmathieu/eNAD/fastq/horse_02_205.fastq 
