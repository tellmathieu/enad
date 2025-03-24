#!/bin/sh

module load mamba

mamba activate enad4

ragtag.py correct /home/tmathieu/eNAD/resources/GCF_002863925.1_EquCab3.0_genomic.fna outputGenomes/horse_01_200/horse_01_200.asm.bp.p_ctg.fa -o horse_01_200 -u
