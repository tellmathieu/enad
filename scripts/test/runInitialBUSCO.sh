#!/bin/bash

module load mamba

mamba activate busco


busco -c 20 --metaeuk --metaeuk_parameters="--disk-space-limit=1000M,--remove-tmp-files=1" --metaeuk_rerun_parameters="--disk-space-limit=1000M,--remove-tmp-files=1" -i outputGenomes/horse_01_200/horse_01_200.asm.bp.p_ctg.fa -o horse_01_200 -m genome -f -l laurasiatheria_odb10
