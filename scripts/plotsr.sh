#!/bin/bash

module load mamba

mamba activate enad4

#--tracks /home/tmathieu/eNAD/resources/tracksPOI.txt \
mkdir -p /group/ctbrowngrp/finnolab_shared/eNAD_PacBio/plotsrGraphsPOI

plotsr --sr /group/ctbrowngrp/finnolab_shared/eNAD_PacBio/horse_05_608/syri/reference_horse_05_608syri.out \
--sr /group/ctbrowngrp/finnolab_shared/eNAD_PacBio/chainAlign/horse_05_608-horse_04_56syri/horse_05_608-horse_04_56syri.out \
--sr /group/ctbrowngrp/finnolab_shared/eNAD_PacBio/chainAlign/horse_04_56-horse_06_1212syri/horse_04_56-horse_06_1212syri.out \
--sr /group/ctbrowngrp/finnolab_shared/eNAD_PacBio/chainAlign/horse_06_1212-horse_03_238syri/horse_06_1212-horse_03_238syri.out \
--sr /group/ctbrowngrp/finnolab_shared/eNAD_PacBio/chainAlign/horse_03_238-horse_01_200syri/horse_03_238-horse_01_200syri.out \
--sr /group/ctbrowngrp/finnolab_shared/eNAD_PacBio/chainAlign/horse_01_200-horse_02_205syri/horse_01_200-horse_02_205syri.out \
--genomes /group/ctbrowngrp/finnolab_shared/eNAD_PacBio/syri_genomes.txt \
--markers /home/tmathieu/eNAD/resources/points_of_interest_historical.bed \
-H 8 \
-W 12 \
--chr chr7 \
-o /group/ctbrowngrp/finnolab_shared/eNAD_PacBio/plotsrGraphsPOI/syri_chr16_combined_output_plot.png