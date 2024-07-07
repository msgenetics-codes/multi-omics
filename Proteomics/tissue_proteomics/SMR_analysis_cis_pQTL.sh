#!/bin/bash

#PBS -l mem=250g
#PBS -l ncpus=1
#PBS -l walltime=224:00:00

### INFORMATION
### DATE: 2020-APRIL-14
### AUTHOR: YUAN ZHOU
### LAST MODIFIED BY XIN LIN: 2021-FEB-23

cd /u/xlin8/pQTL_final/results

####################### S2: this is to use pQTL file

/u/xlin8/software/SMR_V1.03 \
--bfile /data/menzies_projects/ms-epi/final/1000G_ref/1000G_20130502_v5_full_EUR \
--gwas-summary /data/menzies_projects/ms-epi/final/MS_GWAS/MS_GWAS_V3_add1000G_EURfreq_SMR.txt \
--beqtl-summary /data/menzies_projects/ms-epi/final/pQTL/data/pQTL_Nature_2976_SOMAMERID \
--exclude-probe /data/menzies_projects/ms-epi/final/pQTL/data/MHC_pQTL.probe \
--diff-freq 0.4 \
--diff-freq-prop 0.05 \
--out smr_cis_pQTL_all

