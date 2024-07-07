#!/bin/bash

#PBS -l mem=40g
#PBS -l ncpus=1
#PBS -l walltime=224:00:00

### INFORMATION
### DATE: 2022-MAR-10
### AUTHOR: XIN LIN

####################### S1: this is to excluede MHC region
##awk '$1==6 && $4>=24000000 && $4<=36000000 {print $2}' /data/menzies_projects/ms-epi/final/eQTL/brain_eQTL/Brain-eMeta/Brain-eMeta.epi > /u/xlin8/metabolism/smr/data/MHC_brain_eQTLGEN.probe


####################### S2: Perform SMR eQTL analysis

cd /u/xlin8/proteomics/results

/u/xlin8/software/SMR_V1.03 \
--bfile /data/menzies_projects/ms-epi/final/1000G_ref/1000G_20130502_v5_full_EUR \
--gwas-summary /data/menzies_projects/ms-epi/final/MS_GWAS/MS_GWAS_V3_add1000G_EURfreq_SMR.txt \
--beqtl-summary /data/menzies_projects/ms-epi/final/eQTL/brain_eQTL/Brain-eMeta/Brain-eMeta \
--exclude-probe /u/xlin8/metabolism/smr/data/MHC_brain_eQTLGEN.probe \
--diff-freq 0.1 \
--diff-freq-prop 0.05 \
--out smr_brain_cis_eQTL_all
