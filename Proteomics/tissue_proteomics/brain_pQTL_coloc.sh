#!/bin/bash

#PBS -l mem=10g
#PBS -l ncpus=1
#PBS -l walltime=24:00:00

### INFORMATION
### DATE: 2021-FEBRUARY-13
### AUTHOR: XIN LIN

#### FUSION to do brain pQTL analysis

## Download and unpack the  FUSION software package from github:
##wget https://github.com/gusevlab/fusion_twas/archive/master.zip
##unzip master.zip
##cd fusion_twas-master

## Download and unpack the (1000 Genomes)  LD reference data:
##wget https://data.broadinstitute.org/alkesgroup/FUSION/LDREF.tar.bz2
##tar xjvf LDREF.tar.bz2

## Launch R and install required libraries:
## module load R/3.6.2


###  discovery
for i in {1..22}
do

Rscript /u/xlin8/software/fusion_twas-master/FUSION.assoc_test.R \
--sumstats /u/xlin8/pQTL_final/data/MS_GWAS_FUSION.txt \
--weights /u/xlin8/pQTL_final/data/train_weights.pos_376 \
--weights_dir /data/menzies_projects/ms-epi/final/pQTL/brain_pQTL/ROSMAP.n376.fusion.WEIGHTS \
--ref_ld_chr /u/xlin8/software/fusion_twas-master/LDREF/1000G.EUR. \
--chr $i \
--coloc_P 1e-03 \
--GWASN 41505 \
--out /u/xlin8/pQTL_final/results/brain_pQTL/brain_pQTL.chr$i.dat_coloc

done

###  discovery results
cd /u/xlin8/pQTL_final/results/brain_pQTL
cat *.dat_coloc > /u/xlin8/pQTL_final/results/brain_pQTL/all.discovery_coloc.txt

###  replication
for i in {1..22}
do

Rscript /u/xlin8/software/fusion_twas-master/FUSION.assoc_test.R \
--sumstats /u/xlin8/pQTL_final/data/MS_GWAS_FUSION.txt \
--weights /u/xlin8/pQTL_final/data/train_weights.pos_152 \
--weights_dir /data/menzies_projects/ms-epi/final/pQTL/brain_pQTL/Banner.n152.fusion.WEIGHTS \
--ref_ld_chr /u/xlin8/software/fusion_twas-master/LDREF/1000G.EUR. \
--chr $i \
--coloc_P 1e-03 \     ###added this command to perform COLOC analysis
--GWASN 41505 \
--out /u/xlin8/pQTL_final/results_valid/brain_pQTL.chr$i.dat_coloc

done

###  replication results
cd /u/xlin8/pQTL_final/results_valid
cat *.dat_coloc > /u/xlin8/pQTL_final/results_valid/all.replicate_coloc.txt


