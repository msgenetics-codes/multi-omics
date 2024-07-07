### G-TSMR_pipeline_caller for multiple exposure & outcome phenotypes
### DATE: 2021-OCT-15
### AUTHOR: XIN LIN


## Set directory
cd /u/xlin8/metabolism/pipeline2


for k in aa ab ac ad ae
do

### prepare shell scripts for jobs
cat << __EOT__ > G-TSMR_MS_run_list${k}.sh


#!/bin/bash

#PBS -l mem=20g
#PBS -l ncpus=1
#PBS -l walltime=48:00:00

cd /u/xlin8/metabolism/pipeline2

## load packages
module load rosalind/1.0
module load gcc-env/6.4.0
module load openmpi-env/3.1.1
module load R/3.6.1

## get file lists (splitting 174 metabolites into five file lists "aa"-"ae")
##awk '{print $3}' /data/menzies_projects/ms-epi/final/metabolism/data/file_downloaded.infor_all_summary > file.list

##split -l 35 file.list file.list

### exposures #####
for j in \`cat file.list${k}\`
do

## usage: Rscript cause_model.R arg1 arg2 arg3 arg4 arg5 arg6 arg7 arg8
## arg1: working directory for pipeline
## arg2: exposure trait name
## arg3: outcome trait name
## arg4: exposure datafile path include filename
## arg5: outcome datafile path include filename
## arg6: file path for R script matching A1&A2 alleles across datasets
## arg7: file path for reference dataset
## arg8: sample size for reference dataset

### Required data columns
##exp: Chr, BP, SNP, A1, A2, a1_freq, BETA, SE, P, N
##out: Chr, BP, SNP, A1, A2, a1_freq, BETA, SE, P, N


Rscript /u/xlin8/metabolism/pipeline2/G-TSMR_pipeline.R \
/u/xlin8/metabolism/pipeline2 \
\${j} \
MS \
/data/menzies_projects/ms-epi/final/metabolism/data/data_processed/\${j}.use.gz \
/u/xlin8/MS_V5matched/N15/MS_1000G_N15_forMR.txt \
/u/xlin8/metabolism/pipeline/alleles_match_noATCG_frequency.R \
/data/menzies_projects/ms-epi/final/1000G_ref/1000G_20130502_v5_full_EUR \
503

done

__EOT__

### Submit the jobs
qsub G-TSMR_MS_run_list${k}.sh

done
