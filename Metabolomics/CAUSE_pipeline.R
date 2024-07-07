### INFORMATION
### DATE: 2021-OCTOBER-15
### AUTHOR: XIN LIN

##
rm(list=ls())

## get the args
args = commandArgs(trailingOnly=TRUE)


### BEGIN ANALYSIS PIPELINE #########################

## create analysis and result paths for CAUSE
result_path  <- paste(args[1], "/CAUSE/results/", args[3], sep = "") ##directory for storing results

mkdir1 <- paste("mkdir -p ", result_path, sep="")

system(mkdir1)

## load library
library(cause)
library(data.table)
library(dplyr)

## define working directory
setwd(args[1])

## load data
exp <- fread(args[4],header=T)
out <- fread(args[5],header=T)

bfile_path <- c(args[6])

## merge the two data
exp_out_merged <- gwas_merge(exp, out, snp_name_cols = c("SNP", "SNP"),
beta_hat_cols = c("BETA", "BETA"),
se_cols = c("SE", "SE"),
A1_cols = c("A1", "A1"),
A2_cols = c("A2", "A2"))

### calculate nuisance parameters with a random subset of 1million variants
set.seed(100)
varlist <- with(exp_out_merged, sample(snp, size=1000000, replace=FALSE))
params <- est_cause_params(exp_out_merged, varlist)

variants <- exp_out_merged %>% mutate(pval1 = 2*pnorm(abs(beta_hat_1/seb1), lower.tail=FALSE))

### ld pruning
## devtools::install_github("explodecomputer/genetics.binaRies")

exp_out_clumped <- variants %>% rename(rsid = snp,
                  pval = pval1) %>%
           ieugwasr::ld_clump(dat = .,
                     clump_r2 = 0.05,
                     clump_p = args[7],               ########default clump_p = 1e-03###
                     plink_bin = genetics.binaRies::get_plink_binary(),
                     bfile = bfile_path)
## cause
keep_snps <- exp_out_clumped$rsid
res <- cause(X=exp_out_merged, variants = keep_snps, param_ests = params)

## filename
filename1 <- paste(result_path, "/", args[2], "_", args[3], "_res.RData",sep="")
filename2 <- paste(result_path, "/", args[2], "_", args[3], "_results.txt",sep="")
## save the results
save(res, file = filename1)

## prepare a data file
dt <- summary(res)
df <- data.frame(dt$tab)
use <- data.frame(p_fit=dt$p,gamma_sharing=df$gamma[1],eta_sharing=df$eta[1],q_sharing=df$q[1],gamma_causal=df$gamma[2],eta_causal=df$eta[2],q_causal=df$q[2],nsnp=length(res$data$snp))

## write the output
write.table(use,file=filename2,quote=T,sep=",",col.names=T,row.names=F)


