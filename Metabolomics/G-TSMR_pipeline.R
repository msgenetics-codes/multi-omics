### GSMR + TSMR pipeline
### DATE: 2021-OCTOBER-04
### AUTHOR: XIN LIN

##
rm(list=ls())


## get the args
args = commandArgs(trailingOnly=TRUE)


### BEGIN ANALYSIS PIPELINE #########################

## create analysis and result paths for GSMR
analysis_path  <- paste(args[1], "/GSMR/analysis/", args[3], "/", args[2], "_", args[3], sep = "") ##directory for analysis
result_path  <- paste(args[1], "/GSMR/results/", args[3], sep = "") ##directory for storing results

## create result path for TSMR
result_path2  <- paste(args[1], "/TSMR/results/", args[3], "/", args[2], "_", args[3], sep = "") ##directory for storing results using clumped SNPs
result_path3  <- paste(args[1], "/TSMR_gsmrSNPs/results/", args[3], "/", args[2], "_", args[3], sep = "") ##directory for storing results using GSMR SNPs

mkdir1 <- paste("mkdir -p ", analysis_path, sep="")
mkdir2 <- paste("mkdir -p ", result_path, sep="")
mkdir3 <- paste("mkdir -p ", result_path2, sep="")
mkdir4 <- paste("mkdir -p ", result_path3, sep="")

system(mkdir1)
system(mkdir2)
system(mkdir3)
system(mkdir4)

## define reference dataset
ref_path  <- c(args[7])  ##reference datafile path include filename
n_ref_path <- c(args[8]) ## Sample size for the reference dataset
n_ref_path <- as.numeric(n_ref_path)

### load library
## install gsmr
##install.packages(c('survey'));
##install.packages("http://cnsgenomics.com/software/gsmr/static/gsmr_1.0.9.tar.gz",repos=NULL,type="source")
library("gsmr")

## see here how to install Two-sample MR: https://mrcieu.github.io/TwoSampleMR/
##remotes::install_github("MRCIEU/TwoSampleMR@0.4.26")
myPaths <- c("/u/xlin8/R/x86_64-pc-linux-gnu-library/3.6")
.libPaths(myPaths)
library(TwoSampleMR)

library(data.table)
library(dplyr)
library(ggplot2)

## define working directory
setwd(args[1])

#####################################################
########## Step1: prepare the data ################
#####################################################

### Required data columns
##exp: Chr, BP, SNP, A1, A2, a1_freq, BETA, SE, P, N
##out: Chr, BP, SNP, A1, A2, a1_freq, BETA, SE, P, N

## load data
exp <- fread(args[4],header=T)
out <- fread(args[5],header=T)

### ensure processed exposure instruments with p<5e-08 & MAF>0.05
exp <- exp[exp$P < 5e-08 & (exp$a1_freq > 0.05 & exp$a1_freq < 0.95), ]
d1 <- exp[nchar(exp$A1)==1 & nchar(exp$A2)==1,]
d1$alleles <- paste(d1$A1, d1$A2, sep="")
exp_use <- d1[d1$alleles!="AT" & d1$alleles!="TA" & d1$alleles!="GC" & d1$alleles!="CG", ]

### merge for GSMR
exp_out_merged <- merge(exp_use,out,by.x=c("Chr", "BP", "SNP"),by.y=c("Chr", "BP", "SNP"))

### flip the beta of outcome (BETA.y) based on the A1&A2 of exposure (denoted here as "A3" and "A4")
## note: AF2 is the frequecy for A3; "B1_correct" is the beta of outcome after flip
source(args[6])

f1 <- alleles_freq("exp_out_merged",A1="A1.y",A2="A2.y",B1="BETA.y",A3="A1.x",A4="A2.x",AF2="a1_freq.x")

dt <- f1[,c("Chr","BP","SNP","A1.x","A2.x","a1_freq.x","BETA.x","SE.x","P.x","N.x","B1_correct","SE.y","P.y","N.y")]
names(dt) <- c("Chr","BP","SNP","A1","A2","a1_freq","bzx","bzx_se","P","bzx_n","bzy","bzy_se","bzy_pval","bzy_n")

## save datafile for plink and gsmr
dtname <- paste(analysis_path, "/" ,args[2], "_", args[3], ".txt", sep = "")
write.table(dt, dtname, col.names=T, row.names=F, quote=F)

#####################################################
################# Step2: Clump #####################
#####################################################
## define working directory
setwd(analysis_path)

## define file paths
snpclump_path <- paste(analysis_path, "/" ,args[2], "_", args[3], "_clump", sep="")
snpclumped_path <- paste(analysis_path, "/" ,args[2], "_", args[3], "_clumped", sep="")

## generate commands for plink clumping
if (nrow(dt) > 0) {
    ## run the command1 and 2 in R
    ## cmd1: plink to clump
    cmd1 <- paste("/u/xlin8/software/plink_1.90 --bfile ", ref_path, " --clump ", dtname, " --clump-field P --clump-kb 1000 --clump-p1 1 --clump-p2 1 --clump-r2 0.05 --out ", snpclump_path, sep="")

    ## cmd2: get clumped SNPs
    cmd2 <- paste("awk '{print $3}' ", snpclump_path, ".clumped > ", snpclumped_path, ".SNPs", sep="")
    ## run the command
    system(cmd1)
    system(cmd2)
    
    ## then check how many SNPs after pruning
    ## note: "awk" printed empty spaces from the clumped file, so reading the clumped file instead to solve this issue
    clumpedfile <- paste(snpclump_path, ".clumped", sep="")
    dt_snps <- fread(clumpedfile)
    snpN <- nrow(dt_snps)
    snpN
    ##
    if (snpN >= 10) {
        ## if there are 10 independent SNPs, then caculate r2
        ## cmd3: caculate r2
        cmd3 <- paste("/u/xlin8/software/plink_1.90 --bfile ", ref_path, " --extract ", snpclumped_path, ".SNPs --r2 square --write-snplist --out ", snpclumped_path, "_SNPs_ld",sep="")
        system(cmd3)
    } else {
        ## print error information
        cmd4 <- paste("echo -n Default independent SNPs is 10. There are only ", snpN," independent SNPs!! > ", result_path, "/", args[2], "_", args[3], "_lessSNPs.file",sep="")
        system(cmd4)
        ## quit the job
        quit()
    }
} else {
    ## print error information
    cmd5 <- paste("echo -n There are no GWAS significant SNPs!! > ", result_path, "/", args[2], "_", args[3], "_noSNPs.file",sep="")
    system(cmd5)
    ## quit the job
    quit()
}

#####################################################
################### Step3: GSMR #####################
#####################################################

## START R script for GSMR ##
ldrho_path <- paste(snpclumped_path, "_SNPs_ld.ld", sep="")
temp_path <- paste(snpclumped_path, "_SNPs_ld.snplist", sep="")

ldrho <- fread(ldrho_path, header=F)
temp <- fread(temp_path, header=F)

snp_coeff_id <- temp$V1

colnames(ldrho) = rownames(ldrho) = snp_coeff_id
dim(ldrho)
ldrho[1:5,1:5]

dt_gsmr = dt[match(snp_coeff_id, dt$SNP),]


### GSMR Analysis
bzx = dt_gsmr$bzx    # SNP effects on the risk factor
bzx_se = dt_gsmr$bzx_se    # standard errors of bzx
bzx_pval = dt_gsmr$P   # p-values for bzx
bzy = dt_gsmr$bzy    # SNP effects on the disease
bzy_se = dt_gsmr$bzy_se    # standard errors of bzy
bzy_pval = dt_gsmr$bzy_pval    # p-values for bzy
n_ref = n_ref_path    # Sample size of the reference sample
gwas_thresh = 5e-8    # GWAS threshold to select SNPs as the instruments for the GSMR analysis
single_snp_heidi_thresh = 0.01    # p-value threshold for single-SNP-based HEIDI-outlier analysis
multi_snp_heidi_thresh = 0.01    # p-value threshold for multi-SNP-based HEIDI-outlier analysis
nsnps_thresh = 10   # the minimum number of instruments required for the GSMR analysis; default is 10
heidi_outlier_flag = T    # flag for HEIDI-outlier analysis
ld_r2_thresh = 0.05    # LD r2 threshold to remove SNPs in high LD
ld_fdr_thresh = 0.05   # FDR threshold to remove the chance correlations between the SNP instruments
gsmr2_beta = 0     # 0 - the original HEIDI-outlier method; 1 - the new HEIDI-outlier method that is currently under development
gsmr_results = gsmr(bzx, bzx_se, bzx_pval, bzy, bzy_se, bzy_pval, ldrho, snp_coeff_id, n_ref, heidi_outlier_flag, gwas_thresh, single_snp_heidi_thresh, multi_snp_heidi_thresh, nsnps_thresh, ld_r2_thresh, ld_fdr_thresh, gsmr2_beta)    # GSMR analysis
filtered_index=gsmr_results$used_index

gsmr_results

## bzx and bzy results for each SNP included in final gsmr analysis
SNP_results <- dt_gsmr[filtered_index]

snp_res_path <- paste(result_path, "/", args[2], "_", args[3], "_SNPs.txt", sep = "")
write.table(SNP_results, file=snp_res_path, col.names=T, row.names=F, quote=F)

## bxy results
Beta_results <- data.frame(gsmr_results$bxy, gsmr_results$bxy_se, gsmr_results$bxy_pval)
Beta_results$nsnp_clumped <- nrow(dt_gsmr)
Beta_results$nsnp_GSMR <- nrow(SNP_results)

beta_res_path <- paste(result_path, "/", args[2], "_", args[3], "_betas.txt", sep = "")
write.table(Beta_results, file=beta_res_path, col.names=T, row.names=F, quote=F)

### Visualisation
#Note: 1.96 is the z-value for calculating 95%CI limits
effect_col = colors()[75]
vals = c(bzx[filtered_index]-1.96*bzx_se[filtered_index], bzx[filtered_index]+1.96*bzx_se[filtered_index])
xmin = min(vals); xmax = max(vals)
vals = c(bzy[filtered_index]-1.96*bzy_se[filtered_index], bzy[filtered_index]+1.96*bzy_se[filtered_index])
ymin = min(vals); ymax = max(vals)

plot_res_path <- paste(result_path, "/", args[2], "_", args[3], "_bxyplot.png", sep = "")
png(filename=plot_res_path, type="cairo", width=4, height=3, units='in', pointsize=5, res=300)
par(mar=c(5,5,4,2))

plot(bzx[filtered_index], bzy[filtered_index], pch=20, cex=1.2, bty="n", cex.axis=1.5, cex.lab=1.5,
        col=effect_col, xlim=c(xmin, xmax), ylim=c(ymin, ymax),
        xlab=eval(bquote(expression(.(args[2])~(italic(b[zx]))))),
        ylab=eval(bquote(expression(.(args[3])~(italic(b[zy]))))))
abline(0, gsmr_results$bxy, lwd=1.5, lty=2, col="dim grey")

nsnps = length(bzx[filtered_index])
for( i in 1:nsnps ) {
    # x axis
    xstart = bzx[filtered_index [i]] - bzx_se[filtered_index[i]]; xend = bzx[filtered_index[i]] + bzx_se[filtered_index[i]]
    ystart = bzy[filtered_index[i]]; yend = bzy[filtered_index[i]]
    segments(xstart, ystart, xend, yend, lwd=1.5, col=effect_col)
    # y axis
    xstart = bzx[filtered_index[i]]; xend = bzx[filtered_index[i]]
    ystart = bzy[filtered_index[i]] - bzy_se[filtered_index[i]]; yend = bzy[filtered_index[i]] + bzy_se[filtered_index[i]]
    segments(xstart, ystart, xend, yend, lwd=1.5, col=effect_col)
}
dev.off()


#####################################################
################### Step4: TSMR #####################
#####################################################

## define exposure and outcome parameters
exposuretraitid = args[2]
exposuretrait = paste(args[2], " || ", args[2], sep="")
outcometraitid = args[3]
outcometrait = paste(args[3], " || ", args[3], sep="")

## define methods to be used
tsmr_methods <- c("mr_egger_regression", "mr_ivw", "mr_weighted_median", "mr_weighted_mode") ##add/remove methods to be used in analysis

## define paths for output files using GSMR SNPs
resultfile = paste(result_path3, "/", args[2], "_", args[3], "_res.txt", sep = "")
resultfile2 = paste(result_path3, "/", args[2], "_", args[3], "_SNP.txt", sep = "")
hetresult=paste(result_path3, "/", args[2], "_", args[3], "_het_stat.txt", sep = "")
pleioresult=paste(result_path3, "/", args[2], "_", args[3], "_pleio_test.txt", sep = "")
singleresult=paste(result_path3, "/", args[2], "_", args[3], "_res_single.txt", sep = "")
looresult=paste(result_path3, "/", args[2], "_", args[3], "_res_loo.txt", sep = "")
steigerresult=paste(result_path3, "/", args[2], "_", args[3], "_steiger_test.txt", sep = "")

### load metabolism GWAS data
use <- dt_gsmr[filtered_index]

### prepare data for TSMR
exp_use <- select(use, SNP, effect_allele.exposure = A1, other_allele.exposure = A2, eaf.exposure = a1_freq, beta.exposure = bzx, se.exposure = bzx_se, pval.exposure = P, samplesize.exposure = bzx_n)
exp_use$id.exposure <- exposuretraitid
exp_use$exposure <- exposuretrait

out_use <- select(use, SNP, effect_allele.outcome = A1, other_allele.outcome = A2, eaf.outcome = a1_freq, beta.outcome = bzy, se.outcome = bzy_se, pval.outcome = bzy_pval, samplesize.outcome = bzy_n)
out_use$id.outcome <- outcometraitid
out_use$outcome <- outcometrait

### TSMR analysis
dat = harmonise_data(exp_use, out_use, action = 2)

### Set seed for reproducibility; weighted mode and weighted median methods produced same beta but slightly different se and p-values
set.seed(111)
res = mr(dat, method_list=tsmr_methods)

### TSMR results
res_use <- select(res, id.exposure, id.outcome, method, nsnp, b, se, pval)

## save results in "resultfile" and clumped snp datafile in "resultfile2"
write.table(res_use, resultfile, col.names=T, row.names=F, quote=T)
write.table(use, resultfile2, col.names=T, row.names=F, quote=T)


### Sensitivity analyses

### (1) Obtain heterogeneity statistics
het_stat <- mr_heterogeneity(dat)
# mr_heterogeneity(dat, method_list=c("mr_egger_regression", "mr_ivw"))

## save heterogeneity statistics
write.table(het_stat, hetresult, col.names=T, row.names=F, quote=T)

### (2) Test directional horizontal pleiotropy using intercept term in MR-Egger
pleio_test <- mr_pleiotropy_test(dat)

## save heterogeneity statistics
write.table(pleio_test, pleioresult, col.names=T, row.names=F, quote=T)

### (3) Single SNP analysis (by default, Wald ratio is used for single-SNP MR; IVW and Egger are used for full MR)
res_single <- mr_singlesnp(dat)
# res_single <- mr_singlesnp(dat, single_method="mr_meta_fixed")
# res_single <- mr_singlesnp(dat, all_method="mr_two_sample_ml")

## save single SNP analysis results
write.table(res_single, singleresult, col.names=T, row.names=F, quote=T)

### (4) Leave-one-out analysis to identify if a single SNP is driving the association (default uses IVW)
res_loo <- mr_leaveoneout(dat)

## save leave-one-out analysis results
write.table(res_loo, looresult, col.names=T, row.names=F, quote=T)

### (5) MR Steiger directionality test of assumed causal direction (instrument->exposure; instrument->exposure->outcome)
steiger_test <- directionality_test(dat)

## Note: test results are obtained from the following codes
#mr_steiger(
#    p_exp = dat$pval.exposure,
#    p_out = dat$pval.outcome,
#    n_exp = dat$samplesize.exposure,
#    n_out = dat$samplesize.outcome,
#    r_xxo = 1,
#    r_yyo = 1,
#    r_exp=0
#)

## save MR Steiger directionality test results
write.table(steiger_test, steigerresult, col.names=T, row.names=F, quote=T)


#####################################################
################### Step4: TSMR #####################
#####################################################

## define exposure and outcome parameters
exposuretraitid = args[2]
exposuretrait = paste(args[2], " || ", args[2], sep="")
outcometraitid = args[3]
outcometrait = paste(args[3], " || ", args[3], sep="")

## define methods to be used
tsmr_methods <- c("mr_egger_regression", "mr_ivw", "mr_weighted_median", "mr_weighted_mode") ##add/remove methods to be used in analysis

## define paths for output files using GSMR SNPs
resultfile = paste(result_path2, "/", args[2], "_", args[3], "_res.txt", sep = "")
resultfile2 = paste(result_path2, "/", args[2], "_", args[3], "_SNP.txt", sep = "")
hetresult=paste(result_path2, "/", args[2], "_", args[3], "_het_stat.txt", sep = "")
pleioresult=paste(result_path2, "/", args[2], "_", args[3], "_pleio_test.txt", sep = "")
singleresult=paste(result_path2, "/", args[2], "_", args[3], "_res_single.txt", sep = "")
looresult=paste(result_path2, "/", args[2], "_", args[3], "_res_loo.txt", sep = "")
steigerresult=paste(result_path2, "/", args[2], "_", args[3], "_steiger_test.txt", sep = "")

### load metabolism GWAS data
use <- dt[match(snp_coeff_id, dt$SNP),]


### prepare data for TSMR
exp_use <- select(use, SNP, effect_allele.exposure = A1, other_allele.exposure = A2, eaf.exposure = a1_freq, beta.exposure = bzx, se.exposure = bzx_se, pval.exposure = P, samplesize.exposure = bzx_n)
exp_use$id.exposure <- exposuretraitid
exp_use$exposure <- exposuretrait

out_use <- select(use, SNP, effect_allele.outcome = A1, other_allele.outcome = A2, eaf.outcome = a1_freq, beta.outcome = bzy, se.outcome = bzy_se, pval.outcome = bzy_pval, samplesize.outcome = bzy_n)
out_use$id.outcome <- outcometraitid
out_use$outcome <- outcometrait

### TSMR analysis
dat = harmonise_data(exp_use, out_use, action = 2)

### Set seed for reproducibility; weighted mode and weighted median methods produced same beta but slightly different se and p-values
set.seed(111)
res = mr(dat, method_list=tsmr_methods)

### TSMR results
res_use <- select(res, id.exposure, id.outcome, method, nsnp, b, se, pval)

## save results in "resultfile" and clumped snp datafile in "resultfile2"
write.table(res_use, resultfile, col.names=T, row.names=F, quote=T)
write.table(use, resultfile2, col.names=T, row.names=F, quote=T)


### Sensitivity analyses

### (1) Obtain heterogeneity statistics
het_stat <- mr_heterogeneity(dat)
# mr_heterogeneity(dat, method_list=c("mr_egger_regression", "mr_ivw"))

## save heterogeneity statistics
write.table(het_stat, hetresult, col.names=T, row.names=F, quote=T)

### (2) Test directional horizontal pleiotropy using intercept term in MR-Egger
pleio_test <- mr_pleiotropy_test(dat)

## save heterogeneity statistics
write.table(pleio_test, pleioresult, col.names=T, row.names=F, quote=T)

### (3) Single SNP analysis (by default, Wald ratio is used for single-SNP MR; IVW and Egger are used for full MR)
res_single <- mr_singlesnp(dat)
# res_single <- mr_singlesnp(dat, single_method="mr_meta_fixed")
# res_single <- mr_singlesnp(dat, all_method="mr_two_sample_ml")

## save single SNP analysis results
write.table(res_single, singleresult, col.names=T, row.names=F, quote=T)

### (4) Leave-one-out analysis to identify if a single SNP is driving the association (default uses IVW)
res_loo <- mr_leaveoneout(dat)

## save leave-one-out analysis results
write.table(res_loo, looresult, col.names=T, row.names=F, quote=T)

### (5) MR Steiger directionality test of assumed causal direction (instrument->exposure; instrument->exposure->outcome)
steiger_test <- directionality_test(dat)

## Note: test results are obtained from the following codes
#mr_steiger(
#    p_exp = dat$pval.exposure,
#    p_out = dat$pval.outcome,
#    n_exp = dat$samplesize.exposure,
#    n_out = dat$samplesize.outcome,
#    r_xxo = 1,
#    r_yyo = 1,
#    r_exp=0
#)

## save MR Steiger directionality test results
write.table(steiger_test, steigerresult, col.names=T, row.names=F, quote=T)

