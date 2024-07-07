### this is an R function to match the allele frequencies of data A using the reference data (e.g.1000G data)

### datause is the merged data of data A and reference data B

### A1 refers to allele in data A
### A2 refers to allele in data A
### B1 refers to the beta in data A

### A3 refers to allele in reference data B (e.g 1000G)
### A4 refers to allele in reference data B (e.g 1000G)
### AF2 refers to allele A3 frequency in reference data B (e.g 1000G)

## exclude the ATCG

alleles_freq <- function(datause,A1,A2,B1,A3,A4,AF2) {
         dt <- get(datause)
         dt$comb <- as.character(paste(dt[[A1]],dt[[A2]],dt[[A3]],dt[[A4]],sep=""))
         dt[[A1]] <- as.character(dt[[A1]])
         dt[[A2]] <- as.character(dt[[A2]])
         dt[[A3]] <- as.character(dt[[A3]])
         dt[[A4]] <- as.character(dt[[A4]])
         
         ## using the references data (e.g. 1000G)
         dt$B1_correct <- NA
   
         #### A) SNPs: subset the data for SNPs
         
         ## dt1: SNPs with the same alleles and direction between data A and reference data
         dt1 <- dt[dt$comb %in% c("AAAA","ACAC","AGAG","TTTT","TCTC","TGTG","CACA","CTCT","CCCC","GAGA","GTGT","GGGG"),]
         dt1$B1_correct <- dt1[[B1]]
         print(paste("SNPs with the same alleles and direction between data A and reference data N = ",dim(dt1)[1]))
        
         ## dt2: SNPs with different alleles but after flip, the alleles and direction are the same
         dt2 <- dt[dt$comb %in% c("AATT","ACTG","AGTC","TTAA","TCAG","TGAC","CAGT","CTGA","CCGG","GACT","GTCA","GGCC"),]
         dt2$B1_correct <- dt2[[B1]]
         print(paste("SNPs with different alleles but after flip, the alleles and direction are the same N = ",dim(dt2)[1]))

         ## dt3: SNPs with the same alleles and but are reversed between data A and reference data
         dt3 <- dt[dt$comb %in% c("ACCA","AGGA","TCCT","TGGT","CAAC","CTTC","GAAG","GTTG"),]
         dt3$B1_correct <- dt3[[B1]] * (-1)
         print(paste("SNPs with the same alleles and but are reversed between data A and reference data N = ",dim(dt3)[1]))

         ## dt4: SNPs with the different alleles and  are reversed between data A and reference data
         dt4 <- dt[dt$comb %in% c("ACGT","AGCT","TCGA","TGCA","CATG","CTAG","GATC","GTAC"),]
         dt4$B1_correct <- dt4[[B1]] * (-1)
         print(paste("SNPs with the different alleles and  are reversed between data A and reference data N = ",dim(dt4)[1]))
         
         ## merge the data
         matched <- rbind(dt1,dt2,dt3,dt4)
         print(paste("SNPs that are keeped after checking, N = ",dim(matched)[1]))
                  
         ## return the data
         return(matched)
}
