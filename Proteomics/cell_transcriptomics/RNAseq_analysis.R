####  single cell RNA QCed UMI matrix: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118257/
### GSE118257_MSCtr_snRNA_ExpressionMatrix_R.txt.gz
### GSE118257_MSCtr_snRNA_FinalAnnotationTable.txt.gz

#### analysis
### 1 ) as the data was QCed, no further QC are needed and directy use the matrix for downstream analysis. We do DEX using [Seurat 3.0] (https://satijalab.org/seurat/)
### 2) as we have the UMI matrix, I run DEX using the nagative binomial regression, adjusted by several covariates indcluding sample ID, number of RNA, and cell cycle score.

#### Read the inputs ####
rm(list=ls())
setwd("~/Desktop/PhD/pQTL_MS/final_results/brain_RNAseq_pQTL/analysis/negb")

##load library
library(Seurat)

## read data
mat = read.table("~/Desktop/PhD/pQTL_MS/final_results/brain_RNAseq_pQTL/data/MSCtr_snRNA_ExpressionMatrix_R.txt",header=T)    # nrow=21581 genes; ncol= 17799 cells #
metadata = read.table("~/Desktop/PhD/pQTL_MS/final_results/brain_RNAseq_pQTL/data/MSCtr_snRNA_FinalAnnotationTable.txt", header=T)

## look at cell number by case status
Freqs <- table(metadata$Celltypes, metadata$Condition) 
Freqs  

#### We focus on cell types with case cell number >5% of total sample and control cell number >5% of total sample; 20 cell types left excluding Immune cells (Ctrl:0; MS:423), Macrophages (Ctrl:1; MS:367) and Vasc_smooth_muscle (Ctrl:3; MS:109).

#### using both case&control cells ####
metadata$Detected = as.character(metadata$Detected)
metadata$Detected = sub(":",".",metadata$Detected)
metadata$Detected = sub("/",".",metadata$Detected)
mat2 = mat[,colnames(mat) %in% metadata$Detected]                            # nrow=21581 genes; ncol= 17799 cells in total#

#### make a seurat object ####
obj = CreateSeuratObject(counts = mat2, project = "obj", min.cells = 0, min.features = 0)

#### calculate cell cycle score ####
# cell cycle score will be used as covariates in DEX analysis
s.genes = cc.genes$s.genes
g2m.genes = cc.genes$g2m.genes
obj = CellCycleScoring(object = obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
obj@meta.data$Detected = row.names(obj@meta.data)

metadata_allsample = merge(obj@meta.data, metadata, by='Detected')
dim(metadata_allsample)       ## 17799

#### for-loop to loop though each cell type and compare expressions in case&control ####

celllist = unique(metadata$Celltypes)
numlist <- c(1:18, 21, 23)

for (i in numlist){
    cellname <- celllist[i]
    metadata_temp <- metadata_allsample[metadata_allsample$Celltypes==cellname,]
    mat_temp = mat[,colnames(mat) %in% metadata_temp$Detected]  # testing numlist==6: nrow=21581 genes; ncol= 1046 cells for Astrocytes#
    
    #### make object ####
    obj = CreateSeuratObject(counts = mat_temp, project = "obj", min.cells = 0, min.features = 0)
    obj@meta.data$Detected = row.names(obj@meta.data)
    obj_meta_upd = merge(obj@meta.data, metadata_temp, by='Detected')

    #### update the meta_data and set case status as indents ####
    obj@meta.data = obj_meta_upd
    row.names(obj@meta.data) = obj@meta.data$Detected
    uu=as.factor(obj@meta.data$Condition)    ## setting case status as indents
    names(uu)=obj@meta.data$Condition
    Idents(obj)=uu

    #### compile list of names for the 21581 genes
    gene_list <- unlist(obj@assays$RNA@counts@Dimnames[1])
    

    #### DEX analysis, for loop to do the analysis comparing expressions in cases and controls for each cell type
    #### the analysis for each cell type will take ~ 2-7 mins
    
    ## If choosing "MS" as ident.1; A positive logFC means that there is overexpression of the gene in MS cells compared to control, and vice versa
    ## Different models such as 'MAST' and 'LR' can be used instead simply substituting 'negbinom' in test.use
    ## pct.1: The percentage of cells where the gene is detected in the cases
    ## pct.2: The percentage of cells where the gene is detected in the control
    ## 

    cluster1.markers = FindMarkers(obj, ident.1 = "MS", test.use = 'negbinom', min.pct = 0.05, logfc.threshold = 0.1, latent.vars = c("Sample","nCount_RNA","S.Score","G2M.Score"))
    cluster1.markers$gene=row.names(cluster1.markers)
    markers_res = cluster1.markers[(cluster1.markers$gene %in% gene_list),]
    print(markers_res)
    filename <- paste("RNAseq_markers_res_",i,"_updated.RData", sep="")
    save(markers_res, file = filename)
}

