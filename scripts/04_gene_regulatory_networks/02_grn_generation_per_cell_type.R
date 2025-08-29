#************************************************#
# GRN generation                                 #
# script: 02_grn_generation_per_cell_ttype.R     #
#************************************************#

# GRN generation

# Load required libraries

# Define the directory were input data from OSF has been downloaded
data_wd <- setwd("../GitHub_code/")

##>>>>>>>>>>>>>>>>>>>>> Inputs required

# 1- Normalized RNA data
# Gene Expression data
rna <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/rna_norm_data.rds"))
rna_metadata <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/rna_metadata.rds"))
rna_anno <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/rna_gene_annotation.rds"))

# Chromatin accessibility data
atac <- readRDS(paste0(data_wd,"/osfstorage-archive/01_atac_seq_data/atac_norm_data.rds"))
atac_metadata <- readRDS(paste0(data_wd,"/osfstorage-archive/01_atac_seq_data/atac_metadata.rds"))
atac_anno <- readRDS(paste0(data_wd,"/osfstorage-archive/01_atac_seq_data/atac_OCR_annotation.rds"))

# TF-OCR-gene associations
tf_ocr_gene <- readRDS(paste0(data_wd,"/osfstorage-archive/04_gene_regulatory_networks/tf_ocr_gene_interactions_curated.rds"))
#-------------------------------------


## >> Get the GRN for each individual cell-type ------------------------------------------------------------------##

# 1. Compute the median gene expression and atac accessibility per group

median_atac_per_group <- t(apply(atac,1,FUN=function(x){return(by(x,atac_metadata$CellType, median))}))

median_rna_per_group <- t(apply(rna,1,FUN=function(x){return(by(x,rna_metadata$CellType, median))}))


# 2. Merge the TF-OCR_gene information with the median rna and atac data

tmp1 <- median_atac_per_group[tf_ocr_gene$ocr,]
colnames(tmp1) <- paste0("atac_",colnames(tmp1))

tmp2 <- median_rna_per_group[tf_ocr_gene$tf_ensembl,]
colnames(tmp2) <- paste0("tf_",colnames(tmp2))

tmp3 <- median_rna_per_group[tf_ocr_gene$gene_ensembl,]
colnames(tmp3) <- paste0("gene_",colnames(tmp3))

tf_ocr_gene_extended <- cbind(tf_ocr_gene,tmp1,tmp2,tmp3) 


# 3. Check if the TF-binding is possible in each cell type
tf_ocr_gene_extended$HSC_interaction <- rep("yes",nrow(tf_ocr_gene_extended))
tf_ocr_gene_extended$HSC_interaction[tf_ocr_gene_extended$atac_HSC<1] <- "no"
tf_ocr_gene_extended$CLP_interaction <- rep("yes",nrow(tf_ocr_gene_extended))
tf_ocr_gene_extended$CLP_interaction[tf_ocr_gene_extended$atac_CLP<1] <- "no"
tf_ocr_gene_extended$proB_interaction <- rep("yes",nrow(tf_ocr_gene_extended))
tf_ocr_gene_extended$proB_interaction[tf_ocr_gene_extended$atac_proB<1] <- "no"
tf_ocr_gene_extended$preB_interaction <- rep("yes",nrow(tf_ocr_gene_extended))
tf_ocr_gene_extended$preB_interaction[tf_ocr_gene_extended$atac_preB<1] <- "no"
tf_ocr_gene_extended$ImmatureB_interaction <- rep("yes",nrow(tf_ocr_gene_extended))
tf_ocr_gene_extended$ImmatureB_interaction[tf_ocr_gene_extended$atac_ImmatureB<1] <- "no"
tf_ocr_gene_extended$Transitional_B_interaction <- rep("yes",nrow(tf_ocr_gene_extended))
tf_ocr_gene_extended$Transitional_B_interaction[tf_ocr_gene_extended$atac_Transitional_B<1] <- "no"
tf_ocr_gene_extended$Naive_CD5pos_interaction <- rep("yes",nrow(tf_ocr_gene_extended))
tf_ocr_gene_extended$Naive_CD5pos_interaction[tf_ocr_gene_extended$atac_Naive_CD5pos<1] <- "no"
tf_ocr_gene_extended$Naive_CD5neg_interaction <- rep("yes",nrow(tf_ocr_gene_extended))
tf_ocr_gene_extended$Naive_CD5pos_interaction[tf_ocr_gene_extended$atac_Naive_CD5pos<1] <- "no"

saveRDS(tf_ocr_gene_extended,"tf_ocr_gene_interactions_curated_extended.rds")
#-------------------------------------
