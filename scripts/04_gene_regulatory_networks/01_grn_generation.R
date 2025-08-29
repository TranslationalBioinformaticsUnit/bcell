#*********************************#
# GRN generation                  #
# script: 01_grn_generation.R     #
#*********************************#

# GRN generation

# Load required libraries
library(reshape2)

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

# OCR-gene associations
ocr_gene <- readRDS(paste0(data_wd,"/osfstorage-archive/02_cCREs_OCR_gene_links/candidate_CREs.rds"))

# OCR-TF (motif) binary matrix
ocr_tf_matrix <- read.table(paste0(data_wd,"/osfstorage-archive/03_TFs_OCR_links/ocr_tf_matrix.txt"), sep="\t", header=TRUE, check.names = FALSE)
# Remove those OCRs regions with no TF binding
sel_ocr <- which(rowSums(ocr_tf_matrix)==0)
ocr_tf_matrix <- ocr_tf_matrix[-sel_ocr,]

# 4 - Differential Analysis clustering
cluster_rna <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/differential_expression_analysis/rna_deg_clustering.rds"))

# 5 - Differential accessibility clustering
cluster_atac <- readRDS(paste0(data_wd,"/osfstorage-archive/01_atac_seq_data/differential_binding_analysis/atac_dar_clustering.rds"))
#-------------------------------------


## >> TF-OCR-gene data frame ------------------------------------------------------------------##

# Identify common OCRs between TF-OCR and OCR-gene pairs.
target_OCRs <- intersect(rownames(ocr_tf_matrix),unique(ocr_gene$peak)) #24,244

# Subset ocr-tf matrix and ocr-gene data frame
ocr_tf_matrix2 <- ocr_tf_matrix[target_OCRs,] #26,052
ocr_gene2 <- ocr_gene[ocr_gene$peak%in%target_OCRs,]

# Which TF-OCRs are linked to a gene?
# General:
dim(atac_anno[rownames(atac_anno)%in%rownames(ocr_tf_matrix) & atac_anno$annotation_simplified=="Promoter",])
length(unique(atac_anno$ENSEMBL[rownames(atac_anno)%in%rownames(ocr_tf_matrix) & atac_anno$annotation_simplified=="Promoter"]))

# Within genes exressed in rna:
length(unique(ocr_gene2$peak[ocr_gene2$peak_type=="Promoter"]))
length(unique(ocr_gene2$gene[ocr_gene2$peak_type=="Promoter"]))

length(unique(ocr_gene2$peak[ocr_gene2$peak_type=="enhancer"]))
length(unique(ocr_gene2$gene[ocr_gene2$peak_type=="enhancer"]))

# 1. Transform OCR-TF matrix to data frame
df <- as.data.frame(ocr_tf_matrix2)
df$ID <- rownames(df)
ocr_tf_df <- melt(df, id.vars = "ID")
ocr_tf_df <- ocr_tf_df[ocr_tf_df$value==1,]
names(ocr_tf_df)[2] <- "TF"
for(i in 1:nrow(ocr_tf_df)){
  ocr_tf_df$tf_ensembl[i] <- rna_anno$ensembl_gene_id[rna_anno$hgnc_symbol==ocr_tf_df$TF[i]]
}
ocr_tf_df$value <- NULL

# TF-OCR correlation
common_ids <- intersect(colnames(rna),colnames(atac))
rna2 <- rna[,common_ids]
atac2 <- atac[,common_ids]

for(i in 1:nrow(ocr_tf_df)){
  rna_target <- rna2[ocr_tf_df$tf_ensembl[i],] 
  atac_target <- atac2[ocr_tf_df$ID[i],]
  rr <- cor.test(as.numeric(atac_target), rna_target, method="spearman")
  ocr_tf_df$TF_OCR_spearman_pval[i] <- rr$p.value
  ocr_tf_df$TF_OCR_spearman_rho[i] <- rr$estimate
}

# 2. Merge TF-OCR and OCR-gene dataframes
ocr_tf_gene <- merge(ocr_tf_df, ocr_gene2, by.x="ID", by.y="peak")

# 3. Include gene-TF correlation
tmp <- unique(data.frame(gene=ocr_tf_gene$gene,tf_ensembl=ocr_tf_gene$tf_ensembl, tf=ocr_tf_gene$TF))
for(i in 1:nrow(tmp)){
  rr <- cor.test(rna[tmp$gene[i],],rna[tmp$tf_ensembl[i],], method="spearman")
  tmp$gene_tf_spearman_pval[i] <- rr$p.value
  tmp$gene_tf_spearman_rho[i] <- rr$estimate
}

ocr_tf_gene$gene_tf_pair <- paste0(ocr_tf_gene$TF, "_", ocr_tf_gene$gene)
tmp$gene_tf_pair <- paste0(tmp$tf,"_",tmp$gene)

ocr_tf_gene <- merge(ocr_tf_gene, tmp, by.x="gene_tf_pair", by.y="gene_tf_pair")


# 4. Get a summary table of: TF-OCR-Gene links:

ocr_tf_gene_summary_table <- data.frame(ocr = ocr_tf_gene$ID,
                                gene_ensembl = ocr_tf_gene$gene.x,
                                tf_ensembl = ocr_tf_gene$tf_ensembl.x,
                                ocr_gene_pval = ocr_tf_gene$spearman_pval,
                                ocr_gene_rho = ocr_tf_gene$spearman_rho,
                                ocr_gene_adj.pval = ocr_tf_gene$spearman_adj.pval,
                                ocr_tf_pval = ocr_tf_gene$TF_OCR_spearman_pval,
                                ocr_tf_rho = ocr_tf_gene$TF_OCR_spearman_rho,
                                gene_tf_pval = ocr_tf_gene$gene_tf_spearman_pval,
                                gene_tf_rho = ocr_tf_gene$gene_tf_spearman_rho,
                                gene_symbol = ocr_tf_gene$hgnc_symbol,
                                tf_symbol = ocr_tf_gene$TF,
                                peak_type = ocr_tf_gene$peak_type)



# 5. Table filtering:
# i) keep only TF-OCR-Genes from DEGs over b-cell diff
# ii) exclude enhancer OCRs not showing DA.

ocr_tf_gene_summary_table_filtered <- ocr_tf_gene_summary_table[ocr_tf_gene_summary_table$gene_ensembl%in%cluster_rna$gene,]
sel <- ocr_tf_gene_summary_table$ocr[ocr_tf_gene_summary_table$peak_type=="enhancer"][which(!ocr_tf_gene_summary_table$ocr[ocr_tf_gene_summary_table$peak_type=="enhancer"]%in%cluster_atac$peak)]
ocr_tf_gene_summary_table_filtered <- ocr_tf_gene_summary_table_filtered[!ocr_tf_gene_summary_table_filtered$ocr%in%sel,]

# 6. Identify patterns

ocr_tf_gene_summary_table_filtered$likelihood <- rep("yes",nrow(ocr_tf_gene_summary_table_filtered))
no_tf_change <- which(ocr_tf_gene_summary_table_filtered$likelihood=="yes" & !ocr_tf_gene_summary_table_filtered$tf_ensembl%in%cluster_rna$gene)
ocr_tf_gene_summary_table_filtered$likelihood[no_tf_change] <- "no_tf_change" #469,073

no_sig_gene_tf <- which(ocr_tf_gene_summary_table_filtered$likelihood=="yes" & ocr_tf_gene_summary_table_filtered$gene_tf_pval>0.05) 
ocr_tf_gene_summary_table_filtered$likelihood[no_sig_gene_tf] <- "no_tf-gene sig corr" #186,837

meaningless <- which(ocr_tf_gene_summary_table_filtered$likelihood=="yes" & ocr_tf_gene_summary_table_filtered$gene_tf_rho<0 & ocr_tf_gene_summary_table_filtered$ocr_tf_rho>0 & ocr_tf_gene_summary_table_filtered$ocr_gene_rho>0)
ocr_tf_gene_summary_table_filtered$likelihood[meaningless] <- "meaningless"
meaningless <- which(ocr_tf_gene_summary_table_filtered$likelihood=="yes" & ocr_tf_gene_summary_table_filtered$gene_tf_rho>0 & ocr_tf_gene_summary_table_filtered$ocr_tf_rho<0 & ocr_tf_gene_summary_table_filtered$ocr_gene_rho<0)
ocr_tf_gene_summary_table_filtered$likelihood[meaningless] <- "meaningless"
meaningless <- which(ocr_tf_gene_summary_table_filtered$likelihood=="yes" & ocr_tf_gene_summary_table_filtered$gene_tf_rho>0 & ocr_tf_gene_summary_table_filtered$ocr_tf_rho>0 & ocr_tf_gene_summary_table_filtered$ocr_gene_rho<0)
ocr_tf_gene_summary_table_filtered$likelihood[meaningless] <- "meaningless"
meaningless <- which(ocr_tf_gene_summary_table_filtered$likelihood=="yes" & ocr_tf_gene_summary_table_filtered$gene_tf_rho<0 & ocr_tf_gene_summary_table_filtered$ocr_tf_rho<0 & ocr_tf_gene_summary_table_filtered$ocr_gene_rho<0)
ocr_tf_gene_summary_table_filtered$likelihood[meaningless] <- "meaningless"
meaningless <- which(ocr_tf_gene_summary_table_filtered$likelihood=="yes" & ocr_tf_gene_summary_table_filtered$gene_tf_rho>0 & ocr_tf_gene_summary_table_filtered$ocr_tf_rho<0 & ocr_tf_gene_summary_table_filtered$ocr_gene_rho>0)
ocr_tf_gene_summary_table_filtered$likelihood[meaningless] <- "meaningless"
#117,846

repressor <- which(ocr_tf_gene_summary_table_filtered$likelihood=="yes" & ocr_tf_gene_summary_table_filtered$gene_tf_rho<0)
ocr_tf_gene_summary_table_filtered$likelihood[repressor] <- "repressor"

activator <- which(ocr_tf_gene_summary_table_filtered$likelihood=="yes" & ocr_tf_gene_summary_table_filtered$gene_tf_rho>0)
ocr_tf_gene_summary_table_filtered$likelihood[activator] <- "activator"

saveRDS(ocr_tf_gene_summary_table_filtered,"tf_ocr_gene_interactions.rds")
# 1,151,024 interactions
#-------------------------------------


## >> TF-OCR-gene data frame curation ------------------------------------------------------------------##
ocr_tf_gene_summary_table_filtered_meaningless <- ocr_tf_gene_summary_table_filtered[which(ocr_tf_gene_summary_table_filtered$likelihood=="activator" | ocr_tf_gene_summary_table_filtered$likelihood=="repressor"),]
#377,268 interactions

saveRDS(ocr_tf_gene_summary_table_filtered_meaningless,"tf_ocr_gene_interactions_curated.rds")
#-------------------------------------
