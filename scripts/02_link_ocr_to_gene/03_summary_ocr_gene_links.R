#**********************************************#
# LINK OCRs TO GENES                           #
# script: 03_summary_ocr_gene_links.R          #
#**********************************************#




# Required packages
require(RColorBrewer)
require(grDevices)
require(ggplot2)


##>>>>>>>>>>>>>>>>>>>>> Summary: GENE-PEAK-TYPE-SIG INFO

# Significant Enhancers
adj_peak_gene_association_intergenic_data_frame <- readRDS("distal_ocr_gene_corr.rds")

output_enhancer <- adj_peak_gene_association_intergenic_data_frame[adj_peak_gene_association_intergenic_data_frame$spearman_adj.pval<0.05,]
output_enhancer$peak_type <- rep("enhancer",nrow(output_enhancer))
#50,312 x 10
#Unique peaks: 22,123


# Promoters
output_promoter <- readRDS("promoter_ocr_gene_corr.rds")

output_promoter2 <- cbind(output_promoter[,c(1:6)],rep(NA,nrow(output_promoter)),output_promoter[,c(7:8)])
names(output_promoter2)[7] <- "dist"
output_promoter2$peak_type <- rep("Promoter",nrow(output_promoter2))
#12,247 x 10

final_output <- rbind(output_promoter2,output_enhancer)
#62,559 x 10


#Include gene info from RNA data
#Gene annotation
gene_anno <- readRDS("rna_gene_annotation.rds")
gene_anno <- gene_anno[gene_anno$ensembl_gene_id %in% final_output$gene,]

final_output2 <- merge(final_output,gene_anno,by.x="gene",by.y="ensembl_gene_id")

saveRDS(final_output2, paste0(format(Sys.time(), "%Y%m%d"),"_candidate_CREs.rds"))
#-------------------------------------
