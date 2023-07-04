#**********************************************#
# LINK OCRs TO GENES                           #
# script: 01_correlation_gene_ocr_promoter.R   #
#**********************************************#



##>>>>>>>>>>>>>>>>>>>>> Input data

#Normalized ATAC data
atac <- readRDS("atac_norm.rds")
atac_metadata <- readRDS("atac_metadata.rds")
atac_anno <- readRDS("atac_OCR_annotation.rds")

#Normalized RNA data
rna <- readRDS("rna_norm.rds")
rna_metadata <- readRDS("rna_metadata.rds")
rna_anno <- readRDS("rna_gene_annotation.rds")

#Common samples between omics
common_ids <- readRDS("rna_atac_paired_samples.rds")


rna_common <- as.matrix(rna[,common_ids])
atac_common <- as.matrix(log2(atac[,common_ids]+1))

#Gene annotation
rna_anno <- rna_anno[rownames(rna_common),]

#OCR genomic annotation
atac_anno <- atac_anno[rownames(atac_common),]
#-------------------------------------


##>>>>>>>>>>>>>>>>>>>>> Correlation

#---------------------------------------------------------------#
#********************* For promoter OCRs ***********************#
#---------------------------------------------------------------#
## Just check the correlation between each gene expression vs atac promoter (1kb).
# Total: 14,589

atac_anno_promoter <- atac_anno[which(atac_anno$annotation_simplified=="Promoter"),]
atac_promoter <- atac_common[rownames(atac_anno_promoter),] #promoter peaks info

peak_gene_association_promoter <- list()
for(i in 1:nrow(rna_common)){
 sel <- rownames(atac_anno_promoter)[which(atac_anno_promoter$ENSEMBL==rownames(rna_common)[i])]
 cat("gene ",i,"-peaks ",length(sel),"\n")
 cor_p <- vector()
 cor_r <- vector()
 cor2_p <- vector()
 cor2_r <- vector()
 cor_dist <- vector()
 if(length(sel)!=0){
 for(j in 1:length(sel)){
   cc <- cor.test(rna_common[i,],atac_promoter[sel[j],])
   cor_p <- c(cor_p,cc$p.value)
   cor_r <- c(cor_r,cc$estimate)
   cc2 <- cor.test(rna_common[i,],atac_promoter[sel[j],],method="spearman")
   cor2_p <- c(cor2_p,cc2$p.value)
   cor2_r <- c(cor2_r,cc2$estimate)
  }
 cor_output <- data.frame(peak=sel,pearson_pval=cor_p,pearson_rho=cor_r,spearman_pval=cor2_p,spearman_rho=cor2_r)
 }
 if(length(sel)==0){
 cor_output <- NA
 }
 peak_gene_association_promoter[[i]] <- cor_output
}
names(peak_gene_association_promoter) <- rownames(rna_common)


#determine the number of contrasts
n_contrasts_promoters <- unlist(lapply(peak_gene_association_promoter, FUN=function(X){if(is.data.frame(X)){return(nrow(X))}; if(is.na(X)){return(0)}}))
sum(n_contrasts_promoters) #12,247 contrasts 

a <- n_contrasts_promoters
summary(a)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.0000  0.0000  0.0000  0.5222  1.0000  3.0000


###Let's see the signifcance of the promoter and mrna
adj_peak_gene_association_promoter <- lapply(peak_gene_association_promoter,FUN=function(X){
  if(is.data.frame(X)){
    X$pearson_adj.pval <- p.adjust(X$pearson_pval,method="bonferroni", n=sum(n_contrasts_promoters))
    X$spearman_adj.pval <- p.adjust(X$spearman_pval,method="bonferroni", n=sum(n_contrasts_promoters))
  }
  return(X)
})

#Save results
adj_peak_gene_association_promoter_data_frame <- bind_rows(adj_peak_gene_association_promoter[unlist(lapply(adj_peak_gene_association_promoter, is.data.frame))], .id="gene")
saveRDS(adj_peak_gene_association_promoter_data_frame, paste0(format(Sys.time(), "%Y%m%d"),"_promoter_ocr_gene_corr.rds"))
