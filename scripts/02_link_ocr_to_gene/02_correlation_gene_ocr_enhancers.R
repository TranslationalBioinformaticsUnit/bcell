#**********************************************#
# LINK OCRs TO GENES                           #
# script: 02_correlation_gene_ocr_enhancer.R   #
#**********************************************#


##>>>>>>>>>>>>>>>>>>>>> Input data

#Normalized ATAC data
atac <- readRDS("atac_norm.rds")
atac_metadata <- readRDS("atac_metadata.rds")

#Normalized RNA data
rna <- readRDS("rna_norm.rds")
rna_metadata <- readRDS("rna_metadata.rds")

#Common samples between omics
ids_atac_rna <- readRDS("rna_atac_paired_samples.rds")
rna_common <- as.matrix(rna[,ids_atac_rna])
atac_common <- as.matrix(log2(atac[,ids_atac_rna]+1))

#Gene annotation
gene_anno <- readRDS("rna_gene_annotation.rds")
gene_anno <- gene_anno[rownames(rna),]

#OCR genomic annotation
ocr_anno <- readRDS("atac_OCR_annotation.rds")
ocr_anno <- ocr_anno[rownames(atac_common),]
#-------------------------------------


##>>>>>>>>>>>>>>>>>>>>> Correlation

#----------------------------------------------------------------#
#******************** Enhancer OCRs *****************************#
#----------------------------------------------------------------#
### Just check the correlation between each gene expression vs intergenic or intronic peaks 
## Total: 75,744 OCRs

ocr_anno <- ocr_anno[which(ocr_anno$annotation_simplified=="Distal Intergenic" | ocr_anno$annotation_simplified=="Intron" ),]
atac_enhancer <- atac_common[rownames(ocr_anno),] #intronic/distal intergenic peaks info
#Dim: 75,744 x 59 
atac_postion <- strsplit(rownames(atac_enhancer),"_")
atac_ini <- as.numeric(unlist(lapply(atac_postion,FUN=function(x){return(x[2])})))
names(atac_ini) <- rownames(atac_enhancer)
atac_fin <- as.numeric(unlist(lapply(atac_postion,FUN=function(x){return(x[3])})))
names(atac_fin) <- rownames(atac_enhancer)
atac_chr <- gsub("chr","",unlist(lapply(atac_postion,FUN=function(x){return(x[1])})))
names(atac_chr) <- rownames(atac_enhancer)


#Get the TSS and define the window
gene_anno$TSS <- rep(NA, nrow(gene_anno))
for(i in 1:nrow(gene_anno)){
  if(gene_anno$strand[i]==-1){ gene_anno$TSS[i] <- gene_anno$end_position[i]} 
  if(gene_anno$strand[i]==1){ gene_anno$TSS[i] <- gene_anno$start_position[i]}
}

##Define the window of 1Mb
window_up <- gene_anno$TSS-1000000
window_dw <- gene_anno$TSS+1000000


peak_gene_association_intergenic <- list()
for(i in 1:nrow(rna_common)){  #for each gene
 sel <- which(atac_chr==gene_anno$chromosome_name[i] & atac_ini>=window_up[i] & atac_fin<=window_dw[i]) #identify candidate ocrs
 sel_id <- names(sel) #keep instead the postion, the peak id
 cor_p <- vector()
 cor_r <- vector()
 cor2_p <- vector()
 cor2_r <- vector()
 cor_dist <- vector()
 if(length(sel)!=0){
 for(j in 1:length(sel)){
   cc <- cor.test(rna_common[i,],atac_enhancer[sel_id[j],])
   cor_p <- c(cor_p,cc$p.value)
   cor_r <- c(cor_r,cc$estimate)
   cc2 <- cor.test(rna_common[i,],atac_enhancer[sel_id[j],],method="spearman")
   cor2_p <- c(cor2_p,cc2$p.value)
   cor2_r <- c(cor2_r,cc2$estimate)
   mean_peak <- atac_ini[sel_id[j]]+((atac_fin[sel_id[j]]-atac_ini[sel_id[j]])/2)
   if(gene_anno$strand[i]==-1){cor_dist <- c(cor_dist,mean_peak-gene_anno$TSS[i])} 
   if(gene_anno$strand[i]==1){cor_dist <- c(cor_dist,gene_anno$TSS[i]-mean_peak)}
  }
 cor_output <- data.frame(peak=sel_id,pearson_pval=cor_p,pearson_rho=cor_r,spearman_pval=cor2_p,spearman_rho=cor2_r,dist=cor_dist)
 }
 if(length(sel)==0){
 cor_output <- NA
 }
 peak_gene_association_intergenic[[i]] <- cor_output
 cat(i,"-")
}


#determine the number of contrasts
n_contrasts_intergenic <- unlist(lapply(peak_gene_association_intergenic, FUN=function(X){if(is.data.frame(X)){return(nrow(X))}; if(is.na(X)){return(0)}}))

sum(n_contrasts_intergenic) #1,639,268 contrasts 

a <- n_contrasts_intergenic
summary(a)
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   0.00   51.00   70.00   69.89   87.00  188.00


###Let's see the signifcance of the enhancer and mrna
adj_peak_gene_association_intergenic <- lapply(peak_gene_association_intergenic,FUN=function(X){
  if(is.data.frame(X)){
    X$pearson_adj.pval <- p.adjust(X$pearson_pval,method="bonferroni", n=sum(n_contrasts_intergenic))
    X$spearman_adj.pval <- p.adjust(X$spearman_pval,method="bonferroni", n=sum(n_contrasts_intergenic))
  }
  return(X)
})
names(adj_peak_gene_association_intergenic) <- gene_anno$ensembl_gene_id


#Save results
adj_peak_gene_association_intergenic_data_frame <- bind_rows(adj_peak_gene_association_intergenic[unlist(lapply(adj_peak_gene_association_intergenic, is.data.frame))], .id="gene")
saveRDS(adj_peak_gene_association_intergenic_data_frame, paste0(format(Sys.time(), "%Y%m%d"),"_distal_ocr_gene_corr.rds"))
