#****************************************#
# OCR-TF analysis - with chromVAR tool   #
# script: 20_TF_chromVAR.R               #
#****************************************#


# 20. OCR - TF analysis - chromVAR


# Load required libraries
library(GenomicRanges)
library(SummarizedExperiment)
library(chromVAR)
library(TFBSTools)
library(ComplexHeatmap)
library(ggplot2)
library(RColorBrewer)
library(pals)
library(BuenColors)
library(circlize)
library(BSgenome.Hsapiens.UCSC.hg38)

# Load internal function
source("chromVAR_internal_functions.R")



##>>>>>>>>>>>>>>>>>>>>> Inputs required

# 1- Set of OCR
peaks_file="consensus_all_cell_types.bed"

# Read consensus peaks
peaks <- read.table(peaks_file) #105,530
rownames(peaks) <- paste(peaks[,1],peaks[,2],peaks[,3],sep="_")
colnames(peaks) <- c("chr","start","end")

# Create the GRanges of consensus peaks
peak_granges <- makeGRangesFromDataFrame(peaks)


# 2- OCR counts and metadata
# counts
counts_file <- "count_data_atac.txt"
my_counts_matrix <- as.matrix(read.delim(counts_file, check.names = FALSE))
#remove outlier sample
my_counts_matrix <- my_counts_matrix[,-which(colnames(my_counts_matrix)=="ImmatureB_ATAC_D215_S5")]

# metadata
atac_metadata <- read.table("metadata_atac.txt", sep="\t", header=TRUE)
rownames(atac_metadata) <- atac_metadata$SampleID
atac_metadata2 <- atac_metadata[colnames(my_counts_matrix),]
#order cell type labels
atac_metadata2$Tissue <- factor(as.character(atac_metadata2$Tissue),levels=c("HSC","CLP","proB","preB","ImmatureB","Transitional_B","Naive_CD5pos","Naive_CD5neg"))


# 3 - OCR-TF (motif) binary matrix
ocr_tf_matrix <- read.table("ocr_tf_matrix.txt", sep="\t", header=TRUE, check.names = FALSE)
#-------------------------------------



##>>>>>>>>>>>>>>>>>>>>> Run chromVAR with internal function chromVAR_assay (Figure 3B)
chromVAR_assay(ocr_granges=peak_granges, ocr_counts=my_counts_matrix, ocr_metadata=atac_metadata2, ocr_motif_anno=ocr_tf_matrix, assay_name="OCR_TF", target_tfs=c("PAX","GATA","Ebox","RUNX","EBF1","TAL","TCF","MYB","IKZF1","GFI"))
#-------------------------------------



# >>>>>>>>>>>> TFBS + RNA heatmap plot (Figure 3C)
# TFBS data
TF_scores <- read.table("TF_TFBS_score.txt", sep="\t", dec=".", header=TRUE)
TF_deviations <- read.table("TF_TFBS_deviations.txt", sep="\t", dec=".", header=TRUE)

variability <- read.table("TF_TFBS_variability.txt", sep="\t", dec=".", header=TRUE)

# RNA data
rna <- read.table("voom_normalized_filtered_data_rna.txt", sep="\t", header=TRUE, check.names = FALSE)
rna_metadata <- read.table("metadata_rna.txt", sep="\t", header=TRUE)
rna_metadata$CellType <- factor(as.character(rna_metadata$CellType),levels=c("HSC","CLP","proB","preB","ImmatureB","Transitional_B","Naive_CD5pos","Naive_CD5neg"))
rna <- rna[,as.character(rna_metadata$Alias)] #23,454 x 79
colnames(rna) <- paste(rna_metadata$Donor,rna_metadata$CellType,sep="_")

anno <- read.table("ensmbl_gene_annotation.txt", sep="\t", header=TRUE, check.names = FALSE)
anno_rna <- anno[anno$ensembl_gene_id %in% rownames(rna),]
rownames(anno_rna) <- anno_rna$ensembl_gene_id
anno_rna <- anno_rna[rownames(rna),]


# Correlation TFBS-RNA
rownames(atac_metadata2) <- paste(atac_metadata2$Factor,atac_metadata2$Tissue,sep="_")
colnames(TF_deviations) <- rownames(atac_metadata2)
colnames(TF_scores) <- rownames(atac_metadata2)

ids_atac_rna <- intersect(rownames(atac_metadata2), colnames(rna)) 
length(ids_atac_rna) #59 share samples

rna2 <- as.matrix(rna[,ids_atac_rna])
rownames(rna2) <- anno_rna$hgnc_symbol

atac2 <- as.matrix(TF_deviations[,ids_atac_rna])

TFs_atac_rna <- intersect(rownames(atac2), rownames(rna2)) 
length(TFs_atac_rna) #275 share TFs (10 complexes excluded)

atac3 <- atac2[TFs_atac_rna,]
rna3 <- rna2[TFs_atac_rna,]
atac_metadata3 <- atac_metadata2[ids_atac_rna,]

#data matrix to plot
tfbs_to_plot <- TF_scores[TFs_atac_rna,]
rna_to_plot <- as.matrix(rna)
rownames(rna_to_plot) <- anno_rna$hgnc_symbol
rna_to_plot <- rna_to_plot[TFs_atac_rna,]


cor_p <- vector()
cor_r <- vector()
cor2_p <- vector()
cor2_r <- vector()

for(i in 1:nrow(rna3)){
  cc <- cor.test(rna3[i,],atac3[i,])
  cor_p <- c(cor_p,cc$p.value)
  cor_r <- c(cor_r,cc$estimate)
  cc2 <- cor.test(rna3[i,],atac3[i,],method="spearman")
  cor2_p <- c(cor2_p,cc2$p.value)
  cor2_r <- c(cor2_r,cc2$estimate)
}

TF_gene_association <- data.frame(TF=rownames(rna3),pearson_pval=cor_p,pearson_rho=cor_r,spearman_pval=cor2_p,spearman_rho=cor2_r)
rownames(TF_gene_association) <- TF_gene_association$TF
TF_gene_association$spearman_pval_adj <- p.adjust(TF_gene_association$spearman_pval,method="bonferroni")
write.table(TF_gene_association,paste(format(Sys.Date(),"%Y%m%d"),"_TF_gene_association.txt",sep=""),sep="\t",dec=".",quote=FALSE)

## All TFs plot
bcell_colors <- jdb_color_map(c("MPP","LMPP","CD8","CD4","GMP-B","Ery","mono","MEP"))
names(bcell_colors) <- levels(atac_metadata2$Tissue)

TFBS_RNA_heatmap(TF_gene_association=TF_gene_association, 
                 atac_metadata=atac_metadata2, 
                 rna_metadata=rna_metadata, 
                 color_group=bcell_colors, 
                 atac_group="Tissue",
                 rna_group="CellType", 
                 tfbs_data=tfbs_to_plot, 
                 rna_data=rna_to_plot, 
                 add_raw_rna=FALSE,
                 tfbs_variability=variability,
                 rna_sd=apply(rna_to_plot,1,sd),
                 target_tfs=c("PAX5","EBF1","GATA1","RUNX1","TAL1","TCF3","MYB","GFI1","ERG","ETV6","SPIB","ATF3","NFE2","IKZF1","ETS2"),
                 assay_name="275_TFs")


##Significant TFs (TFBS-RNA) (Figure 3C)
sel_sig <- rownames(TF_gene_association)[TF_gene_association$spearman_pval_adj<0.05] #83

TFBS_RNA_heatmap(TF_gene_association=TF_gene_association[sel_sig,], 
                 atac_metadata=atac_metadata2, 
                 rna_metadata=rna_metadata, 
                 color_group=bcell_colors, 
                 atac_group="Tissue",
                 rna_group="CellType", 
                 tfbs_data=tfbs_to_plot[sel_sig,], 
                 rna_data=rna_to_plot[sel_sig,],
                 add_raw_rna=FALSE,
                 tfbs_variability=variability[sel_sig,],
                 rna_sd=apply(rna_to_plot,1,sd)[sel_sig],
                 target_tfs=c("PAX5","EBF1","GATA1","RUNX1","TAL1","TCF3","MYB","GFI1","ERG","ETS2","ETV6","SPIB","ATF3","NFE2","IKZF1"),
                 assay_name="83_bonferroni_TFs")
#-------------------------------------



# >>>>>>>>>>>> Heatmap of sorted TF for each cell-type based on TFBS deviation (Figure 3D)
TF_deviations_sig <- TF_deviations[sel_sig,]

order_TFBS <- matrix(ncol=length(cell_types), nrow=nrow(TF_deviations_sig))
rownames(order_TFBS) <- rownames(TF_deviations_sig)

for(i in 1:length(cell_types)){
  mean_tfbs <- rowMeans(TF_deviations_sig[,atac_metadata2$Tissue==cell_types[i]])
  oo <- order(mean_tfbs, decreasing = FALSE)
  for(j in 1:length(mean_tfbs)){
    order_TFBS[oo[j],i] <- j
  }
}

ha = HeatmapAnnotation(stage = cell_types,
                       col = list(stage =bcell_colors))

pdf(paste(format(Sys.Date(),"%Y%m%d"),"_sorted_sig_TFBS_heatmap.pdf",sep=""))
draw(Heatmap(order_TFBS, top_annotation = ha, 
             show_row_names = TRUE, show_column_names =FALSE, 
             col=jdb_palette("flame_light"), cluster_columns = FALSE,
             name="Position", row_names_gp = gpar(fontsize = 6),
             column_title = "TFBS accessibility score"))
dev.off()
#-------------------------------------
