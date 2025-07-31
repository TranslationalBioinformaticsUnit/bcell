#***************************************************
## Paper: "Uncovering the regulatory landscape of early human B-cell lymphopoiesis and its implications in the pathogenesis of B-ALL"
## Authors: Planell N et al.
## Date: 2025
## Code: Supplementary Figures (S3-S14) & Tables (S6, S8-S9)
## Input data available at: https://osf.io/gswpy/?view_only=39cfd50d5a7e44e5a2d141869dcb1315
#***************************************************

# Define the directory were input data from OSF has been downloaded
data_wd <- setwd("../GitHub_code/")

### Supplementary Figures & Tables ###
# Figure S2 -----------------------------------------------------------------
# RNA-seq: DEGs over differentiation ------------------------------------ ###
# >>> Required packages
library(ComplexHeatmap)
library(circlize)

# >>> Input data
# Gene Expression data
rna <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/rna_norm_data.rds"))
rna_metadata <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/rna_metadata.rds"))

# DEA by transitions
deg <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/differential_expression_analysis/rna_deg_clustering.rds"))

# >>> Plotting
# Plot colors
bcell_colors <- readRDS(paste0(data_wd,"/osfstorage-archive/00_plot_parameters/bcell_colors.rds"))
rna_colors <- readRDS(paste0(data_wd,"/osfstorage-archive/00_plot_parameters/rna_expression_colors.rds"))

# Plot
data_to_plot <- t(scale(t(rna[rownames(deg),]),center=TRUE, scale=TRUE))

exprs_pattern <- Heatmap(data_to_plot, 
                         name = "Gene expression (z-score)",
                         top_annotation = HeatmapAnnotation(stage = rna_metadata$CellType,
                                                            col = list(stage = bcell_colors)), 
                         column_order=order(rna_metadata$CellType), 
                         col=colorRamp2(seq(-3,3, length.out=length(rna_colors)),rna_colors),
                         show_column_names = FALSE, 
                         show_row_names = FALSE, 
                         cluster_rows = FALSE,  
                         cluster_columns = FALSE, 
                         row_order = order(deg$pattern_id),
                         use_raster=TRUE, 
                         width = unit(4, "cm"))
### --------------------------------------------------------------------- ###

# ATAC-seq: OCR Genomic location annotation ----------------------------- ###
# >>> Required packages
library(ggplot2)

# >>> Input data
# Chromatin accessibility data
atac_anno <- readRDS(paste0(data_wd,"/osfstorage-archive/01_atac_seq_data/atac_OCR_annotation.rds"))

# >>> Plotting
# Plot
data_to_plot <- data.frame(peak=unique(atac_anno$annotation_simplified), value=table(atac_anno$annotation_simplified))
p <- ggplot(data_to_plot, aes(fill=value.Var1, y=value.Freq, x=1)) + 
  geom_bar(position="stack", stat="identity")+
  scale_fill_viridis(discrete = T)

ggsave(plot=p, dpi=100, filename = paste0(format(Sys.time(), "%Y%m%d"),"_OCR_annotation_distribution.pdf"))
### --------------------------------------------------------------------- ###

# ATAC-seq: DARs over differentiation ----------------------------------- ###
# >>> Required packages
library(ComplexHeatmap)
library(circlize)

# >>> Input data
# Chromatin accessibility data
atac <- readRDS(paste0(data_wd,"/osfstorage-archive/01_atac_seq_data/atac_norm_data.rds"))
atac_metadata <- readRDS(paste0(data_wd,"/osfstorage-archive/01_atac_seq_data/atac_metadata.rds"))

# DARs by transitions
dar <- readRDS(paste0(data_wd,"/osfstorage-archive/01_atac_seq_data/differential_binding_analysis/atac_dar_clustering.rds"))

# >>> Plotting
# Plot colors
bcell_colors <- readRDS(paste0(data_wd,"/osfstorage-archive/00_plot_parameters/bcell_colors.rds"))
atac_colors <- readRDS(paste0(data_wd,"/osfstorage-archive/00_plot_parameters/atac_colors.rds"))

# Plot
data_to_plot <- t(scale(t(atac[rownames(dar),]),center=TRUE, scale=TRUE))

atac_pattern <- Heatmap(data_to_plot, 
                         name = "Gene expression (z-score)",
                         top_annotation = HeatmapAnnotation(stage = atac_metadata$CellType,
                                                            col = list(stage = bcell_colors)), 
                         column_order=order(atac_metadata$CellType), 
                         col=colorRamp2(seq(-3,3, length.out=length(atac_colors)),atac_colors),
                         show_column_names = FALSE, 
                         show_row_names = FALSE, 
                         cluster_rows = FALSE,  
                         cluster_columns = FALSE, 
                         row_order = order(dar$pattern_id),
                         use_raster=TRUE, 
                         width = unit(4, "cm"))
### --------------------------------------------------------------------- ###

# ATAC-seq: TFBS score -------------------------------------------------- ###
# >>> Required packages
library(ComplexHeatmap)
library(circlize)
library(BuenColors)

# >>> Input data
# TFBS data
tfbs <- read.table(paste0(data_wd,"/osfstorage-archive/03_TFs_OCR_links/chromVAR_tf/TFBS_scores_from_ocr_tf_matrix.txt"))
atac_metadata <- readRDS(paste0(data_wd,"/osfstorage-archive/01_atac_seq_data/atac_metadata.rds"))

# >>> Plotting
# Plot colors
bcell_colors <- readRDS(paste0(data_wd,"/osfstorage-archive/00_plot_parameters/bcell_colors.rds"))

# Plot
tfbs_plot <- Heatmap(tfbs,
                     top_annotation = HeatmapAnnotation(Subpopulation = atac_metadata$CellType, 
                                                        col = list(Subpopulation = bcell_colors),
                                                        show_annotation_name = FALSE),
                     col = jdb_palette("solar_extra"),
                     column_order=order(atac_metadata$CellType),
                     show_column_names = FALSE, 
                     show_row_names = FALSE,
                     cluster_rows = TRUE, 
                     show_row_dend = FALSE,  
                     cluster_columns = FALSE, 
                     use_raster=TRUE, 
                     width = unit(4, "cm"),
                     name="TFBS score"
)
### --------------------------------------------------------------------- ###

# cis-Regulatory Elements: Promoter OCR- gene --------------------------- ###
# >>> Required packages
library(ggplot2)

# >>> Input data
# Chromatin accessibility data
atac_anno <- readRDS(paste0(data_wd,"/osfstorage-archive/01_atac_seq_data/atac_OCR_annotation.rds"))

# >>> Plotting
promoter_ocr_gene <- table(atac_anno$ENSEMBL[atac_anno$annotation_simplified=="Promoter"])

data_to_plot <- data.frame(n_genes=c(1:3), 
                           value=c(sum(promoter_ocr_gene==1),
                                   sum(promoter_ocr_gene==2),
                                   sum(promoter_ocr_gene==3)))

axis_size <- element_text(size=25)
p <- ggplot(data=data_to_plot, aes(x=n_genes, y=value)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=value), vjust=-0.3, size=10)+
  theme_classic()+
  theme(axis.text=axis_size, 
        axis.title.x=axis_size,
        axis.title.y=axis_size,
        plot.title=axis_size) +
  ylim(0,max(data_to_plot$value)+400) +
  labs(title="Promoter-Gene associations", 
       x="Number of promoter OCRs", y = "Number of genes")

ggsave(plot = p, filename = paste0(format(Sys.time(), "%Y%m%d"),"_barplot_OCR_promoter_gene_annotation_1kb.pdf"), width=10,height=10)
### --------------------------------------------------------------------- ###

# cis-Regulatory Elements: Enhancer OCR- gene --------------------------- ###
# >>> Required packages
library(ggplot2)
library(ggbreak)

# >>> Input data
# CREs
ocr_gene <- readRDS(paste0(data_wd,"/osfstorage-archive/02_cCREs_OCR_gene_links/candidate_CREs.rds"))

# >>> Plotting
enhancer_ocr_gene <- table(ocr_gene$gene[ocr_gene$peak_type=="enhancer"])

data_to_plot <- data.frame(peaks=factor(c("1","2-5","6-10","11-20",">20"), levels=c("1","2-5","6-10","11-20",">20")),
                           value=c(sum(enhancer_ocr_gene==1),
                                   sum(enhancer_ocr_gene>1 & enhancer_ocr_gene<6, na.rm=TRUE),
                                   sum(enhancer_ocr_gene>5 & enhancer_ocr_gene<11, na.rm=TRUE),
                                   sum(enhancer_ocr_gene>10 & enhancer_ocr_gene<21),
                                   sum(enhancer_ocr_gene>20)))


axis_size <- element_text(size=25)
# Basic barplot
p <- ggplot(data=data_to_plot, aes(x=peaks, y=value)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=value), vjust=-0.3, size=10)+
  theme_classic() +
  theme(axis.text=axis_size, 
        axis.title.x=axis_size,
        axis.title.y=axis_size,
        plot.title=axis_size) +
  ylim(0,max(data_to_plot$value)+400) +
  labs(title="Signifincat correlation (Spearman)", 
       x="Number of enhancer OCRs", y = "Number of genes")

ggsave(plot = p, filename = paste0(format(Sys.time(), "%Y%m%d"),"_barplot_OCR_enhancer_gene_sig_association_spearman.pdf"), width=10, height=10)
###-----------------------------------------------------------------------------------###
# ------------------------------------------------------------------------------------ #



# Figure S3 -----------------------------------------------------------------
# TF activity profile across B-cell differentiation. 
# Heatmap representation of the TFBS activity score profiling for all the TF identified (n=275). 
# Samples are shown by columns (n=78) and grouped by cell subpopulation, and TFs by rows. 

# >>> Required packages
library(BuenColors)
library(ComplexHeatmap)

# >>> Input data
# TFBS data
tfbs <- read.table(paste0(data_wd,"/osfstorage-archive/03_TFs_OCR_links/chromVAR_tf/TFBS_scores_from_ocr_tf_matrix.txt"))
atac_metadata <- readRDS(paste0(data_wd,"/osfstorage-archive/01_atac_seq_data/atac_metadata.rds"))

# >>> Plotting
# Plot colors
bcell_colors <- readRDS(paste0(data_wd,"/osfstorage-archive/00_plot_parameters/bcell_colors.rds"))

# Plot
tfbs_plot <- Heatmap(tfbs,
                     top_annotation = HeatmapAnnotation(Subpopulation = atac_metadata$CellType, 
                                                        col = list(Subpopulation = bcell_colors),
                                                        show_annotation_name = FALSE),
                     col = jdb_palette("solar_extra"),
                     cluster_columns = TRUE,
                     cluster_rows = FALSE,
                     column_gap = unit(0.5, "mm"),
                     show_column_names = FALSE,
                     width = unit(4, "cm"),
                     height = unit(1, "mm")*nrow(tfbs),
                     row_names_gp = gpar(fontsize = 4),
                     row_title=NULL,
                     column_title = NULL,
                     name="TFBS score"
)

pdf(paste0(format(Sys.time(), "%Y%m%d"),"_heatmap_tfbs_zscore_clustered.pdf"), height = 14)
tfbs_plot
dev.off()
# ------------------------------------------------------------------------------------ #




# Figure S4 - Panel A -----------------------------------------------------------------
# B-cell subpopulations signatures
# Heatmap representation of 8,642 genes defining the signature profile of the B-cell subpopulations. 
# Biomarker genes for each B-cell subpopulation were derived from 
# significant genes (logFC≥1.5 and adjusted p-value≤0.05) in at least 4 contrasts (see methods).

# >>> Required packages
library(ComplexHeatmap)
library(BuenColors)

# >>> Input data
# Gene Expression data
rna <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/rna_norm_data.rds"))
rna_metadata <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/rna_metadata.rds"))
rna_anno <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/rna_gene_annotation.rds"))
# Human B-cell subpopulations biomarkers
bcell_biomarkers <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/differential_expression_analysis/Bcell_biomarkers.rds"))

# >>> Plotting
# Plot colors
bcell_colors <- readRDS(paste0(data_wd,"/osfstorage-archive/00_plot_parameters/bcell_colors.rds"))

# Plot
rna_biomarkers_heatmap <- Heatmap(t(scale(t(rna[unique(unlist(bcell_biomarkers)),]))),
                                  top_annotation = HeatmapAnnotation(Subpopulation = rna_metadata$CellType, 
                                                                     col = list(Subpopulation = bcell_colors),
                                                                     show_annotation_name = FALSE),
                                  cluster_rows = TRUE,
                                  column_order = order(rna_metadata$CellType),
                                  show_column_names = FALSE,
                                  show_row_names = FALSE,
                                  width = unit(4, "cm"),
                                  row_names_gp = gpar(fontsize = 4),
                                  show_row_dend = FALSE,
                                  row_title=NULL,
                                  column_title = NULL,
                                  name="Gene expression"
)


pdf(paste0(format(Sys.time(), "%Y%m%d"),"_heatmap_bcell_biomarkers.pdf"))
rna_biomarkers_heatmap
dev.off()
# ------------------------------------------------------------------------------------ #



# Figure S4 - Panel B -----------------------------------------------------------------
# B-cell signature upset plot
# Upset plot showing the overlap between the gene expression profile of each B-cell subpopulation.

# >>> Required packages
library(ComplexHeatmap)

# >>> Input data
# Human B-cell subpopulations biomarkers
bcell_biomarkers <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/differential_expression_analysis/Bcell_biomarkers.rds"))

# >>> Plotting
#UpSet plot
m = make_comb_mat(Bcell_4sig_biomarkers)
p <- UpSet(m, comb_order = order(comb_size(m), decreasing = TRUE))

pdf(paste0(format(Sys.time(), "%Y%m%d"),"_upset_bcell_biomarkers.pdf"))
print(p)
dev.off()
# ------------------------------------------------------------------------------------ #



# Figure S5 - Panel A -----------------------------------------------------------------
# PCA RNA-Seq mice
# Principal Component Analysis (PCA) of gene expression from public mice B-cell 
# subpopulations (17 samples) showing the distribution of cell subpopulations. 

# >>> Required packages
library(ggplot2)
library(BuenColors)

# >>> Input data
# Gene Expression data from mice (GSE109125)
rna_mice <- readRDS(paste0(data_wd,"/osfstorage-archive/06_mice/rna_norm_data_mice.rds"))
rna_metadata_mice <- readRDS(paste0(data_wd,"/osfstorage-archive/06_mice/rna_metadata_mice.rds"))

# >>> Computing
# PCA computation
pca_output <- prcomp(t(rna_mice))
var_prcomp <- pca_output$sdev^2

# PCs variance
pcvar <- data.frame(var=var_prcomp/sum(var_prcomp), 
                    pc=c(1:length(var_prcomp)))

# >>> Plotting
# Plot colors
bcell_colors_mice <- readRDS(paste0(data_wd,"/osfstorage-archive/00_plot_parameters/bcell_colors_mice.rds"))

# Plot
data_to_plot_pca <- data.frame(PC1=pca_output$x[,1],PC2=pca_output$x[,2], CellType=rna_metadata_mice$Bcell_subtype)

pca_mice <- ggplot(data_to_plot_pca, aes(x = PC1, y = PC2, color = CellType)) +
  geom_point(size = 3, show.legend=FALSE) + 
  coord_fixed(ratio=1.35) + 
  scale_color_manual(values = bcell_colors_mice) + 
  directlabels::geom_dl(aes(label = CellType), method = list("smart.grid", cex=1)) + 
  labs(title="PCA: RNA - Mice", 
       x = paste("PC1 (",round(pcvar[1,1]*100,0),"%)", sep=""), 
       y = paste("PC2 (",round(pcvar[2,1]*100,0),"%)", sep=""))+
  theme_classic()

pdf(paste(format(Sys.time(), "%Y%m%d"),"_pca_mice.pdf",sep=""))
pca_mice
dev.off()
# ------------------------------------------------------------------------------------ #



# Figure S5 - Panel B -----------------------------------------------------------------
# Heatmap representation of gene expression profile for 54 well-known genes 
# associated with mice B-cell lymphopoiesis. 

# >>> Required packages
library(ComplexHeatmap)
library(BuenColors)

# >>> Input data
# Gene Expression data from mice (GSE109125)
rna_mice <- readRDS(paste0(data_wd,"/osfstorage-archive/06_mice/rna_norm_data_mice.rds"))
rna_metadata_mice <- readRDS(paste0(data_wd,"/osfstorage-archive/06_mice/rna_metadata_mice.rds"))

# >>> Ploting
# Plot colors
bcell_colors_mice <- readRDS(paste0(data_wd,"/osfstorage-archive/00_plot_parameters/bcell_colors_mice.rds"))

# Define biomarkers to plot
key_bcell_genes <- data.frame(gene=c("Crhbp","Mllt3","Hlf",
                                     "Cd34","Prom1","Kit","Smim24","C1qtnf4","Flt3",
                                     "Cd38","Cd19","Ms4a1","Cd5","Cd24a",
                                     "Rag1","Rag2","Dntt",
                                     "Vpreb1","Igll1","Cd79a","Cd79b","Ighm","Ighd",
                                     "Hmgb1","Hmgb2","Hmgb3","Mki67","Top2a","Cks2","Nasp","Pttg1","Stmn1","Smc4","Nusap1",
                                     "Mcl1","Bcl2","Bcl2l1","Socs1","Akt1","Tnfrsf13c",
                                     "Cxcr4","Il7r",
                                     "Lef1","Ikzf3","Irf4","Irf8","Foxo1","Ebf1","Tcf3","Pax5","Myb","Ikzf1","Spi1","Sox4"),
                              class=factor(c(rep("Stemness",3),
                                             rep("Lymphoid Progenitors",6),
                                             rep("B-cell lineage",5),                                      
                                             rep("Ig-recombination",3),
                                             rep("pre-BCR/BCR",6),
                                             rep("Proliferation",11),
                                             rep("Survival",6),
                                             rep("Differentiation",2),
                                             rep("Transcription factors",12)),levels=c("Stemness","Lymphoid Progenitors","B-cell lineage",
                                                                                       "Ig-recombination","pre-BCR/BCR",
                                                                                       "Differentiation","Proliferation","Survival","Transcription factors"))
)

rownames(key_bcell_genes) <- key_bcell_genes$gene

# Plot
data_to_plot <- rna_mice[key_bcell_genes$gene,]

mice_exprs_pattern <- Heatmap(data_to_plot, 
                              top_annotation = HeatmapAnnotation(Subpopulation = rna_metadata_mice$Bcell_subtype,
                                                                 col = list(Subpopulation = bcell_colors_mice),
                                                                 show_annotation_name = FALSE), 
                              row_title_gp = gpar(fontsize = 7),
                              cluster_columns = FALSE,
                              column_split = rna_metadata_mice$Bcell_subtype,
                              col=jdb_palette("wolfgang_extra"),
                              row_split = key_bcell_genes$class,
                              column_gap = unit(0.5, "mm"),
                              show_row_dend = FALSE,
                              show_column_names = FALSE, 
                              show_row_names = TRUE, 
                              cluster_rows = TRUE,
                              cluster_row_slices = FALSE,
                              use_raster=TRUE, 
                              row_title_rot = 0,
                              width = unit(4, "cm"),
                              height = unit(3, "mm")*nrow(data_to_plot),
                              row_names_gp = gpar(fontsize = 6),
                              column_title = NULL,
                              name="Gene expression")


pdf(paste0(format(Sys.time(), "%Y%m%d"),"_bcell_mice_expression_biomarkers.pdf"), height = 10)
mice_exprs_pattern
dev.off()
# ------------------------------------------------------------------------------------ #

# Figure S6 - Panel A -----------------------------------------------------------------
# Heatmap representation of genes included in the GRN (n=7,310 genes).

# >>> Required packages
library(ComplexHeatmap)
library(circlize)
library(magick)

# >>> Input data
# Gene Expression data
rna <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/rna_norm_data.rds"))
rna_metadata <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/rna_metadata.rds"))
rna_anno <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/rna_gene_annotation.rds"))
# DEGs clusters
deg_clust <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/differential_expression_analysis/rna_deg_clustering.rds"))
# GRN data
tf_ocr_gene <- read.table(paste0(data_wd,"/osfstorage-archive/04_gene_regulatory_networks/tf_ocr_gene_interactions_curated.txt"),
                          sep="\t",
                          dec=",",
                          header=TRUE)

# >>> Plotting
data_to_plot <- t(scale(t(rna[unique(tf_ocr_gene$gene_ensembl),]),center=TRUE, scale=TRUE))
data_to_plot2 <- deg_clust[rownames(data_to_plot),]
row_order <- order(data_to_plot2$pattern_id)

data_to_plot <- data_to_plot[row_order,]
data_to_plot2 <- data_to_plot2[row_order,]

# Genes to highlight
target_genes <- c("PAX5","HOXB6","BCL2","NOTCH1","CDH2","MEIS1","CCND1",
                  "EBF1","TAL1","GATA1","TCF3","CD19","MS4A1","MME","CD79A","CD79B",
                  "CD38","CD34","CD5","IGHM","DNTT","HMGB2","EZH2","HMGB1","MYBL2","BIRC5",
                  "CD83","BTG1")

sel_pos <- vector()
target_genes2 <- vector()
for(i in 1:length(target_genes)){
  ss <- which(rownames(data_to_plot)==rna_anno$ensembl_gene_id[rna_anno$hgnc_symbol==target_genes[i]])
  if(length(ss)>0){
    sel_pos <- c(sel_pos,ss)
    target_genes2 <- c(target_genes2,target_genes[i])
  }
}


# Plot colors
bcell_colors <- readRDS(paste0(data_wd,"/osfstorage-archive/00_plot_parameters/bcell_colors.rds"))
rna_colors <- readRDS(paste0(data_wd,"/osfstorage-archive/00_plot_parameters/rna_expression_colors.rds"))

# Plot
exprs_pattern <- Heatmap(data_to_plot, 
                         name = "Expression patterns",
                         top_annotation = HeatmapAnnotation(stage = rna_metadata$CellType,
                                                            col = list(stage = bcell_colors)), 
                         column_order=order(rna_metadata$CellType), 
                         col=colorRamp2(seq(-3,3, length.out=length(rna_colors)),rna_colors),
                         show_column_names = FALSE, 
                         show_row_names = FALSE, 
                         row_names_gp = gpar(fontsize = 3), 
                         cluster_rows = FALSE,
                         row_split=data_to_plot2$pattern_id,
                         cluster_columns = FALSE, 
                         row_gap = unit(0.5, "mm"), 
                         use_raster=TRUE, 
                         width = unit(4, "cm"))

binary_pattern <- Heatmap(data_to_plot2[,5:12], 
                          name = "Peaks vs Cell Type",
                          col = colorRamp2(c(0, 1), c("beige","darkgoldenrod")),
                          top_annotation = HeatmapAnnotation(stage = labels(bcell_colors),
                                                             col = list(stage = bcell_colors)),
                          right_annotation = rowAnnotation(genes=anno_mark(at = sel_pos, 
                                                                           labels = target_genes2)), 
                          show_column_names = TRUE, 
                          show_row_names = FALSE, 
                          row_names_gp = gpar(fontsize = 3), 
                          cluster_rows = FALSE, 
                          row_split=data_to_plot2$pattern_id,
                          cluster_columns = FALSE, 
                          row_gap = unit(0.5, "mm"), 
                          column_gap = unit(1, "mm"), 
                          use_raster=TRUE, 
                          border="grey", 
                          column_split=factor(labels(bcell_colors), levels=labels(bcell_colors)), 
                          width = unit(3, "cm"))

pdf(paste0(format(Sys.time(), "%Y%m%d"),"_heatmap_genes_grn.pdf"), width=10)
exprs_pattern + binary_pattern
dev.off()
# ------------------------------------------------------------------------------------ #



# Figure S6 - Panel B -----------------------------------------------------------------
# Heatmap representation of OCRs included in the GRN (n=16,074 OCRs).

# >>> Required packages
library(ComplexHeatmap)
library(circlize)
library(magick)

# >>> Input data
# Chromatin accessibility data
atac <- readRDS(paste0(data_wd,"/osfstorage-archive/01_atac_seq_data/atac_norm_data.rds"))
atac_metadata <- readRDS(paste0(data_wd,"/osfstorage-archive/01_atac_seq_data/atac_metadata.rds"))
atac_anno <- readRDS(paste0(data_wd,"/osfstorage-archive/01_atac_seq_data/atac_OCR_annotation.rds"))
# DARs clusters
dar_clust <- readRDS(paste0(data_wd,"/osfstorage-archive/01_atac_seq_data/differential_binding_analysis/atac_dar_clustering.rds"))
# GRN data
tf_ocr_gene <- read.table(paste0(data_wd,"/osfstorage-archive/04_gene_regulatory_networks/tf_ocr_gene_interactions_curated.txt"),
                          sep="\t",
                          dec=",",
                          header=TRUE)

# >>> Plotting
data_to_plot <- t(scale(t(atac[unique(tf_ocr_gene$ocr),]),center=TRUE, scale=TRUE))
data_to_plot2 <- dar_clust[rownames(data_to_plot),]
row_order <- order(data_to_plot2$pattern_id)

data_to_plot <- data_to_plot[row_order,]
data_to_plot2 <- data_to_plot2[row_order,]
atac_anno2 <- atac_anno[rownames(data_to_plot),]

# Genes to highlight
target_genes <- c("PAX5","HOXB6","BCL2","NOTCH1","CDH2","MEIS1","CCND1",
                  "EBF1","TAL1","GATA1","TCF3","CD19","MS4A1","MME","CD79A","CD79B",
                  "CD38","CD34","CD5","IGHM","DNTT","HMGB2","EZH2","HMGB1","MYBL2","BIRC5",
                  "CD83","BTG1")

sel_pos <- vector()
target_genes2 <- vector()
for(i in 1:length(target_genes)){
  ss <- which(rownames(data_to_plot)%in%rownames(atac_anno2)[which(atac_anno2$annotation_simplified=="Promoter" & atac_anno2$SYMBOL==target_genes[i])])
  if(length(ss)>0){
    sel_pos <- c(sel_pos,ss)
    target_genes2 <- c(target_genes2,paste0(target_genes[i],"_",rownames(data_to_plot)[ss]))
  }
}

# Plot colors
bcell_colors <- readRDS(paste0(data_wd,"/osfstorage-archive/00_plot_parameters/bcell_colors.rds"))
rna_colors <- readRDS(paste0(data_wd,"/osfstorage-archive/00_plot_parameters/rna_expression_colors.rds"))
atac_colors <- readRDS(paste0(data_wd,"/osfstorage-archive/00_plot_parameters/atac_colors.rds"))
peak_color <- c("Promoter"="#F6313E","Intron"="#46A040","Distal Intergenic"="#490C65")

# Plot
accessibility_pattern <- Heatmap(data_to_plot, 
                                 name = "Accessibility pattern",
                                 top_annotation = HeatmapAnnotation(stage = atac_metadata$CellType,
                                                                    col = list(stage = bcell_colors)), 
                                 column_order=order(atac_metadata$CellType),
                                 right_annotation = rowAnnotation(peak = atac_anno2$annotation_simplified,
                                                                  col = list(peak = peak_color)),		
                                 col=colorRamp2(seq(-3,3, length.out=length(atac_colors)),atac_colors),
                                 show_column_names = FALSE, 
                                 show_row_names = FALSE, 
                                 row_names_gp = gpar(fontsize = 3), 
                                 cluster_rows = FALSE,
                                 row_split=data_to_plot2$pattern_id,
                                 cluster_columns = FALSE, 
                                 row_gap = unit(0.5, "mm"), 
                                 use_raster=TRUE, 
                                 width = unit(4, "cm"))

binary_pattern <- Heatmap(data_to_plot2[,5:12], 
                          name = "Peaks vs Cell Type",
                          col = colorRamp2(c(0, 1), c("beige","darkgoldenrod")),
                          top_annotation = HeatmapAnnotation(stage = labels(bcell_colors),
                                                             col = list(stage = bcell_colors)),
                          right_annotation = rowAnnotation(genes=anno_mark(at = sel_pos, 
                                                                           labels = target_genes2)), 
                          show_column_names = TRUE, 
                          show_row_names = FALSE, 
                          row_names_gp = gpar(fontsize = 3), 
                          cluster_rows = FALSE, 
                          row_split=data_to_plot2$pattern_id,
                          cluster_columns = FALSE, 
                          row_gap = unit(0.5, "mm"), 
                          column_gap = unit(1, "mm"), 
                          use_raster=TRUE, 
                          border="grey", 
                          column_split=factor(labels(bcell_colors), 
                                              levels=labels(bcell_colors)), 
                          width = unit(3, "cm"))

pdf(paste0(format(Sys.time(), "%Y%m%d"),"_heatmap_ocr_grn.pdf"), width=10)
accessibility_pattern + binary_pattern
dev.off()
# ------------------------------------------------------------------------------------ #



# Figure S6 - Panel C -----------------------------------------------------------------
## Heatmap representation of the TFBS activity score of those TFs included the GRN (n=169 TFs). 

# >>> Required packages
library(ComplexHeatmap)
library(BuenColors)
library(circlize)
library(magick)

# >>> Input data
# TFBS data
tfbs <- read.table(paste0(data_wd,"/osfstorage-archive/03_TFs_OCR_links/chromVAR_tf_cCREs/TFBS_scores_from_ocr_tf_matrix_curated.txt"))
atac_metadata <- readRDS(paste0(data_wd,"/osfstorage-archive/01_atac_seq_data/atac_metadata.rds"))
# GRN data
tf_ocr_gene <- read.table(paste0(data_wd,"/osfstorage-archive/04_gene_regulatory_networks/tf_ocr_gene_interactions_curated.txt"),
                          sep="\t",
                          dec=",",
                          header=TRUE)

# >>> Plotting
data_to_plot <- tfbs[unique(tf_ocr_gene$tf_symbol),]

# TFs to highlight
target_TF <- c("SPIB","IRF4","SPI1","BCL11A","JUND","FLI1","PAX5","FOXO4","FOXO1",
               "GFI1","IKZF1","TAL1","GATA1","MYB","MLXIP","RUNX3","ELK3","ETV6",
               "ETS2","ERG","E2F3","RUNX1","EBF1","TCF3","TCF4","LEF1")

target_TF_position <- vector()
for(i in 1:length(target_TF)){
  target_TF_position <- c(target_TF_position,grep(target_TF[i],rownames(data_to_plot)))
}

# Plot colors
bcell_colors <- readRDS(paste0(data_wd,"/osfstorage-archive/00_plot_parameters/bcell_colors.rds"))
tfbs_colors <- readRDS(paste0(data_wd,"/osfstorage-archive/00_plot_parameters/tfbs_colors.rds"))

# Plot
exprs_pattern <- Heatmap(data_to_plot, 
                         name = "TFBS patterns",
                         top_annotation = HeatmapAnnotation(stage = atac_metadata$CellType,
                                                            col = list(stage = bcell_colors)), 
                         column_order=order(atac_metadata$CellType),
                         right_annotation = rowAnnotation(foo = anno_mark(at = target_TF_position, 
                                                                          labels = rownames(data_to_plot)[target_TF_position],
                                                                          labels_gp = gpar(fontsize = 6))),
                         col=colorRamp2(seq(-1*min(abs(min(data_to_plot)),max(data_to_plot)), min(abs(min(data_to_plot)),max(data_to_plot)), length.out=9),tfbs_colors),
                         show_column_names = FALSE, 
                         show_row_names = FALSE, 
                         cluster_rows = TRUE,
                         width = unit(4, "cm"))

pdf(paste0(format(Sys.time(), "%Y%m%d"),"_heatmap_tfbs_grn.pdf"))
exprs_pattern
dev.off()
# ------------------------------------------------------------------------------------ #



# Figure S6 - Panel D -----------------------------------------------------------------
# Heatmap showing the TFs that are included in each cell-specific GRN (csGRN). 

# >>> Required packages
library(ComplexHeatmap)
library(igraph)

# >>> Input data
# Gene Expression data
rna_anno <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/rna_gene_annotation.rds"))
# GRN data
tf_ocr_gene <- read.table(paste0(data_wd,"/osfstorage-archive/04_gene_regulatory_networks/tf_ocr_gene_interactions_curated.txt"),
                          sep="\t",
                          dec=",",
                          header=TRUE)
grn <- readRDS(paste0(data_wd,"/osfstorage-archive/04_gene_regulatory_networks/bcell_grn_by_cell_type.rds"))

# >>> Computing
# Binary matrix showing if TF is included or not in a csGRN.
list_tfs <- unique(tf_ocr_gene$tf_ensembl)
binary_tfs <- matrix(0,nrow=length(list_tfs),ncol=length(grn),
                     dimnames=list(list_tfs,names(grn)))

for(i in 1:length(grn)){
  binary_tfs[rownames(binary_tfs)%in%V(grn[[i]])$name,i] <- 1
}

# >>> Plotting
rna_anno_tf <- rna_anno[rownames(binary_tfs),]
data_to_plot <- binary_tfs
rownames(data_to_plot) <- rna_anno_tf$hgnc_symbol

# Plot colors
bcell_colors <- readRDS(paste0(data_wd,"/osfstorage-archive/00_plot_parameters/bcell_colors.rds"))

# Plot
TFs_grns <- Heatmap(data_to_plot, 
                    top_annotation = HeatmapAnnotation(stage = colnames(data_to_plot),
                                                       col = list(stage = bcell_colors)),
                    col=c("lightgray","orange"),
                    rect_gp = gpar(col = "white"),
                    show_column_names = FALSE, 
                    show_row_names = TRUE, 
                    cluster_rows = TRUE,
                    cluster_columns = FALSE, 
                    use_raster=TRUE, 
                    row_title_rot = 0,
                    width = unit(4, "cm"),
                    row_names_gp = gpar(fontsize = 5))

pdf(paste0(format(Sys.time(), "%Y%m%d"),"_bcell_GRN_TFs_snapshot.pdf"),height=10, width=6)
TFs_grns
dev.off()
# ------------------------------------------------------------------------------------ #



# Figure S7 - Panel A -----------------------------------------------------------------
#  Regulon activity, based on RNA-Seq expression score based on GSVA and TF expression.
# (Similar to plot from SCENIC+)
# Heat map/dot-plot showing, for each TF regulon with an activator profile, 
# the TF expression on a color scale and the Regulon Specificity Score (RSS) of target regions on a size scale.

# >>> Required packages
library(GSVA)
library(ComplexHeatmap)
library(BuenColors)
library(ggplot2)
library(reshape2)

# >>> Required auxiliar functions
source(paste0(data_wd,"/auxiliar_functions.R"))

# >>> Input data
# Gene Expression data
rna <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/rna_norm_data.rds"))
rna_metadata <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/rna_metadata.rds"))
rna_anno <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/rna_gene_annotation.rds"))
# GRN data
tf_ocr_gene <- read.table(paste0(data_wd,"/osfstorage-archive/04_gene_regulatory_networks/tf_ocr_gene_interactions_curated.txt"),
                          sep="\t",
                          dec=",",
                          header=TRUE)

# >>> Computing
# Get the list of TFs in the B-cell GRN
tf_regulons <- unique(tf_ocr_gene[,c("tf_symbol","tf_ensembl")])

# Compute the accessibility of the atac for each regulon spliting by positive/negative correlation
gsva_tf_regulons_pos_mean <- matrix(NA, nrow=nrow(tf_regulons), ncol=length(levels(rna_metadata$CellType)),
                                    dimnames = list(tf_regulons$tf_ensembl, levels(rna_metadata$CellType)))
gsva_tf_regulons_neg_mean <- matrix(NA, nrow=nrow(tf_regulons), ncol=length(levels(rna_metadata$CellType)),
                                    dimnames = list(tf_regulons$tf_ensembl, levels(rna_metadata$CellType)))
gsva_tf_regulons_pos <- matrix(NA, nrow=nrow(tf_regulons), ncol=ncol(rna),
                               dimnames = list(tf_regulons$tf_ensembl, colnames(rna)))
gsva_tf_regulons_neg <- matrix(NA, nrow=nrow(tf_regulons), ncol=ncol(rna),
                               dimnames = list(tf_regulons$tf_ensembl, colnames(rna)))
for(i in tf_regulons$tf_ensembl){
  #Identify the TF-regulon
  target_tf_regulon <- tf_ocr_gene[tf_ocr_gene$tf_ensembl==i,]
  dim(target_tf_regulon)
  # Select positive regulation
  target_tf_regulon_pos <- target_tf_regulon[target_tf_regulon$gene_tf_rho>0,]
  dim(target_tf_regulon_pos)
  # Select negative regulation
  target_tf_regulon_neg <- target_tf_regulon[target_tf_regulon$gene_tf_rho<0,]
  dim(target_tf_regulon_neg)
  # Compute the GSVA on CREs for pos/neg regulons
  gs <- list(pos=unique(target_tf_regulon_pos$gene_ensembl),
             neg=unique(target_tf_regulon_neg$gene_ensembl))
  gsvaPar <- gsvaParam(as.matrix(rna), gs)
  gsva.es <- gsva(gsvaPar, verbose=FALSE)
  gsva_mean_by_celltype <- apply(gsva.es,1,FUN=function(x){
    mm <- by(x,rna_metadata$CellType,mean)
    return(as.numeric(mm))
  })
  rownames(gsva_mean_by_celltype) <- levels(rna_metadata$CellType)
  gsva_tf_regulons_pos_mean[i,] <- gsva_mean_by_celltype[,"pos"]
  gsva_tf_regulons_neg_mean[i,] <- gsva_mean_by_celltype[,"neg"]
  gsva_tf_regulons_pos[i,] <- gsva.es["pos",]
  gsva_tf_regulons_neg[i,] <- gsva.es["neg",]
}

# Defined thresholds: 0.2
## Binarize the rna_seq matrix

gsva_tf_regulons_pos_binary <- matrix(0, nrow=nrow(tf_regulons), ncol=ncol(rna),
                                      dimnames = list(tf_regulons$tf_ensembl, colnames(rna)))

gsva_tf_regulons_pos_binary[gsva_tf_regulons_pos>0.2] <- 1
rownames(gsva_tf_regulons_pos_binary) <- tf_regulons$tf_symbol

gsva_tf_regulons_neg_binary <- matrix(0, nrow=nrow(tf_regulons), ncol=ncol(rna),
                                      dimnames = list(tf_regulons$tf_ensembl, colnames(rna)))

gsva_tf_regulons_neg_binary[gsva_tf_regulons_neg>0.2] <- 1
rownames(gsva_tf_regulons_neg_binary) <- tf_regulons$tf_symbol

# Calculates Regulon specificity score (RSS)
rrs_pos <- calculate_rss(metadata=rna_metadata, binary_regulons=gsva_tf_regulons_pos_binary, cell_type_column="CellType")
rrs_neg <- calculate_rss(metadata=rna_metadata, binary_regulons=gsva_tf_regulons_neg_binary, cell_type_column="CellType")

# >>> Plotting
# Z-score norm RNA
data_to_plot <- t(scale(t(rna)))
rna_anno2 <- rna_anno[rownames(data_to_plot),]
table(rownames(data_to_plot)==rna_anno2$ensembl_gene_id)
rownames(data_to_plot) <- rna_anno2$hgnc_symbol
data_to_plot <- data_to_plot[rownames(data_to_plot)%in%tf_regulons$tf_symbol,]

# Mean z-score norm RNA
data_to_plot2 <- apply(data_to_plot,1,FUN=function(x){a <- by(x,rna_metadata$CellType,mean); return(as.numeric(a))})
rownames(data_to_plot2) <- levels(rna_metadata$CellType)
data_to_plot2 <- t(data_to_plot2)

data1 <- melt(data_to_plot2)

for(i in 1:nrow(data1)){
  data1$rrs_pos[i] <- rrs_pos$RSS[which(rrs_pos$regulon==data1$Var1[i] & rrs_pos$cell_type==data1$Var2[i])]
  data1$rrs_neg[i] <- rrs_neg$RSS[which(rrs_neg$regulon==data1$Var1[i] & rrs_neg$cell_type==data1$Var2[i])]
}

# Plot colors
bcell_colors <- readRDS(paste0(data_wd,"/osfstorage-archive/00_plot_parameters/bcell_colors.rds"))

# Get the data clustering based on positive RSS
regulons_group <- readRDS(paste0(data_wd,"/osfstorage-archive/04_gene_regulatory_networks/regulons/regulons_clusters.rds"))
regulons_group_plot <- c(rep("1",length(regulons_list[[1]])),rep("2",length(regulons_list[[2]])),
                         rep("3",length(regulons_list[[3]])),rep("4",length(regulons_list[[4]])),rep("5",length(regulons_list[[5]])))
names(regulons_group_plot) <- c(regulons_list[[1]],regulons_list[[2]],regulons_list[[3]],regulons_list[[4]],regulons_list[[5]])
regulons_color <- c("1"="#92c5de","2"="#2166ac","3"="#91cf60","4"="#dd1c77","5"="black")


# Include the cluster info in the formatted data to plot (data1)
data1$Var1 <- factor(data1$Var1, levels=rev(c(rev(regulons_list[[1]]),rev(regulons_list[[2]]),rev(regulons_list[[3]]),rev(regulons_list[[4]]),rev(regulons_list[[5]]))))

# Plot
# RSS tf-gene positively correlated (putative activator TFs)
ht_pos <- ggplot(data=data1, mapping=aes(x=Var2, y=Var1)) + 
  geom_raster(aes(fill=value)) +
  scale_fill_gradientn(colors=jdb_palette("brewer_celsius")) +
  geom_point(mapping=aes(size=rrs_pos), color="black") +
  scale_radius(range=c(min.size=0.1, max.size=2)) +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.title.y=element_blank(), 
        axis.text.x=element_text(angle=90, hjust=1, size=6),
        axis.text.y=element_text(size=6))

ggsave(plot = ht_pos, width = 4, height = 14, filename = paste0(format(Sys.time(), "%Y%m%d"),"_gene_RSS_tf_exprss_pos.pdf")) 
# ------------------------------------------------------------------------------------ #







# Figure S7 - Panel B -----------------------------------------------------------------
# Heatmap representation of gene expression for TF regulators of EBF1.

# >>> Required packages
library(ComplexHeatmap)
library(circlize)
library(BuenColors)

# >>> Input data
# Gene Expression data
rna <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/rna_norm_data.rds"))
rna_metadata <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/rna_metadata.rds"))
rna_anno <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/rna_gene_annotation.rds"))
# GRN data
tf_ocr_gene <- read.table(paste0(data_wd,"/osfstorage-archive/04_gene_regulatory_networks/tf_ocr_gene_interactions_curated.txt"),
                          sep="\t",
                          dec=",",
                          header=TRUE)

# >>> Computing
# Upstream regulators
target_gene <- "EBF1"
target_grn <- unique(tf_ocr_gene[tf_ocr_gene$gene_symbol==target_gene,c("gene_ensembl","tf_ensembl","gene_tf_rho")])

# Data to plot
M <- rna[target_grn$tf_ensembl,]
rownames(M) <- rna_anno[rownames(M),"hgnc_symbol"]

# >>> Plotting
# Plot colors
bcell_colors <- readRDS(paste0(data_wd,"/osfstorage-archive/00_plot_parameters/bcell_colors.rds"))
col_fun1 = colorRamp2(c(-0.9,-0.8,-0.7,-0.6, 0, 0.6, 0.7, 0.8, 0.9), jdb_palette("solar_flare",9,"discrete"))

# Plot
# Row annotation
target_genes <- c("PAX5", "EBF1", "TCF3", "FOXO1", "SPIB", "RUNX1", "CEBPA", "LEF1", "MEF2A", "IRF4", "TCF4", "POU4F1", "FOXO6", "MLXIP")
sel_pos <- vector("numeric",length=length(target_genes))
for(i in 1:length(target_genes)){
  sel_pos[i] <- which(rownames(M)==target_genes[i])
}

plot <- Heatmap(M, 
                top_annotation = HeatmapAnnotation(stage = rna_metadata$"CellType", 
                                                   col = list(stage = bcell_colors),
                                                   rna = anno_barplot(rna[rna_anno$ensembl_gene_id[rna_anno$hgnc_symbol==target_gene],])), 
                right_annotation=rowAnnotation(rho=target_grn$gene_tf_rho,
                                               col = list(rho = col_fun1),
                                               genes=anno_mark(at = sel_pos, labels = target_genes)), 
                show_row_names = TRUE, 
                show_column_names =FALSE, 
                column_names_gp = gpar(fontsize = 5), 
                row_names_gp = gpar(fontsize = 6),
                col=jdb_palette("wolfgang_extra"), 
                heatmap_legend_param = list(title = "Gene expression"),
                column_order = order(rna_metadata[,"CellType"]),
                width = unit(4, "cm"),
                height = unit(3, "mm")*nrow(M),
                row_names_side = "left",
                show_row_dend = FALSE)

pdf(paste0(format(Sys.time(), "%Y%m%d"),"_",target_gene,"_upstream_TFs_heatmap.pdf"), height = 13)
draw(plot)
dev.off()
# ------------------------------------------------------------------------------------ #



# Figure S7 - Panel C -----------------------------------------------------------------
# Heatmap representation of gene expression for TF regulated by EBF1.

# >>> Required packages
library(ComplexHeatmap)
library(circlize)
library(BuenColors)

# >>> Input data
# Gene Expression data
rna <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/rna_norm_data.rds"))
rna_metadata <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/rna_metadata.rds"))
rna_anno <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/rna_gene_annotation.rds"))
# TF data
ocr_tf_matrix <- read.table(paste0(data_wd,"/osfstorage-archive/03_TFs_OCR_links/ocr_tf_matrix.txt"),
                            sep="\t",
                            header=TRUE)
# GRN data
tf_ocr_gene <- read.table(paste0(data_wd,"/osfstorage-archive/04_gene_regulatory_networks/tf_ocr_gene_interactions_curated.txt"),
                          sep="\t",
                          dec=",",
                          header=TRUE)

# >>> Computing
#Genes regulated by EBF1 according to our network.
target_gene <- "EBF1"
list_target_genes_EBF1 <- unique(tf_ocr_gene[which(tf_ocr_gene$tf_symbol==target_gene),"gene_ensembl"])
list_target_tfs_EBF1 <- list_target_genes_EBF1[list_target_genes_EBF1%in%unique(rna_anno$ensembl_gene_id[rna_anno$hgnc_symbol%in%colnames(ocr_tf_matrix)])]
target_grn <- unique(tf_ocr_gene[which(tf_ocr_gene$tf_symbol==target_gene & tf_ocr_gene$gene_ensembl%in%list_target_tfs_EBF1),c("gene_ensembl","tf_ensembl","gene_tf_rho")])


# Data to plot
M <- rna[target_grn$gene_ensembl,]
rownames(M) <- rna_anno[rownames(M),"hgnc_symbol"]

# >>> Plotting
# Plot colors
bcell_colors <- readRDS(paste0(data_wd,"/osfstorage-archive/00_plot_parameters/bcell_colors.rds"))
col_fun1 = colorRamp2(c(-0.9,-0.8,-0.7,-0.6, 0, 0.6, 0.7, 0.8, 0.9), jdb_palette("solar_flare",9,"discrete"))

# Plot
# Row annotation
target_genes <- c("PAX5", "TCF3", "FOXO1", "MEF2A", "IRF4", "TCF4", "POU2F1", "EBF1","STAT5A")
sel_pos <- vector("numeric",length=length(target_genes))
for(i in 1:length(target_genes)){
  sel_pos[i] <- which(rownames(M)==target_genes[i])
}

plot <- Heatmap(M, 
                top_annotation = HeatmapAnnotation(stage = rna_metadata$"CellType", 
                                                   col = list(stage = bcell_colors),
                                                   rna = anno_barplot(rna[rna_anno$ensembl_gene_id[rna_anno$hgnc_symbol==target_gene],])), 
                right_annotation=rowAnnotation(rho=target_grn$gene_tf_rho,
                                               col = list(rho = col_fun1),
                                               genes=anno_mark(at = sel_pos, labels = target_genes)), 
                show_row_names = TRUE, 
                show_column_names =FALSE, 
                column_names_gp = gpar(fontsize = 5), 
                row_names_gp = gpar(fontsize = 6),
                col=jdb_palette("wolfgang_extra"), 
                heatmap_legend_param = list(title = "Gene expression"),
                column_order = order(rna_metadata[,"CellType"]),
                width = unit(4, "cm"),
                height = unit(3, "mm")*nrow(M),
                row_names_side = "left",
                show_row_dend = FALSE)

pdf(paste0(format(Sys.time(), "%Y%m%d"),"_",target_gene,"_downstream_TFs_heatmap.pdf"))
draw(plot)
dev.off()
# ------------------------------------------------------------------------------------ #



# Figure S7 - Panel D,E & Table S7 -----------------------------------------------------------------
# Representation of the enriched GO terms derived from the pathway analysis of the EBF1-regulon-activated genes and EBF1-regulon-repressed genes for each csGRN.

# >>> Required packages
library(igraph)
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)

# >>> Input data
# Gene Expression data
rna <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/rna_norm_data.rds"))
rna_anno <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/rna_gene_annotation.rds"))
# GRN data
tf_ocr_gene <- read.table(paste0(data_wd,"/osfstorage-archive/04_gene_regulatory_networks/tf_ocr_gene_interactions_curated.txt"),
                          sep="\t",
                          dec=",",
                          header=TRUE)
bcell_csGRN <- readRDS(paste0(data_wd,"/osfstorage-archive/04_gene_regulatory_networks/bcell_grn_by_cell_type.rds"))

# >>> Computing
#Genes regulated by EBF1 according to our network.
target_gene <- "EBF1"
list_target_genes_EBF1 <- unique(tf_ocr_gene[which(tf_ocr_gene$tf_symbol==target_gene),])

# Split by positive or negative regulation
list_target_genes_EBF1_pos <- unique(list_target_genes_EBF1[list_target_genes_EBF1$gene_tf_rho>0,"gene_ensembl"])
list_target_genes_EBF1_neg <- unique(list_target_genes_EBF1[list_target_genes_EBF1$gene_tf_rho<0,"gene_ensembl"])

# EBF1-regulon ORA analysis from csGRN
list_target_genes_activator_EBF1_by_cell_type <- vector(mode="list", length=length(bcell_csGRN))
names(list_target_genes_activator_EBF1_by_cell_type) <- names(bcell_csGRN)

list_target_genes_repressor_EBF1_by_cell_type <- vector(mode="list", length=length(bcell_csGRN))
names(list_target_genes_repressor_EBF1_by_cell_type) <- names(bcell_csGRN)

for(i in names(bcell_csGRN)){
  aa <- as_edgelist(bcell_csGRN[[i]])
  bb <- unique(aa[aa[,1]==rna_anno$ensembl_gene_id[rna_anno$hgnc_symbol==target_gene],2])
  if(length(bb)>0){
    list_target_genes_activator_EBF1_by_cell_type[[i]] <- bb[bb%in%list_target_genes_EBF1_pos]
    list_target_genes_repressor_EBF1_by_cell_type[[i]] <- bb[bb%in%list_target_genes_EBF1_neg]
  }
}

str(list_target_genes_activator_EBF1_by_cell_type)
str(list_target_genes_repressor_EBF1_by_cell_type)

# activator EBF1-regulon per csGRN (Table S7)
ck_go_005_activator <- compareCluster(geneCluster = list_target_genes_activator_EBF1_by_cell_type, 
                                      fun = enrichGO,
                                      universe      = rownames(rna),
                                      OrgDb         = org.Hs.eg.db,
                                      keyType       = 'ENSEMBL',
                                      ont           = "BP",
                                      pAdjustMethod = "BH",
                                      minGSSize = 20,
                                      maxGSSize = 500,
                                      pvalueCutoff  = 0.05,
                                      readable      = TRUE)

# repressor EBF1-regulon per csGRN ((Table S7))
ck_go_005_repressor <- compareCluster(geneCluster = list_target_genes_repressor_EBF1_by_cell_type, 
                                      fun = enrichGO,
                                      universe      = rownames(rna),
                                      OrgDb         = org.Hs.eg.db,
                                      keyType       = 'ENSEMBL',
                                      ont           = "BP",
                                      pAdjustMethod = "BH",
                                      minGSSize = 20,
                                      maxGSSize = 500,
                                      pvalueCutoff  = 0.05,
                                      readable      = TRUE)

# >>> Plotting
pdf(paste0(format(Sys.time(), "%Y%m%d"),"_ora_GO_activator_EBF1_regulon_from_csGRNs_dotplot_targeted.pdf"))
dotplot(ck_go_005_activator, size="Count", showCategory = c("chromosome segregation",
                                                            "B cell activation",
                                                            "B cell differentiation",
                                                            "B cell receptor signaling pathway",
                                                            "mature B cell differentiation",
                                                            "B cell proliferation",
                                                            "antigen receptor-mediated signaling pathway",
                                                            "lymphocyte differentiation",
                                                            "regulation of antigen receptor-mediated signaling pathway")) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf(paste0(format(Sys.time(), "%Y%m%d"),"_ora_GO_repressor_EBF1_regulon_from_csGRNs_dotplot.pdf"))
dotplot(ck_go_005_repressor, size="Count", showCategory =10) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
# ------------------------------------------------------------------------------------ #


# Figure S9 - Panel A-- -----------------------------------------------------------------
# Correlation between TF regulon score (Fig2A) from bulk and single-cell (Fig4C,D).

# >>> Required packages
library(GSVA)
library(ComplexHeatmap)
library(BuenColors)
library(ggplot2)
library(reshape2)
library(circlize)
library(ggrepel)

# >>> Input data
# Gene Expression data
rna <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/rna_norm_data.rds"))
rna_metadata <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/rna_metadata.rds"))

# TF regulon targets
df_regulons_list_pos <- readRDS(paste0(data_wd,"/osfstorage-archive/04_gene_regulatory_networks/regulons/TF_regulon_activated_gene_set_ensembl.rds"))

## AUCell single-cell B-cell differentiation
scpublic_bcell <- readRDS("C:/Users/nplanellpic/Documents/npp/project/ongoing/bcell/submissions/20250127_Scadvance_submission/20250723_revision_2ndRound/regulon_score_comparison/20250724_scPublic_metadata.rds")
dim(scpublic_bcell)

scinhouse_bcell <- readRDS("C:/Users/nplanellpic/Documents/npp/project/ongoing/bcell/submissions/20250127_Scadvance_submission/20250723_revision_2ndRound/regulon_score_comparison/20250724_scInHouse_metadata.rds")
dim(scinhouse_bcell)

# Regulons groups
regulons_group <- readRDS(paste0(data_wd,"/osfstorage-archive/04_gene_regulatory_networks/regulons/TF_regulons_clusters.rds"))
regulons_group_plot <- c(rep("1",length(regulons_group[[1]])),rep("2",length(regulons_group[[2]])),rep("3",length(regulons_group[[3]])),rep("4",length(regulons_group[[4]])),rep("5",length(regulons_group[[5]])))
names(regulons_group_plot) <- c(regulons_group[[1]],regulons_group[[2]],regulons_group[[3]],regulons_group[[4]],regulons_group[[5]])

# >>> Computing
gs_ensembl <- df_regulons_list_pos
gsvaPar <- gsvaParam(rna, gs_ensembl, kcdf="Gaussian")
gsva.Bcell <- gsva(gsvaPar, verbose=FALSE)


## Get the mean activity per cell type:
mean_TFactivity_bulk <- t(apply(gsva.Bcell,1,FUN=function(x){return(as.numeric(by(x,rna_metadata$CellType,mean)))}))
colnames(mean_TFactivity_bulk) <- levels(rna_metadata$CellType)
head(mean_TFactivity_bulk)

mean_TFactivity_scpublic_bcell <- t(apply(scpublic_bcell[,grep("^AUCell_",names(scpublic_bcell))],2,FUN=function(x){return(as.numeric(by(x,scpublic_bcell$cell_annotation,mean)))}))
colnames(mean_TFactivity_scpublic_bcell) <- levels(scpublic_bcell$cell_annotation)
rownames(mean_TFactivity_scpublic_bcell) <- gsub("^AUCell_","",rownames(mean_TFactivity_scpublic_bcell))
head(mean_TFactivity_scpublic_bcell)

mean_TFactivity_scinhouse_bcell <- t(apply(scinhouse_bcell[,grep("^AUCell_",names(scinhouse_bcell))],2,FUN=function(x){return(as.numeric(by(x,scinhouse_bcell$cell_type_reannotation,mean)))}))
colnames(mean_TFactivity_scinhouse_bcell) <- levels(scinhouse_bcell$cell_type_reannotation)
rownames(mean_TFactivity_scinhouse_bcell) <- gsub("^AUCell_","",rownames(mean_TFactivity_scinhouse_bcell))
head(mean_TFactivity_scinhouse_bcell)


## Compute the Spearman correlation and save the p-value and the rho

# Bulk vs scPublic ----------------------------
table(rownames(mean_TFactivity_bulk)==rownames(mean_TFactivity_scpublic_bcell))

TFregulon_bulk_scPublic_correlation <- matrix(nrow=nrow(mean_TFactivity_bulk),ncol=3, 
                                              dimnames = list(rownames(mean_TFactivity_bulk), c("rho","Pvalue","log10P")))

for(i in rownames(mean_TFactivity_bulk)){
  res <- cor.test(mean_TFactivity_bulk[i,c("HSC","CLP","proB","preB","ImmatureB","Transitional_B")],mean_TFactivity_scpublic_bcell[i,c("HSC/MPP","CLP","Pro-B","Pre-B","Immature-B","Mature-B")], method="spearman")
  TFregulon_bulk_scPublic_correlation[i,"rho"] <- res$estimate
  TFregulon_bulk_scPublic_correlation[i,"Pvalue"] <-   res$p.value
  TFregulon_bulk_scPublic_correlation[i,"log10P"] <- -log10(res$p.value)
}

## Plotting
# Plot colors
bcell_colors <- readRDS(paste0(data_wd,"/osfstorage-archive/00_plot_parameters/bcell_colors.rds"))

# Summary correlation plot (Panel A; left)
data_to_plot <- as.data.frame(TFregulon_bulk_scPublic_correlation)

regulons_color <- c("1"="#92c5de","2"="#2166ac","3"="#91cf60","4"="#dd1c77","5"="black")
data_to_plot$regulon_group <- regulons_group_plot[rownames(TFregulon_bulk_scPublic_correlation)]
data_to_plot$TF <- rownames(data_to_plot)

p1 <- ggplot(data_to_plot,aes(rho,log10P, color=regulon_group)) +
  geom_point(position = "jitter") +
  scale_color_manual(values = regulons_color) +
  geom_hline(yintercept = 1.3, linetype = "dashed", color = "black") +
  theme_bw() +
  geom_text_repel(
    data = subset(data_to_plot, Pvalue > 0.05), 
    aes(label = TF),
    box.padding = 0.5, point.padding = 0.5, segment.color = 'grey50',
    size = 3) +
  geom_text_repel(
    data = subset(data_to_plot, TF %in% c("IRX2","ELK3","MECOM","ZSCAN31","MEIS1","EBF1","PAX5","CEBPB")), 
    aes(label = TF),
    box.padding = 0.5, point.padding = 0.5, segment.color = 'grey50',
    size = 3) +
  ggtitle("Spearman correlation of TF-regulon score:\n B-cell (bulk) vs Public - phs002371.v7.2 (single-cell)")

ggsave(filename=paste0(format(Sys.time(), "%Y%m%d"),"_corrplot_bulk_scPublic.pdf"),plot=p1, width = 6, height = 5)

# Representative correlation plots (Panel B)
target_TFs <- c("IRX2","ELK3","MECOM","ZSCAN31","MEIS1","EBF1","PAX5","CEBPB","VDR","ARID5A","OTX1","AIRE","MLXIP","NR2F6","MXI1")
data_to_plot <- melt(mean_TFactivity_bulk[target_TFs,c("HSC","CLP","proB","preB","ImmatureB","Transitional_B")])
names(data_to_plot)[3] <- "bulk"
aa <- melt(mean_TFactivity_scpublic_bcell[target_TFs,c("HSC/MPP","CLP","Pro-B","Pre-B","Immature-B","Mature-B")])
table(data_to_plot$Var1==aa$Var1)
data_to_plot$sc <- aa$value

pp1 <- ggplot(data_to_plot, aes(bulk,sc)) +
  geom_line() +
  geom_point(aes(color=Var2), size=3) +
  facet_wrap(~Var1, scales = "free") +
  theme_bw() +
  ylab("Public dataset (AUCell)") +
  xlab("B cell (bulk - GSVA)") +
  labs(color = "Cell type") +
  scale_color_manual(values=bcell_colors) +
  ggtitle("TF-regulon score:\n B-cell (bulk) vs Public - phs002371.v7.2 (single-cell)")


ggsave(filename=paste0(format(Sys.time(), "%Y%m%d"),"_scatterplot_bulk_scPublic.pdf"),plot=pp1, width=8)


# Bulk vs scInHouse ----------------------------
table(rownames(mean_TFactivity_bulk)==rownames(mean_TFactivity_scinhouse_bcell))

TFregulon_bulk_scinhouse_bcell_correlation <- matrix(nrow=nrow(mean_TFactivity_bulk),ncol=3, 
                                                     dimnames = list(rownames(mean_TFactivity_bulk), c("rho","Pvalue","log10P")))

for(i in rownames(mean_TFactivity_bulk)){
  res <- cor.test(mean_TFactivity_bulk[i,c("HSC","CLP","proB","preB")],mean_TFactivity_scinhouse_bcell[i,c("HSC/MPP","CLP/Pre-pro-B","pro-B","pre-B")], method="spearman")
  TFregulon_bulk_scinhouse_bcell_correlation[i,"rho"] <- res$estimate
  TFregulon_bulk_scinhouse_bcell_correlation[i,"Pvalue"] <-   res$p.value
  TFregulon_bulk_scinhouse_bcell_correlation[i,"log10P"] <- -log10(res$p.value)
}

## Plotting
# Summary correlation plot (Panel A; right)
data_to_plot <- as.data.frame(TFregulon_bulk_scinhouse_bcell_correlation)
data_to_plot$regulon_group <- regulons_group_plot[rownames(TFregulon_bulk_scPublic_correlation)]
data_to_plot$TF <- rownames(data_to_plot)

p2 <- ggplot(data_to_plot,aes(rho,log10P, color=regulon_group)) +
  geom_point(position = "jitter") +
  scale_color_manual(values = regulons_color) +
  theme_bw() +
  geom_text_repel(
    data = subset(data_to_plot, Pvalue > 0.05), 
    aes(label = TF),
    box.padding = 0.5, point.padding = 0.5, segment.color = 'grey50',
    size = 3) +
  geom_text_repel(
    data = subset(data_to_plot, TF %in% c("IRX2","ELK3","MECOM","ZSCAN31","MEIS1","EBF1","PAX5","CEBPB")), 
    aes(label = TF),
    box.padding = 0.5, point.padding = 0.5, segment.color = 'grey50',
    size = 3) +
  ggtitle("Spearman correlation of TF-regulon score:\n B-cell (bulk) vs Inhouse - CD34+ BM (single-cell)")


ggsave(filename=paste0(format(Sys.time(), "%Y%m%d"),"_corrplot_bulk_scInhouse.pdf"),plot=p2, width = 6, height = 5)


# Representative correlation plots (Panel C)
target_TFs <- c("IRX2","ELK3","MECOM","ZSCAN31","MEIS1","EBF1","PAX5","CEBPB","VDR","ARID5A","OTX1","AIRE","MLXIP","NR2F6","MXI1")
data_to_plot <- melt(mean_TFactivity_bulk[target_TFs,c("HSC","CLP","proB","preB")])
names(data_to_plot)[3] <- "bulk"
aa <- melt(mean_TFactivity_scinhouse_bcell[target_TFs,c("HSC/MPP","CLP/Pre-pro-B","pro-B","pre-B")])
table(data_to_plot$Var1==aa$Var1)
data_to_plot$sc <- aa$value

pp2 <- ggplot(data_to_plot, aes(bulk,sc)) +
  geom_line() +
  geom_point(aes(color=Var2), size=3) +
  facet_wrap(~Var1, scales = "free") +
  theme_bw() +
  ylab("Public dataset (AUCell)") +
  xlab("B cell (bulk - GSVA)") +
  labs(color = "Cell type") +
  scale_color_manual(values=bcell_colors) +
  ggtitle("TF-regulon score:\n B-cell (bulk) vs Inhouse - CD34+ BM (single-cell)")


ggsave(filename=paste0(format(Sys.time(), "%Y%m%d"),"_scatterplot_bulk_scInhouse.pdf"),plot=pp2, width=8)

# ------------------------------------------------------------------------------------ #



# Figure S9 - Panel B-C -----------------------------------------------------------------

# ------------------------------------------------------------------------------------ #


# Figure S10 - Panel A-C -----------------------------------------------------------------
# >>> Exploratory analysis Li dataset
# Principal component analysis
pca_output <- prcomp(t(ball_norm_rna_Li))
pca_output_summary <- summary(pca_output)

data_to_plot <- data.frame(PC1=pca_output$x[,1],PC2=pca_output$x[,2], ALL_subgroup=ball_metadata_Li$mutation, batch=ball_metadata_Li$batch, age=ball_metadata_Li$age)

# Plotting
# Plot colors
ball_colors <- jdb_palette("corona", length(levels(ball_metadata_Li$mutation)))
names(ball_colors) <- levels(ball_metadata_Li$mutation)

age_color <- c("Paediatric"="#DAA520","Adult"="#2F4F4F")


# Plot
p_subgroup <- ggplot(data_to_plot, aes(x = PC1, y = PC2, color = ALL_subgroup)) +
  geom_point(size=3) +
  scale_color_manual(values=ball_colors) +
  labs(title="B-ALL Li J. multi-dataset", x = paste0("PC1 (",round(pca_output_summary$importance[2,1]*100,1),")%"), y = paste0("PC2 (",round(pca_output_summary$importance[2,2]*100,1),")%")) +
  theme_classic()

p_batch <- ggplot(data_to_plot, aes(x = PC1, y = PC2, color = batch)) +
  geom_point(size=3) +
  labs(title="B-ALL Li J. multi-dataset", x = paste0("PC1 (",round(pca_output_summary$importance[2,1]*100,1),")%"), y = paste0("PC2 (",round(pca_output_summary$importance[2,2]*100,1),")%")) +
  theme_classic()

p_age <- ggplot(data_to_plot, aes(x = PC1, y = PC2, color = age)) +
  geom_point(size=3) +
  scale_color_manual(values=age_color) +
  labs(title="B-ALL Li J. multi-dataset", x = paste0("PC1 (",round(pca_output_summary$importance[2,1]*100,1),")%"), y = paste0("PC2 (",round(pca_output_summary$importance[2,2]*100,1),")%")) +
  theme_classic()

pdf(paste(format(Sys.time(), "%Y%m%d"),"_pca_BALL_Li.pdf",sep=""), width=6, height = 6)
print(p_subgroup)
print(p_subgroup+theme(legend.position = "none"))
print(p_batch)
print(p_batch+theme(legend.position = "none"))
print(p_age)
print(p_age+theme(legend.position = "none"))
dev.off()
# ------------------------------------------------------------------------------------ #




# Figure S10 - Panel D -----------------------------------------------------------------
# B-ALL subtypes signatures
# Heatmap representation of 5,749 genes defining the signature profile of the B-ALL subtypes. 
# Biomarker genes for each B-ALL subtype were derived from 
# significant genes (logFC≥1.5 and adjusted p-value≤0.05) in at least 6 contrasts (see methods).

# >>> Required packages
library(ComplexHeatmap)
library(BuenColors)
library(tidyverse)

# >>> Input data
# B-ALL Gene Expression data: Li dataset
load("BALL_Li.RData", verbose=TRUE)

# Human B-ALL subtype biomarkers derived from Li dataset
ball_biomarkers <- readRDS(paste0(data_wd,"/osfstorage-archive/05_B-ALL/ALL_subgroup_biomarkers.rds"))
#ball_biomarkers <- readRDS("ALL_subgroup_biomarkers.rds")
str(ball_biomarkers)

data_to_plot <- ball_norm_rna_Li[rownames(ball_norm_rna_Li)%in%unique(unlist(ball_biomarkers)),]

# >>> Plotting
# Plot colors
ball_colors <- jdb_palette("corona", length(levels(ball_metadata_Li$mutation)))
names(ball_colors) <- levels(ball_metadata_Li$mutation)

Zscore <- t(scale(t(ball_norm_rna_Li[biomarkers_BALL,]))) 
dim(Zscore)

cols <- jdb_palette("solar_extra")
col_fun <- colorRamp2(seq(-4, 4, length.out = length(cols)), cols)

h_BALL_biomarkers_scaled <- Heatmap(Zscore,
                                    top_annotation = HeatmapAnnotation("B-ALL subtype" = ball_metadata_Li$mutation,
                                                                       col=list("B-ALL subtype" = ball_colors)),
                                    col = col_fun,
                                    column_split = group,
                                    column_gap = unit(0.5, "mm"),
                                    show_column_names = FALSE,
                                    show_row_names = FALSE,
                                    show_row_dend = FALSE,
                                    show_column_dend = FALSE,
                                    row_title = NULL,
                                    name="expression")

pdf(paste(format(Sys.time(), "%Y%m%d"),"_BALL_subgoups_markers_heatmap.pdf",sep=""))
h_BALL_biomarkers_scaled
dev.off()
# ------------------------------------------------------------------------------------ #




# Figure S10 - Panel E -----------------------------------------------------------------
# B-cell signature upset plot
# Upset plot showing the overlap between the gene expression profile of each B-cell subpopulation.

# >>> Required packages
library(ComplexHeatmap)
library(dplyr)

# >>> Input data
# Human B-ALL subtype biomarkers
ball_biomarkers <- readRDS(paste0(data_wd,"/osfstorage-archive/05_B-ALL/ALL_subgroup_biomarkers.rds"))
#ball_biomarkers <- readRDS("ALL_subgroup_biomarkers.rds")
str(ball_biomarkers)

# >>> Plotting
#UpSet plot
m1 = make_comb_mat(ball_biomarkers)
BALL_subgoups_markers_upsetplot <- UpSet(m1,comb_order = order(comb_size(m1), decreasing = TRUE), 
                                         set_order = c("MEF2D fusions",
                                                       "TCF3-PBX1",
                                                       "ETV6-RUNX1/like",
                                                       "DUX4 fusions",
                                                       "ZNF384 fusions",
                                                       "BCR-ABL1/ph-like",
                                                       "Hyperdiploidy",
                                                       "KTM2A fusions"),
                                         pt_size = unit(2,"mm"))

pdf(paste0(format(Sys.time(), "%Y%m%d"),"_upset_ball_biomarkers.pdf"), height = 5)
print(BALL_subgoups_markers_upsetplot)
dev.off()
# ------------------------------------------------------------------------------------ #


# Figure S10 - Panel F-H -----------------------------------------------------------------
# >>> Exploratory analysis Rainer dataset
# Principal component analysis
pca_output <- prcomp(t(ball_norm_rna_Rainer))
pca_output_summary <- summary(pca_output)

data_to_plot <- data.frame(PC1=pca_output$x[,1],PC2=pca_output$x[,2], 
                           ALL_subgroup=ball_metadata_Rainer$mutation, 
                           timepoint=ball_metadata_Rainer$timepoint, 
                           patient=ball_metadata_Rainer$"source_name_ch1")

# Plotting
# Plot colors
ball_colors <- jdb_palette("corona", length(levels(ball_metadata_Li$mutation)))
names(ball_colors) <- levels(ball_metadata_Li$mutation)

timepoint_color <- c("0"="#7B68EE","6/8"="#FF8C00","24"="#C71585")


# Plot
p_subgroup <- ggplot(data_to_plot, aes(x = PC1, y = PC2, color = ALL_subgroup)) +
  geom_point(size=3) +
  scale_color_manual(values=ball_colors) +
  labs(title="B-ALL Rainer dataset", x = paste0("PC1 (",round(pca_output_summary$importance[2,1]*100,1),")%"), y = paste0("PC2 (",round(pca_output_summary$importance[2,2]*100,1),")%")) +
  theme_classic()

p_timepoint <- ggplot(data_to_plot, aes(x = PC1, y = PC2, color = timepoint)) +
  geom_point(size=3) +
  scale_color_manual(values=timepoint_color) +
  labs(title="B-ALL Rainer dataset", x = paste0("PC1 (",round(pca_output_summary$importance[2,1]*100,1),")%"), y = paste0("PC2 (",round(pca_output_summary$importance[2,2]*100,1),")%")) +
  theme_classic()

p_patient <- ggplot(data_to_plot, aes(x = PC1, y = PC2, color = patient)) +
  geom_point(size=3) +
  labs(title="B-ALL Rainer dataset", x = paste0("PC1 (",round(pca_output_summary$importance[2,1]*100,1),")%"), y = paste0("PC2 (",round(pca_output_summary$importance[2,2]*100,1),")%")) +
  theme_classic()

pdf(paste(format(Sys.time(), "%Y%m%d"),"_pca_BALL_Rainer.pdf",sep=""), width=6, height = 6)
print(p_subgroup)
print(p_subgroup+theme(legend.position = "none"))
print(p_timepoint)
print(p_timepoint+theme(legend.position = "none"))
print(p_patient)
print(p_patient+theme(legend.position = "none"))
dev.off()
# ------------------------------------------------------------------------------------ #


# Figure S10 - Panel I -----------------------------------------------------------------
# >>> Exploratory analysis Rainer dataset
# Principal component analysis
pca_output <- prcomp(t(ball_norm_rna_Iacobucci))
pca_output_summary <- summary(pca_output)

data_to_plot <- data.frame(PC1=pca_output$x[,1],PC2=pca_output$x[,2], 
                           ALL_subgroup=ball_metadata_Iacobucci$mutation)

# Plotting
# Plot colors
ball_colors <- c("BCR-ABL1/ph-like"="#8c564b","Others"="#ADD8E6")

# Plot
p_subgroup <- ggplot(data_to_plot, aes(x = PC1, y = PC2, color = ALL_subgroup)) +
  geom_point(size=3) +
  scale_color_manual(values=ball_colors) +
  labs(title="B-ALL Iacobucci dataset", x = paste0("PC1 (",round(pca_output_summary$importance[2,1]*100,1),")%"), y = paste0("PC2 (",round(pca_output_summary$importance[2,2]*100,1),")%")) +
  theme_classic()

pdf(paste(format(Sys.time(), "%Y%m%d"),"_pca_BALL_Iacobucci.pdf",sep=""), width=6, height = 6)
print(p_subgroup)
print(p_subgroup+theme(legend.position = "none"))
dev.off()
# ------------------------------------------------------------------------------------ #



# Figure S11 ---------------------------------------------------------------------------
# TF regulon activity of pro-B cell regulons across B-ALL subtypes

# >>> Required packages
library(ggplot2)
library(BuenColors)
library(tidyverse)
library(ggpubr)

# >>> Input data
# B-ALL Gene Expression datasets
load(paste0(data_wd,"/osfstorage-archive/05_B-ALL/BALL_Li.RData"), verbose = TRUE)

# Regulons group
regulons_group <- readRDS(paste0(data_wd,"/osfstorage-archive/04_gene_regulatory_networks/regulons/TF_regulons_clusters.rds"))
regulons_group_plot <- c(rep("1",length(regulons_group[[1]])),rep("2",length(regulons_group[[2]])),rep("3",length(regulons_group[[3]])),rep("4",length(regulons_group[[4]])),rep("5",length(regulons_group[[5]])))
names(regulons_group_plot) <- c(regulons_group[[1]],regulons_group[[2]],regulons_group[[3]],regulons_group[[4]],regulons_group[[5]])


# >>> Ploting
# Plot colors
ball_colors <- jdb_palette("corona", length(levels(ball_metadata_Li$mutation)))
names(ball_colors) <- levels(ball_metadata_Li$mutation)

# Plot  
# Boxplots
study <- c("Li")
gsva_data_list <- list("Li"=gsva.BALL_Li)
metadata_list <- list("Li"=ball_metadata_Li)
target_gsva <- regulons_group[[3]]


pdf(paste0(format(Sys.time(), "%Y%m%d"),"_boxplots_GSVA_regulons_BALL_paper_stats.pdf"), width=15, height=15)
i <- "Li"
dd <- gsva_data_list[[i]][rownames(gsva_data_list[[i]])%in%target_gsva,]
data_to_plot <- data.frame(melt(dd),ALL_subgroup=rep(metadata_list[[i]]$mutation, each=nrow(dd)))
ball_colors <- jdb_palette("corona", length(levels(ball_metadata_Li$mutation)))
names(ball_colors) <- levels(ball_metadata_Li$mutation)

p <- ggplot(data_to_plot,aes(x=ALL_subgroup,y=value, fill=ALL_subgroup)) +
  geom_boxplot(outlier.size = 1) +
  facet_wrap(~Var1, ncol=4) +
  scale_fill_manual(values=ball_colors) +
  theme_bw() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  stat_compare_means(method = "wilcox.test", comparisons = list(c("ETV6-RUNX1/like","MEF2D fusions"),
                                                                c("ETV6-RUNX1/like","TCF3-PBX1"),
                                                                c("ETV6-RUNX1/like","DUX4 fusions"),
                                                                c("ETV6-RUNX1/like","ZNF384 fusions"),
                                                                c("ETV6-RUNX1/like","BCR-ABL1/ph-like"),
                                                                c("ETV6-RUNX1/like","Hyperdiploidy"),
                                                                c("ETV6-RUNX1/like","KTM2A fusions")),
                     label = "p.signif") +
  ylab("GSVA score")
print(p)
dev.off()

# ------------------------------------------------------------------------------------ #



# Figure S12 - Panel A ----------------------------------------------------------------
# Heatmap representation of Pearson correlation matrices (rho) between samples based on RNA-seq data, 
# ATAC-seq data, ATAC-seq data at all promoter OCRs, ATAC-seq data at all intronic OCRs, and ATAC-seq data at all distal intergenic OCRs. 

# >>> Required packages
library(ComplexHeatmap)
library(circlize)

# >>> Input data
# Gene Expression data
rna <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/rna_norm_data.rds"))
rna_metadata <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/rna_metadata.rds"))

# Chromatin accessibility data
atac <- readRDS(paste0(data_wd,"/osfstorage-archive/01_atac_seq_data/atac_norm_data.rds"))
atac_metadata <- readRDS(paste0(data_wd,"/osfstorage-archive/01_atac_seq_data/atac_metadata.rds"))
atac_anno <- readRDS(paste0(data_wd,"/osfstorage-archive/01_atac_seq_data/atac_OCR_annotation.rds"))

# >>> Computing
rna_cor <- cor(rna, method = "pearson")
atac_cor <- cor(atac, method = "pearson")
atac_promoter_cor <- cor(atac[rownames(atac)%in%rownames(atac_anno)[atac_anno$annotation_simplified=="Promoter"],], method = "pearson")
atac_intron_cor <- cor(atac[rownames(atac)%in%rownames(atac_anno)[atac_anno$annotation_simplified=="Intron"],], method = "pearson")
atac_distal_cor <- cor(atac[rownames(atac)%in%rownames(atac_anno)[atac_anno$annotation_simplified=="Distal Intergenic"],], method = "pearson")

# >>> Plotting
# Plot colors
bcell_colors <- readRDS(paste0(data_wd,"/osfstorage-archive/00_plot_parameters/bcell_colors.rds"))
correlation_colors <- colorRamp2(seq(0.2, 1, length = 3), c("blue", "orange", "red"))

# Plot
rna_cor_plot <- Heatmap(rna_cor,
                     top_annotation = HeatmapAnnotation(Subpopulation = rna_metadata$CellType, 
                                                        col = list(Subpopulation = bcell_colors),
                                                        show_annotation_name = FALSE),
                     left_annotation = rowAnnotation(Subpopulation = rna_metadata$CellType, 
                                                        col = list(Subpopulation = bcell_colors),
                                                        show_annotation_name = FALSE),
                     col = correlation_colors,
                     cluster_columns = FALSE,
                     cluster_rows = FALSE,
                     column_names_gp = gpar(fontsize = 6),
                     column_order = order(rna_metadata$CellType),
                     row_order = order(rna_metadata$CellType),
                     show_column_names = FALSE,
                     show_row_names = FALSE,
                     column_title = "RNA-seq",
                     name="RNA-seq",
                     width = unit(8, "cm"), 
                     height = unit(8, "cm")
)

atac_cor_plot <- Heatmap(atac_cor,
                        top_annotation = HeatmapAnnotation(Subpopulation = atac_metadata$CellType, 
                                                           col = list(Subpopulation = bcell_colors),
                                                           show_annotation_name = FALSE),
                        left_annotation = rowAnnotation(Subpopulation = atac_metadata$CellType, 
                                                        col = list(Subpopulation = bcell_colors),
                                                        show_annotation_name = FALSE),
                        col = correlation_colors,
                        cluster_columns = FALSE,
                        cluster_rows = FALSE,
                        column_names_gp = gpar(fontsize = 6),
                        column_order = order(atac_metadata$CellType),
                        row_order = order(atac_metadata$CellType),
                        show_column_names = FALSE,
                        show_row_names = FALSE,
                        column_title = "ATAC-seq",
                        name="ATAC-seq",
                        width = unit(8, "cm"), 
                        height = unit(8, "cm")
)

atac_cor_plot_promoter <- Heatmap(atac_promoter_cor,
                         top_annotation = HeatmapAnnotation(Subpopulation = atac_metadata$CellType, 
                                                            col = list(Subpopulation = bcell_colors),
                                                            show_annotation_name = FALSE),
                         left_annotation = rowAnnotation(Subpopulation = atac_metadata$CellType, 
                                                         col = list(Subpopulation = bcell_colors),
                                                         show_annotation_name = FALSE),
                         col = correlation_colors,
                         cluster_columns = FALSE,
                         cluster_rows = FALSE,
                         column_names_gp = gpar(fontsize = 6),
                         column_order = order(atac_metadata$CellType),
                         row_order = order(atac_metadata$CellType),
                         show_column_names = FALSE,
                         show_row_names = FALSE,
                         column_title = paste0("Promoter OCRs (",round(sum(atac_anno$annotation_simplified=="Promoter")/nrow(atac_anno)*100,0),"%)"),
                         width = unit(8, "cm"), 
                         height = unit(8, "cm"))

atac_cor_plot_intron <- Heatmap(atac_intron_cor,
                         top_annotation = HeatmapAnnotation(Subpopulation = atac_metadata$CellType, 
                                                            col = list(Subpopulation = bcell_colors),
                                                            show_annotation_name = FALSE),
                         left_annotation = rowAnnotation(Subpopulation = atac_metadata$CellType, 
                                                         col = list(Subpopulation = bcell_colors),
                                                         show_annotation_name = FALSE),
                         col = correlation_colors,
                         cluster_columns = FALSE,
                         cluster_rows = FALSE,
                         column_names_gp = gpar(fontsize = 6),
                         column_order = order(atac_metadata$CellType),
                         row_order = order(atac_metadata$CellType),
                         show_column_names = FALSE,
                         show_row_names = FALSE,
                         column_title = paste0("Intronic OCRs (",round(sum(atac_anno$annotation_simplified=="Intron")/nrow(atac_anno)*100,0),"%)"),
                         width = unit(8, "cm"), 
                         height = unit(8, "cm"))

atac_cor_plot_distal <- Heatmap(atac_distal_cor,
                         top_annotation = HeatmapAnnotation(Subpopulation = atac_metadata$CellType, 
                                                            col = list(Subpopulation = bcell_colors),
                                                            show_annotation_name = FALSE),
                         left_annotation = rowAnnotation(Subpopulation = atac_metadata$CellType, 
                                                         col = list(Subpopulation = bcell_colors),
                                                         show_annotation_name = FALSE),
                         col = correlation_colors,
                         cluster_columns = FALSE,
                         cluster_rows = FALSE,
                         column_names_gp = gpar(fontsize = 6),
                         column_order = order(atac_metadata$CellType),
                         row_order = order(atac_metadata$CellType),
                         show_column_names = FALSE,
                         show_row_names = FALSE,
                         column_title=paste0("Distal Intergenic OCRs (",round(sum(atac_anno$annotation_simplified=="Distal Intergenic")/nrow(atac_anno)*100,0),"%)"),
                         width = unit(8, "cm"), 
                         height = unit(8, "cm"))

pdf(paste0(format(Sys.time(), "%Y%m%d"),"_rna_atac_correlation_heatmaps.pdf"))
  rna_cor_plot
  atac_cor_plot
  atac_cor_plot_promoter
  atac_cor_plot_intron
  atac_cor_plot_distal
dev.off()
# ------------------------------------------------------------------------------------ #



# Figure S12 - Panel B ----------------------------------------------------------------
# A t-SNE representation of OCRs identified in the study. Gini index is shown.

# >>> Required packages
library(umap)
library(ggplot2)
library(BuenColors)
library(REAT)

# >>> Input data
# Chromatin accessibility data
atac <- readRDS(paste0(data_wd,"/osfstorage-archive/01_atac_seq_data/atac_norm_data.rds"))
atac_metadata <- readRDS(paste0(data_wd,"/osfstorage-archive/01_atac_seq_data/atac_metadata.rds"))
atac_anno <- readRDS(paste0(data_wd,"/osfstorage-archive/01_atac_seq_data/atac_OCR_annotation.rds"))
atac_umap <- readRDS(paste0(data_wd,"/osfstorage-archive/01_atac_seq_data/atac_umap.rds"))

# >>> Computing
# UMAP
#set.seed(123)
#atac_umap <- umap(log2(atac+1))

# Gini Index
gini_output <- vector()
for(i in 1:nrow(atac)){
  dd <- data.frame(x=as.numeric(atac[i,]), y=atac_metadata$CellType)
  a <- by(dd$x,dd$y,max)
  r <- gini(as.numeric(a)) 
  gini_output <- c(gini_output,r)
}

# >>> Plotting
# Gini Index UMAP
m <- data.frame(umap1=atac_umap$layout[,1],umap2=atac_umap$layout[,2],gini_index=gini_output, peak_type=atac_anno[rownames(atac_umap$layout),"annotation_simplified"])

gg <- ggplot(m, aes(x=umap1, y=umap2, color=gini_index)) + 
  geom_point(size=0.1) + 
  scale_color_gradientn(colours = rainbow(3)) +
  labs(subtitle="ATAC-Seq peaks", 
       y="UMAP2", 
       x="UMAP1", 
       title="UMAP") + 
  theme_void()

ggsave("umap_gini_index.pdf")


# Promoter density UMAP
gg <- ggplot(m[m$peak_type=="Promoter",], aes(x=umap1, y=umap2)) + 
  stat_density_2d(aes(fill=..density..),geom="raster",contour=FALSE, n = 500, h=1) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  scale_fill_gradientn(colors = jdb_palette("horizon")) + 
  labs(subtitle="ATAC-Seq peaks - Promoters", 
       y="UMAP2", 
       x="UMAP1", 
       title="UMAP") +
  guides(colour = guide_legend(override.aes = list(size=2)))

ggsave("umap_promoters.pdf")

# Intronic density UMAP
gg <- ggplot(m[m$peak_type=="Intron",], aes(x=umap1, y=umap2)) + 
  stat_density_2d(aes(fill=..density..),geom="raster",contour=FALSE, n = 500, h=1) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  scale_fill_gradientn(colors = jdb_palette("horizon")) + 
  labs(subtitle="ATAC-Seq peaks - Introns", 
       y="UMAP2", 
       x="UMAP1", 
       title="UMAP") +
  guides(colour = guide_legend(override.aes = list(size=2)))

ggsave("umap_intron.pdf")

# Distal intergenic density UMAP
gg <- ggplot(m[m$peak_type=="Distal Intergenic",], aes(x=umap1, y=umap2)) + 
  stat_density_2d(aes(fill=..density..),geom="raster",contour=FALSE, n = 500, h=1) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  scale_fill_gradientn(colors = jdb_palette("horizon")) + 
  labs(subtitle="ATAC-Seq peaks - Distal intergenic", 
       y="UMAP2", 
       x="UMAP1", 
       title="UMAP") +
  guides(colour = guide_legend(override.aes = list(size=2)))

ggsave("umap_intergenic.pdf")
# ------------------------------------------------------------------------------------ #



# Figure S12 - Panel C ----------------------------------------------------------------
# Variance component decomposition of the mRNA expression for every gene (as column), 
# in a variance component model that discretizes the explanatory power of promoter OCRs (in green),
# distal enhancer OCRs (in blue), and unexplained variance (in red). 
# The limix python library was used for fitting the regression models (https://github.com/limix/limix).

# >>> Required packages
library(ggplot2)
library(ggsci)
library(tidyverse)
library(Ternary)

# >>> Input data
# Variance component model
sigmas2 <- read.csv(paste0(data_wd,"/osfstorage-archive/08_atac_rna_variance_decomposition/sigmas2.csv"), row.names = 1)
sigmas2_permuted <- read.csv(paste0(data_wd,"/osfstorage-archive/08_atac_rna_variance_decomposition/sigmas2_permuted.csv"), row.names = 1)

# DEGs by transitions
deg <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/differential_expression_analysis/rna_deg_clustering.rds"))

# >>> Computing
# removing the failed ones
toRemove <- which(sigmas2$pro_sigma2 == -1)
if(length(toRemove) > 0){
  print(paste0('To remove: ', length(toRemove)))
  sigmas2 <- sigmas2[-toRemove, ]
}

toRemove <- which(sigmas2_permuted$pro_sigma2 == -1)
if(length(toRemove) > 0){
  print(paste0('To remove: ', length(toRemove)))
  sigmas2_permuted <- sigmas2_permuted[-toRemove, ]
}

# computing the explained variance
explained_vars <- t(apply(sigmas2,1,FUN=function(x){return((x/sum(x))*100)}))
explained_vars <- data.frame(explained_vars,sig=rownames(explained_vars)%in%rownames(deg))
explained_vars_premuted <- t(apply(sigmas2_permuted,1,FUN=function(x){return((x/sum(x)))}))
explained_vars_premuted <- data.frame(explained_vars_premuted,sig=rownames(explained_vars_premuted)%in%rownames(deg))

# >>> Plotting
# Order according to enhancer
explained_vars_ord <- explained_vars[order(explained_vars$enh_sigma2, explained_vars$pro_sigma2, decreasing=TRUE),]
explained_vars_ord[which(explained_vars_ord$enh_sigma2<0.005),] <- explained_vars_ord[which(explained_vars_ord$enh_sigma2<0.005)[order(explained_vars_ord$pro_sigma2[which(explained_vars_ord$enh_sigma2<0.005)], decreasing=TRUE)],]

data_to_plot <- data.frame(value=c(explained_vars_ord$pro_sigma2,explained_vars_ord$enh_sigma2,explained_vars_ord$noise_sigma2),
                           type=factor(rep(c("Promoter","Enhancer","Noise"),each=nrow(explained_vars_ord)), levels=c("Noise","Promoter","Enhancer")),
                           rank=rep(1:nrow(explained_vars_ord),3))

p <- ggplot(data_to_plot, aes(x = rank, y = value, fill=type)) +
  geom_area(alpha=1) + 
  scale_fill_manual(values=c("#e41a1c","#4daf4a","#377eb8")) +
  scale_x_continuous(name = 'Genes') + 
  #scale_y_reverse() + 
  #scale_fill_npg() + 
  theme_bw()
ggsave(paste0(format(Sys.time(), "%Y%m%d"),"_explained_variance_ordered_by_enhancers.pdf"))


# Significant genes
explained_vars_sig_genes <- explained_vars[rownames(explained_vars)%in%rownames(deg),]

# order according to enhancer
explained_vars_ord <- explained_vars_sig_genes[order(explained_vars_sig_genes$enh_sigma2,explained_vars_sig_genes$pro_sigma2,decreasing=TRUE),]
explained_vars_ord[which(explained_vars_ord$enh_sigma2<0.005),] <- explained_vars_ord[which(explained_vars_ord$enh_sigma2<0.005)[order(explained_vars_ord$pro_sigma2[which(explained_vars_ord$enh_sigma2<0.005)],decreasing=TRUE)],]

data_to_plot <- data.frame(value=c(explained_vars_ord$pro_sigma2,explained_vars_ord$enh_sigma2,explained_vars_ord$noise_sigma2),
                           type=factor(rep(c("Promoter","Enhancer","Noise"),each=nrow(explained_vars_ord)), levels=c("Noise","Promoter","Enhancer")),
                           rank=rep(1:nrow(explained_vars_ord),3))

p <- ggplot(data_to_plot, aes(x = rank, y = value, fill=type)) +
  geom_area(alpha=1) + 
  scale_fill_manual(values=c("#e41a1c","#4daf4a","#377eb8")) +
  scale_x_continuous(name = 'Genes') + 
  #scale_y_reverse() + 
  #scale_fill_npg() + 
  theme_bw()
ggsave(paste0(format(Sys.time(), "%Y%m%d"),"_explained_variance_ordered_by_enhancers_degs.pdf"))


# No significant genes
explained_vars_nosig_genes <- explained_vars[!rownames(explained_vars)%in%rownames(deg),]

# order according to enhancer
explained_vars_ord <- explained_vars_nosig_genes[order(explained_vars_nosig_genes$enh_sigma2,
                                                       explained_vars_nosig_genes$pro_sigma2,decreasing=TRUE),]
explained_vars_ord[which(explained_vars_ord$enh_sigma2<0.005),] <- explained_vars_ord[which(explained_vars_ord$enh_sigma2<0.005)[order(explained_vars_ord$pro_sigma2[which(explained_vars_ord$enh_sigma2<0.005)],decreasing=TRUE)],]

data_to_plot <- data.frame(value=c(explained_vars_ord$pro_sigma2,explained_vars_ord$enh_sigma2,explained_vars_ord$noise_sigma2),
                           type=factor(rep(c("Promoter","Enhancer","Noise"),each=nrow(explained_vars_ord)), levels=c("Noise","Promoter","Enhancer")),
                           rank=rep(1:nrow(explained_vars_ord),3))

p <- ggplot(data_to_plot, aes(x = rank, y = value, fill=type)) +
  geom_area(alpha=1) + 
  scale_fill_manual(values=c("#e41a1c","#4daf4a","#377eb8")) +
  scale_x_continuous(name = 'Genes') + 
  #scale_y_reverse() + 
  #scale_fill_npg() + 
  theme_bw()
ggsave(paste0(format(Sys.time(), "%Y%m%d"),"_explained_variance_ordered_by_enhancers_no_degs.pdf"))


# Permuted values
# computing the explained variance order according to enhancer
explained_vars_ord_permuted <- explained_vars_premuted[order(explained_vars_premuted$enh_sigma2,explained_vars_premuted$pro_sigma2, decreasing=TRUE),]
explained_vars_ord_permuted[which(explained_vars_ord_permuted$enh_sigma2<0.001),] <- explained_vars_ord_permuted[which(explained_vars_ord_permuted$enh_sigma2<0.001)[order(explained_vars_ord_permuted$pro_sigma2[which(explained_vars_ord_permuted$enh_sigma2<0.001)], decreasing=TRUE)],]

data_to_plot <- data.frame(value=c(explained_vars_ord_permuted$pro_sigma2,explained_vars_ord_permuted$enh_sigma2,explained_vars_ord_permuted$noise_sigma2),
                           type=factor(rep(c("Promoter","Enhancer","Noise"),each=nrow(explained_vars_ord_permuted)), levels=c("Noise","Promoter","Enhancer")),
                           rank=rep(1:nrow(explained_vars_ord_permuted),3))

# plotting
p <- ggplot(data_to_plot, aes(x = rank, y = value, fill=type)) +
  geom_area(alpha=1) + 
  scale_fill_manual(values=c("#e41a1c","#4daf4a","#377eb8")) +
  scale_x_continuous(name = 'Genes') + 
  #scale_y_reverse() + 
  #scale_fill_npg() + 
  theme_bw()
ggsave(paste0(format(Sys.time(), "%Y%m%d"),"_explained_variance_ordered_by_enhancers_permuted.pdf"))
# ------------------------------------------------------------------------------------ #



# Figure S13 - Panel A ----------------------------------------------------------------
# Histone compartment distribution; overlap with out OCRs
# >>> Required packages
library(ggplot2)
library(GenomicRanges)

# >>> Input data
#HiC compartments
ocr_gene <- readRDS(paste0(data_wd,"/osfstorage-archive/02_cCREs_OCR_gene_links/public_epigenetic_data_for_validation/atac_OCR_annotation_extended.rds"))
ocr_gene$peak <- paste0(ocr_gene$chr,"_",ocr_gene$start,"_",ocr_gene$end)
# >>> Plotting
# Plot colors

# Plot
#Barplot output summary

data_to_plot <- ocr_gene[ocr_gene$annotation_simplified%in%c("Promoter","Intron","Distal Intergenic"),] #34,370

# Subset for All
p1 <- ggplot(data_to_plot, aes(x=annotation_simplified, fill=histone)) +
  geom_bar(position="fill") +
  theme_classic() +
  xlab("OCR")

# Subset for cCREs
p2 <- ggplot(data_to_plot[data_to_plot$CRE=="yes",], aes(x=annotation_simplified, fill=histone)) +
  geom_bar(position="fill") +
  theme_classic() +
  xlab("OCR")

# Subset for no cCREs
p3 <- ggplot(data_to_plot[data_to_plot$CRE=="no",], aes(x=annotation_simplified, fill=histone)) +
  geom_bar(position="fill") +
  theme_classic() +
  xlab("OCR")

p <- p1+p2+p3

ggsave(filename = paste(format(Sys.Date(),"%Y%m%d"),"_histone_marks.pdf",sep=""), plot = p)
# ------------------------------------------------------------------------------------ #



# Figure S13 - Panel B ----------------------------------------------------------------
# HiC compartment distribution; overlap with out OCRs
# >>> Required packages
library(ggplot2)
library(GenomicRanges)

# >>> Input data
#HiC compartments
ocr_gene <- readRDS(paste0(data_wd,"/osfstorage-archive/02_cCREs_OCR_gene_links/public_epigenetic_data_for_validation/atac_OCR_annotation_extended.rds"))
ocr_gene$peak <- paste0(ocr_gene$chr,"_",ocr_gene$start,"_",ocr_gene$end)
# >>> Plotting
# Plot colors

# Plot
#Barplot output summary

data_to_plot <- ocr_gene[ocr_gene$annotation_simplified%in%c("Promoter","Intron","Distal Intergenic"),] #34,370

# Subset for All

p1 <- ggplot(data_to_plot, aes(x=annotation_simplified, fill=HiC_compartment)) +
  geom_bar(position="fill") +
  theme_classic() +
  xlab("OCR")
  
# Subset for cCREs
p2 <- ggplot(data_to_plot[data_to_plot$CRE=="yes",], aes(x=annotation_simplified, fill=HiC_compartment)) +
  geom_bar(position="fill") +
  theme_classic() +
  xlab("OCR")

# Subset for no cCREs
p3 <- ggplot(data_to_plot[data_to_plot$CRE=="no",], aes(x=annotation_simplified, fill=HiC_compartment)) +
  geom_bar(position="fill") +
  theme_classic() +
  xlab("OCR")

p <- p1+p2+p3

ggsave(filename = paste(format(Sys.Date(),"%Y%m%d"),"_HiC.pdf",sep=""), plot = p)
# ------------------------------------------------------------------------------------ #



# Figure S13 - Panel C ----------------------------------------------------------------
# RNA and ATAC heatmap for PC-HiC validated cCREs

# >>> Required packages
library(ggplot2)
library(ggsci)
library(tidyverse)
library(Ternary)

# >>> Input data
# Gene Expression data
rna <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/rna_norm_data.rds"))
rna_metadata <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/rna_metadata.rds"))
rna_anno <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/rna_gene_annotation.rds"))

# Chromatin accessibility data
atac <- readRDS(paste0(data_wd,"/osfstorage-archive/01_atac_seq_data/atac_norm_data.rds"))
atac_metadata <- readRDS(paste0(data_wd,"/osfstorage-archive/01_atac_seq_data/atac_metadata.rds"))
atac_anno <- readRDS(paste0(data_wd,"/osfstorage-archive/01_atac_seq_data/atac_OCR_annotation.rds"))

# GRN data
tf_ocr_gene <- read.table(paste0(data_wd,"/osfstorage-archive/04_gene_regulatory_networks/tf_ocr_gene_interactions_curated.txt"),
                          sep="\t",
                          dec=",",
                          header=TRUE)

ocr_gene <- readRDS(paste0(data_wd,"/osfstorage-archive/02_cCREs_OCR_gene_links/public_epigenetic_data_for_validation/promoter_capture_HiC/candidate_CREs_plus_PCHiC.rds"))

# >>> Plotting
# Plot colors
bcell_colors <- readRDS(paste0(data_wd,"/osfstorage-archive/00_plot_parameters/bcell_colors.rds"))
bcell_colors <- readRDS(paste0(data_wd,"/osfstorage-archive/00_plot_parameters/bcell_colors.rds"))
my_colors_rna <- readRDS(paste0(data_wd,"/osfstorage-archive/00_plot_parameters/rna_expression_colors.rds"))
my_colors_atac <- readRDS(paste0(data_wd,"/osfstorage-archive/00_plot_parameters/atac_colors.rds"))


enhancer_ocr_gene <- ocr_gene[ocr_gene$overlap_ER_Biola_ourgene_in_baitName=="TRUE",] #2,161
rownames(enhancer_ocr_gene) <- paste0(enhancer_ocr_gene$hgnc_symbol," : ",enhancer_ocr_gene$peak)

#RNA data
data_to_plot <- t(scale(t(rna[enhancer_ocr_gene$gene,]),center=TRUE, scale=TRUE))
rownames(data_to_plot) <- rownames(enhancer_ocr_gene)

ha = HeatmapAnnotation(stage = rna_metadata$CellType,
                       col = list(stage = bcell_colors))


target_genes <- c("PAX5","MS4A1","CD79A","CD79B","IKZF1","TCF3","EBF1")

sel_pos <- vector()
target_genes2 <- vector()
for(i in 1:length(target_genes)){
  ss <- grep(target_genes[i],rownames(data_to_plot))
  if(length(ss)>0){
    sel_pos <- c(sel_pos,ss)
    target_genes2 <- c(target_genes2,rownames(data_to_plot)[ss])
  }
}

rha_text <- rowAnnotation(genes=anno_mark(at = sel_pos, labels = target_genes2))

exprs_pattern <- Heatmap(data_to_plot, name = "Expression patterns",
                         top_annotation = ha, column_order=order(rna_metadata$CellType), 
                         right_annotation = rha_text, 
                         col=colorRamp2(seq(-3,3, length.out=length(my_colors_rna)),my_colors_rna),
                         show_column_names = FALSE, show_row_names = FALSE, 
                         cluster_rows = TRUE,
                         cluster_columns = FALSE, use_raster=TRUE, width = unit(4, "cm"))


#ATAC data
data_to_plot <- t(scale(t(atac[enhancer_ocr_gene$peak,]),center=TRUE, scale=TRUE))
rownames(data_to_plot) <- rownames(enhancer_ocr_gene)

ha = HeatmapAnnotation(stage = atac_metadata$CellType,
                       col = list(stage = bcell_colors))

atac_pattern <- Heatmap(data_to_plot, name = "Expression patterns",
                        top_annotation = ha, column_order=order(atac_metadata$CellType), 
                        col=colorRamp2(seq(-3,3, length.out=length(my_colors_atac)),my_colors_atac),
                        show_column_names = FALSE, show_row_names = FALSE, 
                        cluster_rows = TRUE,
                        cluster_columns = FALSE, use_raster=TRUE, width = unit(4, "cm"))


pdf(paste(format(Sys.time(), "%Y%m%d"),"_heatmap_cCREs_validated_by_PCHiC.pdf",sep=""), width=12)
atac_pattern + exprs_pattern
dev.off()
# ------------------------------------------------------------------------------------ #




# Figure S14 ---------------------------------------------------------------------

# B-ALL Public datasets analysis -------------------------------------------------
###-----------------------------------------------------------------------------------###
# ****** Zinngrebe dataset (GSE140556) ****** #
###  Zinngrebe dataset ----------
# Biomarker Profile for Prediction of Response to SMAC Mimetic Monotherapy in Pediatric Precursor B-Cell Acute Lymphoblastic Leukemia
# Second mitochondria-derived activator of caspase (SMAC) mimetics (SMs) targeting inhibitor of apoptosis proteins (IAPs) activate cell death pathways.
# Microarray data from 8 SM-sensitive and 6 SM-insensitive BCP-ALL samples.

library(Biobase)
library(GEOquery)

# >>> Input data
# B-ALL Gene Expression data
# Get the data from GEO database
gse <- getGEO("GSE140556", GSEMatrix = TRUE)
show(gse)
gse <- gse[[1]]

ball_norm_rna_Zinngrebe <- exprs(gse)
ball_metadata_Zinngrebe <- pData(gse)

# Translate probe ID to gene names:
## Exclude genes with multiple annotations
gene_name_Zinngrebe <- fData(gse)[,"Gene Symbol"]
to_exclude <- grep("/",gene_name_Zinngrebe)
length(to_exclude)
ball_norm_rna_Zinngrebe <- ball_norm_rna_Zinngrebe[-to_exclude,]
rownames(ball_norm_rna_Zinngrebe) <- gene_name_Zinngrebe[-to_exclude]


# >>> Exploratory analysis
# Principal component analysis
pca_output <- prcomp(t(ball_norm_rna_Zinngrebe))
pca_output_summary <- summary(pca_output)

data_to_plot <- data.frame(PC1=pca_output$x[,1],PC2=pca_output$x[,2], 
                           Response=ball_metadata_Zinngrebe$"sm response:ch1")

# Plotting
# Plot colors
response_color <- c("SM-insensitive"="#bebada","SM-sensitive"="#fb8072")

# Plot
p_response <- ggplot(data_to_plot, aes(x = PC1, y = PC2, color = Response)) +
  geom_point(size=3) +
  scale_color_manual(values=response_color) +
  labs(title="B-ALL Zinngrebe dataset", x = paste0("PC1 (",round(pca_output_summary$importance[2,1]*100,1),")%"), y = paste0("PC2 (",round(pca_output_summary$importance[2,2]*100,1),")%")) +
  theme_classic()


pdf(paste(format(Sys.time(), "%Y%m%d"),"_pca_BALL_Zinngrebe.pdf",sep=""), width=6, height = 6)
print(p_response)
print(p_response+theme(legend.position = "none"))
dev.off()


# >>> Computing GSVA
# GRN data
tf_ocr_gene <- read.table(paste0(data_wd,"/osfstorage-archive/04_gene_regulatory_networks/tf_ocr_gene_interactions_curated.txt"),
                          sep="\t",
                          dec=",",
                          header=TRUE)

### GSVA
# Get the TF regulons gene sets for the positively regulated genes.
df_regulons_pos <- unique(tf_ocr_gene[tf_ocr_gene$likelihood=="activator",c("tf_symbol","gene_symbol")])

df_regulons_list_pos <- list()
for(i in unique(df_regulons_pos$tf_symbol)){
  df_regulons_list_pos[[i]] <- df_regulons_pos$gene_symbol[df_regulons_pos$tf_symbol==i]
}

gs_symbol <- df_regulons_list_pos
gsvaPar <- gsvaParam(ball_norm_rna_Zinngrebe, gs_symbol, kcdf="Gaussian")
gsva.BALL_Zinngrebe <- gsva(gsvaPar, verbose=FALSE)


# >>> Save results
save(ball_norm_rna_Zinngrebe,ball_metadata_Zinngrebe,gsva.BALL_Zinngrebe,file="BALL_Zinngrebe.RData")
###-----------------------------------------------------------------------------------###


###-----------------------------------------------------------------------------------###
# ****** Singh dataset (GSE153670) ****** #
###  Singh dataset ----------
# Analyzed sequencing data from bone marrow samples taken from patients enrolled on ENSIGN or ELIANA clinical trials prior to infusion of CART19.

library(Biobase)
library(GEOquery)

# >>> Input data
# B-ALL Gene Expression data
# Get the data from GEO database
gse <- getGEO("GSE153670", GSEMatrix = TRUE)
show(gse)
gse <- gse[[1]]
#getGEOSuppFiles("GSE153670")
ball_rna_Singh <- read.csv("GSE153670/GSE153670_Tisagenlecleucel_patient_RNAseq_Count_matrix.csv", sep=";", check.names = FALSE, row.names = 1)
ball_metadata_Singh <- pData(gse)


# >>>  Data normalization
require(edgeR)
# Create DEGlist object with count data and defining the grouping of the variables
dge <- DGEList(counts=ball_rna_Singh, group=ball_metadata_Singh$"patient response:ch1")

# Filter out low expressed genes
filtering_criteria_edgeR <- filterByExpr(dge, group = ball_metadata_Singh$"patient response:ch1",  min.count = 5, min.total.count = 10)
dge.filt <- dge[filtering_criteria_edgeR,]
dim(dge.filt)

# Define the model matrix
dge.filt <- calcNormFactors(dge.filt) # computes the correction factor based on TMM (default method)
group <- as.factor(ball_metadata_Singh$"patient response:ch1")
design <- model.matrix( ~ 0 + group)
colnames(design) <- gsub("group","", colnames(design))

# Voom transformation including log2 CPM normalization
v <- voom(dge.filt, design = design, plot = TRUE)
ball_norm_rna_Singh <- v$E

# >>> Exploratory analysis
# Principal component analysis
pca_output <- prcomp(t(ball_norm_rna_Singh))
pca_output_summary <- summary(pca_output)

data_to_plot <- data.frame(PC1=pca_output$x[,1],PC2=pca_output$x[,2], 
                           Response=ball_metadata_Singh$"patient response:ch1")

# Plotting
# Plot colors
response_color <- c("NR"="#bebada","CR"="#fb8072")


# Plot
p_response <- ggplot(data_to_plot, aes(x = PC1, y = PC2, color = Response)) +
  geom_point(size=3) +
  scale_color_manual(values=response_color) +
  labs(title="B-ALL Singh multi-dataset", x = paste0("PC1 (",round(pca_output_summary$importance[2,1]*100,1),")%"), y = paste0("PC2 (",round(pca_output_summary$importance[2,2]*100,1),")%")) +
  theme_classic()


pdf(paste(format(Sys.time(), "%Y%m%d"),"_pca_BALL_Singh.pdf",sep=""), width=6, height = 6)
print(p_response)
print(p_response+theme(legend.position = "none"))
dev.off()


# >>> Computing GSVA
# GRN data
tf_ocr_gene <- read.table(paste0(data_wd,"/osfstorage-archive/04_gene_regulatory_networks/tf_ocr_gene_interactions_curated.txt"),
                          sep="\t",
                          dec=",",
                          header=TRUE)

### GSVA
# Get the TF regulons gene sets for the positively regulated genes.
df_regulons_pos <- unique(tf_ocr_gene[tf_ocr_gene$likelihood=="activator",c("tf_symbol","gene_symbol")])

df_regulons_list_pos <- list()
for(i in unique(df_regulons_pos$tf_symbol)){
  df_regulons_list_pos[[i]] <- df_regulons_pos$gene_symbol[df_regulons_pos$tf_symbol==i]
}

gs_symbol <- df_regulons_list_pos
gsvaPar <- gsvaParam(ball_norm_rna_Singh, gs_symbol, kcdf="Gaussian")
gsva.BALL_Singh <- gsva(gsvaPar, verbose=FALSE)


# >>> Save results
save(ball_norm_rna_Singh,ball_metadata_Singh,gsva.BALL_Singh,file="BALL_Singh.RData")
###-----------------------------------------------------------------------------------###


# Figure S14 - Panel A,B ----------------------------------------------------------------
###-----------------------------------------------------------------------------------###
# ****** Singh dataset (GSE153670) ****** #
###  Singh dataset ----------
# Bulk RNA sequencing from patient bone marrow samples

library(Biobase)
library(GEOquery)

# >>> Input data
# B-ALL Gene Expression data
# B-ALL Singh Gene Expression data
load("BALL_Singh.RData", verbose = TRUE)

# >>> Exploratory analysis
# Principal component analysis
pca_output <- prcomp(t(ball_norm_rna_Singh))
pca_output_summary <- summary(pca_output)

data_to_plot <- data.frame(PC1=pca_output$x[,1],PC2=pca_output$x[,2], Response=ball_metadata_Singh$"patient response:ch1")

# Plotting
# Plot colors
response_color <- c("NR"="#bebada","CR"="#fb8072")


# Plot
p_response <- ggplot(data_to_plot, aes(x = PC1, y = PC2, color = Response)) +
  geom_point(size=3) +
  scale_color_manual(values=response_color) +
  labs(title="B-ALL Singh dataset", x = paste0("PC1 (",round(pca_output_summary$importance[2,1]*100,1),"%)"), y = paste0("PC2 (",round(pca_output_summary$importance[2,2]*100,1),"%)")) +
  theme_classic()

pdf(paste(format(Sys.time(), "%Y%m%d"),"_pca_BALL_Singh.pdf",sep=""), width=6, height = 6)
  print(p_response)
  print(p_response+theme(legend.position = "none"))
dev.off()


# >>> Computing GSVA
# GRN data
tf_ocr_gene <- read.table(paste0(data_wd,"/osfstorage-archive/04_gene_regulatory_networks/tf_ocr_gene_interactions_curated.txt"),
                          sep="\t",
                          dec=",",
                          header=TRUE)

### GSVA
# Get the TF regulons gene sets for the positively regulated genes.
df_regulons_pos <- unique(tf_ocr_gene[tf_ocr_gene$likelihood=="activator",c("tf_symbol","gene_symbol")])

df_regulons_list_pos <- list()
for(i in unique(df_regulons_pos$tf_symbol)){
  df_regulons_list_pos[[i]] <- df_regulons_pos$gene_symbol[df_regulons_pos$tf_symbol==i]
}

gs_symbol <- df_regulons_list_pos
gsvaPar <- gsvaParam(ball_norm_rna_Singh, gs_symbol, kcdf="Gaussian")
gsva.BALL_Singh <- gsva(gsvaPar, verbose=FALSE)


# >>> Save results
save(ball_norm_rna_Singh,ball_metadata_Singh,gsva.BALL_Singh,file="BALL_Singh.RData")


# >>> Heatmap plot
# Plot  
# Get the regulons classification
regulons_group <- readRDS(paste0(data_wd,"/osfstorage-archive/04_gene_regulatory_networks/regulons/regulons_clusters.rds"))
regulons_group_plot <- c(rep("1",length(regulons_group[[1]])),rep("2",length(regulons_group[[2]])),rep("3",length(regulons_group[[3]])),rep("4",length(regulons_group[[4]])),rep("5",length(regulons_group[[5]])))
names(regulons_group_plot) <- c(regulons_group[[1]],regulons_group[[2]],regulons_group[[3]],regulons_group[[4]],regulons_group[[5]])
# Plot colors
cols <- jdb_palette("brewer_spectra")
col_fun <- colorRamp2(seq(-2, 2, length.out = length(cols)), cols)
regulons_color <- c("1"="#92c5de","2"="#2166ac","3"="#91cf60","4"="#dd1c77","5"="black")
response_color <- c("NR"="#bebada","CR"="#fb8072")

# Heatmap
set.seed(123)
h_scaled_BALL_Singh <- Heatmap(t(scale(t(gsva.BALL_Singh))),
                                   top_annotation = HeatmapAnnotation("Response" = ball_metadata_Singh$`patient response:ch1`,
                                                                      col=list("Response"= response_color)),
                                   left_annotation = rowAnnotation("TF_subgroup" = regulons_group_plot[rownames(gsva.BALL_Singh)],
                                                                   col=list("TF_subgroup" = regulons_color)),
                                   col = col_fun,
                                   column_km = 2,
                                   row_split = factor(regulons_group_plot[rownames(gsva.BALL_Singh)],levels=c("1","2","3","4","5")),
                                   column_gap = unit(0.5, "mm"),
                                   row_gap = unit(0.5, "mm"),
                                   show_column_names = FALSE,
                                   show_row_names = FALSE,
                                   row_title = NULL,
                                   column_title = NULL,
                                   row_names_gp = gpar(fontsize = 6),
                                   name="GSVA score (scaled)"
)

## Get the sample clusterng
column_clusters <- column_order(h_scaled_BALL_Singh)
cluster_assignments <- rep(NA, ncol(gsva.BALL_Singh))
names(cluster_assignments) <- colnames(gsva.BALL_Singh)

for (i in seq_along(column_clusters)) {
  cluster_assignments[column_clusters[[i]]] <- i
}

head(cluster_assignments)


h_scaled_BALL_Singh_clust <- Heatmap(t(scale(t(gsva.BALL_Singh))),
                               top_annotation = HeatmapAnnotation("km.Clustering" = cluster_assignments,"Response" = ball_metadata_Singh$`patient response:ch1`,
                                                                  col=list("Response"= response_color,
                                                                           "km.Clustering" = c("1"="#DAA520","2"="#2F4F4F"))),
                               left_annotation = rowAnnotation("TF_subgroup" = regulons_group_plot[rownames(gsva.BALL_Singh)],
                                                               col=list("TF_subgroup" = regulons_color)),
                               col = col_fun,
                               column_km = 2,
                               row_split = factor(regulons_group_plot[rownames(gsva.BALL_Singh)],levels=c("1","2","3","4","5")),
                               column_gap = unit(0.5, "mm"),
                               row_gap = unit(0.5, "mm"),
                               show_column_names = FALSE,
                               show_row_names = FALSE,
                               row_title = NULL,
                               column_title = NULL,
                               row_names_gp = gpar(fontsize = 6),
                               name="GSVA score (scaled)"
)


pdf(paste0(format(Sys.time(), "%Y%m%d"),"_heatmap_GSVA_regulons_BALL_Singh.pdf"), width = 8, heigh=5)
  h_scaled_BALL_Singh_clust
dev.off()

## Differential expression at regulon level
library(sva)
library(limma)

group <- ball_metadata_Singh$`patient response:ch1`
## build design matrix of the model to which we fit the data
mod <- model.matrix(~ group, as.data.frame(colnames(gsva.BALL_Singh)))
## build design matrix of the corresponding null model
mod0 <- model.matrix(~ 1, as.data.frame(colnames(gsva.BALL_Singh)))
## estimate surrogate variables (SVs) with SVA
sv <- sva(gsva.BALL_Singh, mod, mod0)
## add SVs to the design matrix of the model of interest
mod <- cbind(mod, sv$sv)
## fit linear models
fit <- lmFit(gsva.BALL_Singh, mod)
## calculate moderated t-statistics using the robust regime
fit.eb <- eBayes(fit, robust=TRUE)
## summarize the extent of differential expression at 5% FDR
res <- decideTests(fit.eb)
summary(res)
tt <- topTable(fit.eb, coef=2, n=Inf)

DEpwys <- rownames(tt)[tt$P.Val <= 0.05]
length(DEpwys)
h_scaled_BALL_Singh_DEGs <- Heatmap(t(scale(t(gsva.BALL_Singh[DEpwys,]))),
                                        top_annotation = HeatmapAnnotation("Response" = ball_metadata_Singh$`patient response:ch1`,
                                                                           col=list("Response"= response_color)),
                                        left_annotation = rowAnnotation("TF_subgroup" = regulons_group_plot[DEpwys],
                                                                        col=list("TF_subgroup" = regulons_color)),
                                        col = col_fun,
                                        column_km = 2,
                                        row_split = factor(regulons_group_plot[DEpwys],levels=c("1","2","3","4","5")),
                                        column_gap = unit(0.5, "mm"),
                                        row_gap = unit(0.5, "mm"),
                                        show_column_names = FALSE,
                                        show_row_names = TRUE,
                                        row_title = NULL,
                                        column_title = NULL,
                                        row_names_gp = gpar(fontsize = 8),
                                        name="GSVA score (scaled)"
)

pdf(paste0(format(Sys.time(), "%Y%m%d"),"_heatmap_GSVA_regulons_BALL_Singh_sig.pdf"), width = 8, heigh=6)
h_scaled_BALL_Singh_DEGs
dev.off()
###-----------------------------------------------------------------------------------###



# Figure S14 - Panel C,D ----------------------------------------------------------------
###-----------------------------------------------------------------------------------###
# ****** Zinngrebe dataset (GSE140556) ****** #
###  Zinngrebe dataset ----------
# microarray data from 8 SM-sensitive and 6 SM-insensitive BCP-ALL samples

# >>> Input data
# B-ALL Zinngrebe Gene Expression data
load("BALL_Zinngrebe.RData", verbose = TRUE)

# >>> Ploting
# Plot colors
response_color <- c("SM-insensitive"="#bebada","SM-sensitive"="#fb8072")

# Plot 
### Exploratory analysis: PCA

# >>> Exploratory analysis
# Principal component analysis
pca_output <- prcomp(t(ball_norm_rna_Zinngrebe))
pca_output_summary <- summary(pca_output)

data_to_plot <- data.frame(PC1=pca_output$x[,1],PC2=pca_output$x[,2], Response=ball_metadata_Zinngrebe$`sm response:ch1`)

# Plot
p_response <- ggplot(data_to_plot, aes(x = PC1, y = PC2, color = Response)) +
  geom_point(size=3) +
  scale_color_manual(values=response_color) +
  labs(title="B-ALL Zinngrebe dataset", x = paste0("PC1 (",round(pca_output_summary$importance[2,1]*100,1),"%)"), y = paste0("PC2 (",round(pca_output_summary$importance[2,2]*100,1),"%)")) +
  theme_classic()

pdf(paste(format(Sys.time(), "%Y%m%d"),"_pca_BALL_Zinngrebe.pdf",sep=""), width=6, height = 6)
  print(p_response)
  print(p_response+theme(legend.position = "none"))
dev.off()


# >>> Heatmap GSVA
# Plot  
# Get the regulons classification
regulons_group <- readRDS(paste0(data_wd,"/osfstorage-archive/04_gene_regulatory_networks/regulons/regulons_clusters.rds"))
regulons_group_plot <- c(rep("1",length(regulons_group[[1]])),rep("2",length(regulons_group[[2]])),rep("3",length(regulons_group[[3]])),rep("4",length(regulons_group[[4]])),rep("5",length(regulons_group[[5]])))
names(regulons_group_plot) <- c(regulons_group[[1]],regulons_group[[2]],regulons_group[[3]],regulons_group[[4]],regulons_group[[5]])
# Plot colors
cols <- jdb_palette("brewer_spectra")
col_fun <- colorRamp2(seq(-2, 2, length.out = length(cols)), cols)
regulons_color <- c("1"="#92c5de","2"="#2166ac","3"="#91cf60","4"="#dd1c77","5"="black")
response_color <- c("SM-insensitive"="#bebada","SM-sensitive"="#fb8072")

# Heatmap
set.seed(123)
h_scaled_BALL_Zinngrebe <- Heatmap(t(scale(t(gsva.BALL_Zinngrebe))),
                               top_annotation = HeatmapAnnotation("Response" = ball_metadata_Zinngrebe$`sm response:ch1`,
                                                                  col=list("Response"= response_color)),
                               left_annotation = rowAnnotation("TF_subgroup" = regulons_group_plot[rownames(gsva.BALL_Zinngrebe)],
                                                               col=list("TF_subgroup" = regulons_color)),
                               col = col_fun,
                               column_km = 2,
                               row_split = factor(regulons_group_plot[rownames(gsva.BALL_Zinngrebe)],levels=c("1","2","3","4","5")),
                               column_gap = unit(0.5, "mm"),
                               row_gap = unit(0.5, "mm"),
                               show_column_names = FALSE,
                               show_row_names = FALSE,
                               row_title = NULL,
                               column_title = NULL,
                               row_names_gp = gpar(fontsize = 6),
                               name="GSVA score (scaled)"
)


## Get the sample clusterng
column_clusters <- column_order(h_scaled_BALL_Zinngrebe)
cluster_assignments <- rep(NA, ncol(gsva.BALL_Zinngrebe))
names(cluster_assignments) <- colnames(gsva.BALL_Zinngrebe)

for (i in seq_along(column_clusters)) {
  cluster_assignments[column_clusters[[i]]] <- i
}

head(cluster_assignments)


h_scaled_BALL_Zinngrebe_clust <- Heatmap(t(scale(t(gsva.BALL_Zinngrebe))),
                                     top_annotation = HeatmapAnnotation("km.Clustering" = cluster_assignments,"Response" = ball_metadata_Zinngrebe$`sm response:ch1`,
                                                                        col=list("Response"= response_color,
                                                                                 "km.Clustering" = c("1"="#DAA520","2"="#2F4F4F"))),
                                     left_annotation = rowAnnotation("TF_subgroup" = regulons_group_plot[rownames(gsva.BALL_Zinngrebe)],
                                                                     col=list("TF_subgroup" = regulons_color)),
                                     col = col_fun,
                                     column_km = 2,
                                     row_split = factor(regulons_group_plot[rownames(gsva.BALL_Zinngrebe)],levels=c("1","2","3","4","5")),
                                     column_gap = unit(0.5, "mm"),
                                     row_gap = unit(0.5, "mm"),
                                     show_column_names = FALSE,
                                     show_row_names = FALSE,
                                     row_title = NULL,
                                     column_title = NULL,
                                     row_names_gp = gpar(fontsize = 6),
                                     name="GSVA score (scaled)"
)


pdf(paste0(format(Sys.time(), "%Y%m%d"),"_heatmap_GSVA_regulons_BALL_Zinngrebe.pdf"), width = 8, heigh=5)
h_scaled_BALL_Zinngrebe_clust
dev.off()


## Differential expression at regulon level
library(sva)
library(limma)

group <- ball_metadata_Zinngrebe$`sm response:ch1`
## build design matrix of the model to which we fit the data
mod <- model.matrix(~ group, as.data.frame(colnames(gsva.BALL_Zinngrebe)))
## build design matrix of the corresponding null model
mod0 <- model.matrix(~ 1, as.data.frame(colnames(gsva.BALL_Zinngrebe)))
## estimate surrogate variables (SVs) with SVA
sv <- sva(gsva.BALL_Zinngrebe, mod, mod0)
## add SVs to the design matrix of the model of interest
mod <- cbind(mod, sv$sv)
## fit linear models
fit <- lmFit(gsva.BALL_Zinngrebe, mod)
## calculate moderated t-statistics using the robust regime
fit.eb <- eBayes(fit, robust=TRUE)
## summarize the extent of differential expression at 5% FDR
res <- decideTests(fit.eb)
summary(res)
tt <- topTable(fit.eb, coef=2, n=Inf)

DEpwys <- rownames(tt)[tt$P.Val <= 0.05]
length(DEpwys)
h_scaled_BALL_Zinngrebe_DEGs <- Heatmap(t(scale(t(gsva.BALL_Zinngrebe[DEpwys,]))),
                                   top_annotation = HeatmapAnnotation("Response" = ball_metadata_Zinngrebe$`sm response:ch1`,
                                                                      col=list("Response"= response_color)),
                                   left_annotation = rowAnnotation("TF_subgroup" = regulons_group_plot[DEpwys],
                                                                   col=list("TF_subgroup" = regulons_color)),
                                   col = col_fun,
                                   column_km = 2,
                                   row_split = factor(regulons_group_plot[DEpwys],levels=c("1","2","3","4","5")),
                                   column_gap = unit(0.5, "mm"),
                                   row_gap = unit(0.5, "mm"),
                                   show_column_names = FALSE,
                                   show_row_names = TRUE,
                                   row_title = NULL,
                                   column_title = NULL,
                                   row_names_gp = gpar(fontsize = 8),
                                   name="GSVA score (scaled)"
)

pdf(paste0(format(Sys.time(), "%Y%m%d"),"_heatmap_GSVA_regulons_BALL_Zinngrebe_sig.pdf"), width = 8, heigh=4)
  h_scaled_BALL_Zinngrebe_DEGs
dev.off()
###-----------------------------------------------------------------------------------###


# Table S6 -----------------------------------------------------------------
# Pathway enrichment analysis of EBF1-regulon. 

# >>> Required packages
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)

# >>> Input data
# Gene Expression data
rna <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/rna_norm_data.rds"))
# GRN data
tf_ocr_gene <- read.table(paste0(data_wd,"/osfstorage-archive/04_gene_regulatory_networks/tf_ocr_gene_interactions_curated.txt"),
                          sep="\t",
                          dec=",",
                          header=TRUE)

# >>> Computing
#Genes regulated by EBF1 according to our network.
target_gene <- "EBF1"
list_target_genes_EBF1 <- unique(tf_ocr_gene[which(tf_ocr_gene$tf_symbol==target_gene),])

# Split by positive or negative regulation
list_target_genes_EBF1_pos <- unique(list_target_genes_EBF1[list_target_genes_EBF1$gene_tf_rho>0,"gene_ensembl"])

list_target_genes_EBF1_neg <- unique(list_target_genes_EBF1[list_target_genes_EBF1$gene_tf_rho<0,"gene_ensembl"])

# ORA based on GO terms with all the genes linked to EBF1
ora_GO_EBF1 <- enrichGO(gene          = list_target_genes_EBF1$gene_ensembl,
                        universe      = rownames(rna),
                        OrgDb         = org.Hs.eg.db,
                        keyType       = 'ENSEMBL',
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        minGSSize = 20,
                        maxGSSize = 500,
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05,
                        readable      = TRUE)

# ORA based on GO terms with all the genes linked to EBF1 positively
ora_GO_EBF1_pos <- enrichGO(gene      = list_target_genes_EBF1_pos,
                            universe      = rownames(rna),
                            OrgDb         = org.Hs.eg.db,
                            keyType       = 'ENSEMBL',
                            ont           = "BP",
                            pAdjustMethod = "BH",
                            minGSSize = 20,
                            maxGSSize = 500,
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 0.05,
                            readable      = TRUE)

# ORA based on GO terms with all the genes linked to EBF1 negatively
ora_GO_EBF1_neg <- enrichGO(gene          = list_target_genes_EBF1_neg,
                            universe      = rownames(rna),
                            OrgDb         = org.Hs.eg.db,
                            keyType       = 'ENSEMBL',
                            ont           = "BP",
                            pAdjustMethod = "BH",
                            minGSSize = 20,
                            maxGSSize = 500,
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 0.05,
                            readable      = TRUE)
# ------------------------------------------------------------------------------------ #



# Table S8, S9 -----------------------------------------------------------------
# Pathway enrichment analysis of ELK3-regulon-activated and ELK3-regulon-repressed genes. 

# >>> Required packages
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)

# >>> Input data
# Gene Expression data
rna <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/rna_norm_data.rds"))
# GRN data
tf_ocr_gene <- read.table(paste0(data_wd,"/osfstorage-archive/04_gene_regulatory_networks/tf_ocr_gene_interactions_curated.txt"),
                          sep="\t",
                          dec=",",
                          header=TRUE)

# >>> Computing
#Genes regulated by ELK3 according to our network.
target_gene <- "ELK3"
list_target_genes_ELK3 <- unique(tf_ocr_gene[which(tf_ocr_gene$tf_symbol==target_gene),])

# Split by positive or negative regulation
list_target_genes_ELK3_pos <- unique(list_target_genes_ELK3[list_target_genes_ELK3$gene_tf_rho>0,"gene_ensembl"])

list_target_genes_ELK3_neg <- unique(list_target_genes_ELK3[list_target_genes_ELK3$gene_tf_rho<0,"gene_ensembl"])

# ORA based on GO terms with all the genes linked to ELK3
ora_GO_ELK3 <- enrichGO(gene          = list_target_genes_ELK3$gene_ensembl,
                        universe      = rownames(rna),
                        OrgDb         = org.Hs.eg.db,
                        keyType       = 'ENSEMBL',
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        minGSSize = 20,
                        maxGSSize = 500,
                        pvalueCutoff  = 0.05,
                        qvalueCutoff  = 0.05,
                        readable      = TRUE)

# ORA based on GO terms with all the genes linked to ELK3 positively
ora_GO_ELK3_pos <- enrichGO(gene      = list_target_genes_ELK3_pos,
                            universe      = rownames(rna),
                            OrgDb         = org.Hs.eg.db,
                            keyType       = 'ENSEMBL',
                            ont           = "BP",
                            pAdjustMethod = "BH",
                            minGSSize = 20,
                            maxGSSize = 500,
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 0.05,
                            readable      = TRUE)

# ORA based on GO terms with all the genes linked to ELK3 negatively
ora_GO_ELK3_neg <- enrichGO(gene          = list_target_genes_ELK3_neg,
                            universe      = rownames(rna),
                            OrgDb         = org.Hs.eg.db,
                            keyType       = 'ENSEMBL',
                            ont           = "BP",
                            pAdjustMethod = "BH",
                            minGSSize = 20,
                            maxGSSize = 500,
                            pvalueCutoff  = 0.05,
                            qvalueCutoff  = 0.05,
                            readable      = TRUE)
# ------------------------------------------------------------------------------------ #



