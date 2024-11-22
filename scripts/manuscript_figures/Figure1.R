#***************************************************
## Paper: "Uncovering the regulatory landscape of early human B-cell lymphopoiesis and its 
## implications in the pathogenesis of B-cell acute lymphoblastic leukemia"
## Authors: Planell N et al.
## Date: 2024
## Code: Figure 1
## Input data available at: https://osf.io/gswpy/?view_only=39cfd50d5a7e44e5a2d141869dcb1315
#***************************************************

# Define the directory were input data from OSF has been downloaded
data_wd <- setwd("../GitHub_code/")

# ****** Figure 1 ****** #

# Figure 1 - Panel A -----------------------------------------------------------------
# B-cell differentiation scheme.
# ------------------------------------------------------------------------------------ #



# Figure 1 - Panel B -----------------------------------------------------------------
# Principal Component Analysis (PCA) of gene expression (79 samples).

# >>> Required packages
library(ggplot2)
library(BuenColors)

# >>> Input data
# Gene Expression data
rna <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/rna_norm_data.rds"))
rna_metadata <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/rna_metadata.rds"))

# >>> Computation
# PCA computation
pca_rna <- prcomp(t(rna))
var_prcomp <- pca_rna$sdev^2

# PCs variance
pcvar <- data.frame(var=var_prcomp/sum(var_prcomp), 
                    pc=c(1:length(var_prcomp)))

# >>> Plotting
# Plot colors
bcell_colors <- readRDS(paste0(data_wd,"/osfstorage-archive/00_plot_parameters/bcell_colors.rds"))

# Plot
data_to_plot <- data.frame(PC1=pca_rna$x[,1],PC2=pca_rna$x[,2], CellType=rna_metadata$CellType)

p <- ggplot(data_to_plot, aes(x = PC1, y = PC2, color = CellType)) +
  geom_point(size = 3, show.legend=FALSE) + 
  coord_fixed(ratio=1) + 
  scale_color_manual(values = bcell_colors) + 
  directlabels::geom_dl(aes(label = CellType), method = list("smart.grid", cex=1)) + 
  labs(title="PCA: RNA-Seq",x = paste("PC1 (",round(pcvar[1,1]*100,0),"%)", sep=""), y = paste("PC2 (",round(pcvar[2,1]*100,0),"%)", sep=""))+
  theme_classic()

pdf(paste(format(Sys.time(), "%Y%m%d"),"_pca_rna.pdf",sep=""), width=6, height = 6)
print(p)
dev.off()
# ------------------------------------------------------------------------------------ #



# Figure 1 - Panel C -----------------------------------------------------------------
# Principal Component Analysis (PCA) of chromatin accessibility (78 samples).

# >>> Required packages
library(ggplot2)
library(BuenColors)

# >>> Input data
# Chromatin accessibility data
atac <- readRDS(paste0(data_wd,"/osfstorage-archive/01_atac_seq_data/atac_norm_data.rds"))
atac_metadata <- readRDS(paste0(data_wd,"/osfstorage-archive/01_atac_seq_data/atac_metadata.rds"))

# >>> Computation
# PCA computation
pca_atac <- prcomp(t(atac))
var_prcomp <- pca_atac$sdev^2

# PCs variance
pcvar <- data.frame(var=var_prcomp/sum(var_prcomp), 
                    pc=c(1:length(var_prcomp)))

# >>> Plotting
# Plot colors
bcell_colors <- readRDS(paste0(data_wd,"/osfstorage-archive/00_plot_parameters/bcell_colors.rds"))

# Plot
data_to_plot <- data.frame(PC1=pca_atac$x[,1],PC2=pca_atac$x[,2], CellType=atac_metadata$CellType)

p <- ggplot(data_to_plot, aes(x = PC1, y = PC2, color = CellType)) +
  geom_point(size = 3, show.legend=FALSE) + 
  coord_fixed(ratio=1) + 
  scale_color_manual(values = bcell_colors) + 
  directlabels::geom_dl(aes(label = CellType), method = list("smart.grid", cex=1)) + 
  labs(title="PCA: ATAC-Seq",x = paste("PC1 (",round(pcvar[1,1]*100,0),"%)", sep=""), y = paste("PC2 (",round(pcvar[2,1]*100,0),"%)", sep=""))+
  theme_classic()

pdf(paste(format(Sys.time(), "%Y%m%d"),"_pca_atac.pdf",sep=""), width=6, height = 6)
print(p)
dev.off()
# ------------------------------------------------------------------------------------ #



# Figure 1 - Panel D -----------------------------------------------------------------
# Principal Component Analysis (PCA) of TFBS activity score (78 samples).

# >>> Required packages
library(ggplot2)
library(BuenColors)

# >>> Input data
# TFBS data
tfbs <- read.table(paste0(data_wd,"/osfstorage-archive/03_TFs_OCR_links/chromVAR_tf/TFBS_deviations_from_ocr_tf_matrix.txt"))
atac_metadata <- readRDS(paste0(data_wd,"/osfstorage-archive/01_atac_seq_data/atac_metadata.rds"))

# >>> Computation
# PCA computation
pca_tfbs <- prcomp(t(tfbs))
var_prcomp <- pca_tfbs$sdev^2

# PCs variance
pcvar <- data.frame(var=var_prcomp/sum(var_prcomp), 
                    pc=c(1:length(var_prcomp)))

# >>> Plotting
# Plot colors
bcell_colors <- readRDS(paste0(data_wd,"/osfstorage-archive/00_plot_parameters/bcell_colors.rds"))

# Plot
data_to_plot <- data.frame(PC1=pca_tfbs$x[,1],PC2=pca_tfbs$x[,2], CellType=atac_metadata$CellType)

p <- ggplot(data_to_plot, aes(x = PC1, y = PC2, color = CellType)) +
  geom_point(size = 3, show.legend=FALSE) + 
  coord_fixed(ratio=1) + 
  scale_color_manual(values = bcell_colors) + 
  directlabels::geom_dl(aes(label = CellType), method = list("smart.grid", cex=1)) + 
  labs(title="PCA: TFBS - Deviation",x = paste("PC1 (",round(pcvar[1,1]*100,0),"%)", sep=""), y = paste("PC2 (",round(pcvar[2,1]*100,0),"%)", sep=""))+
  theme_classic()

pdf(paste(format(Sys.time(), "%Y%m%d"),"_pca_tfbs.pdf",sep=""), width=6, height = 6)
print(p)
dev.off()
# ------------------------------------------------------------------------------------ #



# Figure 1 - Panel E -----------------------------------------------------------------
# Heatmap representation of TFBS activity score for key B-cell regulators. 
# Samples are shown by columns (n=78) and grouped by cell subpopulations, and TFs by rows (n=12). 

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

# Define TFs to plot
target_gene <- c("GATA2","TAL1","MEIS1","FOXO1","IRF4","IRF8","EBF1","PAX5","TCF3","IKZF1","SPIB","POU2F2")

# Plot
tfbs_plot <- Heatmap(tfbs[target_gene,],
        top_annotation = HeatmapAnnotation(Subpopulation = atac_metadata$CellType, 
                                           col = list(Subpopulation = bcell_colors),
                                           show_annotation_name = FALSE),
        col = jdb_palette("solar_extra"),
        cluster_columns = TRUE,
        column_split = atac_metadata$CellType,
        cluster_rows = FALSE,
        column_gap = unit(0.5, "mm"),
        row_split = c(rep("1",3),rep("2",3),rep("3",2),rep("4",2),rep("5",2)),
        show_column_names = FALSE,
        width = unit(4, "cm"),
        height = unit(3, "mm")*length(target_gene),
        row_names_gp = gpar(fontsize = 6),
        row_title=NULL,
        column_title = NULL,
        name="TFBS score"
)

pdf(paste(format(Sys.time(), "%Y%m%d"),"_heatmap_target_tfbs_zscore_clustered.pdf",sep=""))
  tfbs_plot
dev.off()
# ------------------------------------------------------------------------------------ #



# Figure 1 - Panel F -----------------------------------------------------------------
# GSEA over B-cell transitions
# Representation of top ten enriched GO terms derived from the transcriptional changes over B-cell transitions. 
# The normalized enrichment score is shown in a red-blue color range for each term and the –log10 of the adjusted p-value as dot size.
# Biological functions are grouped into four main terms.

# >>> Required packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

# >>> Required auxiliar functions
source("auxiliar_functions.R")

# >>> Input data
# DEA by transitions
dea <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/differential_expression_analysis/rna_dea_transitions.rds"))

# >>> Computation
# Insted of run this section you can upload the following file:
# "/osfstorage-archive/01_rna_seq_data/differential_expression_analysis/rna_dea_transitions_GSEA.rds"
# corresponding to gseaGO_bulk

set.seed(12345)

# 1 - Creating geneList (list with the sorted gene fold-change for each transitions)
comp <- c("HSC_CLP","CLP_ProB","ProB_PreB","PreB_Immature.B",
          "Immature.B_Transitional.B","Transitional.B_Naive.CD5neg",
          "Transitional.B_Naive.CD5pos")

num_comp <- length(comp)
geneList <- vector(length=num_comp, mode = "list")
names(geneList) <- comp

for(i in 1:num_comp){
  fold_change <- (dea[,paste("logFC_",comp[i],sep="")])*-1 
  names(fold_change) <- rownames(dea)
  geneList[[i]] <- fold_change[order(fold_change, decreasing = TRUE)]
}

# 2 - Run GSEA on GO terms
gseaGO_bulk <- compareCluster(geneClusters = geneList,
                              fun = "gseGO",
                              OrgDb = org.Hs.eg.db,
                              keyType = "ENSEMBL",
                              ont = "BP",
                              minGSSize = 20,
                              maxGSSize = 500,
                              pvalueCutoff = 0.01)

# >>> Plotting
gseaGO_to_plot <- pairwise_termsim(gseaGO_bulk)

pdf(paste0(format(Sys.time(), "%Y%m%d"),"_treeplot_gsea_transitions_top10_c4.pdf"))
  plot_treedot(comp_pair_term=gseaGO_to_plot, top_paths=10, clust_num=4)
dev.off()
# ------------------------------------------------------------------------------------ #



# Figure 1 - Panel G -----------------------------------------------------------------
# Heatmap representation of gene expression profile for 58 well-known genes associated with B-cell lymphopoiesis. 

# >>> Required packages
library(ComplexHeatmap)
library(BuenColors)

# >>> Input data
# Gene Expression data
rna <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/rna_norm_data.rds"))
rna_metadata <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/rna_metadata.rds"))
rna_anno <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/rna_gene_annotation.rds"))
# DEA by transitions
dea <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/differential_expression_analysis/rna_dea_transitions.rds"))


# >>> Plotting
# Plot colors
bcell_colors <- readRDS(paste0(data_wd,"/osfstorage-archive/00_plot_parameters/bcell_colors.rds"))

# Select genes to plot
key_bcell_genes <- data.frame(gene=c("CRHBP","HOPX","MLLT3","HLF","AVP",
                                     "CD34","PROM1","KIT","SMIM24","SPINK2","C1QTNF4","FLT3",
                                     "CD38","CD19","MS4A1","CD5","CD24","MME",
                                     "RAG1","RAG2","DNTT",
                                     "VPREB1","IGLL1","CD79A","CD79B","IGHM","IGHD",
                                     "HMGB1","HMGB2","HMGB3","MKI67","TOP2A","CKS2","NASP","PTTG1","STMN1","SMC4","NUSAP1",
                                     "MCL1","BCL2","BCL2L1","SOCS1","AKT1","TNFRSF13C",
                                     "CXCR4","IL7R",
                                     "LEF1","IKZF3","IRF4","IRF8","FOXO1","EBF1","TCF3","PAX5","MYB","IKZF1","SPI1","SOX4"),
                              class=factor(c(rep("Stemness",5),
                                             rep("Lymphoid Progenitors",7),
                                             rep("B-cell lineage",6),                                      
                                             rep("Ig-recombination",3),
                                             rep("pre-BCR/BCR",6),
                                             rep("Proliferation",11),
                                             rep("Survival",6),
                                             rep("Differentiation",2),
                                             rep("Transcription factors",12)), 
                                           levels=c("Stemness","Lymphoid Progenitors","B-cell lineage",
                                                    "Ig-recombination","pre-BCR/BCR",
                                                    "Differentiation","Proliferation","Survival","Transcription factors"))
)

rownames(key_bcell_genes) <- key_bcell_genes$gene

# Prepare data to plot
rna_anno2 <- rna_anno[rna_anno$hgnc_symbol%in%key_bcell_genes$gene,]
data_to_plot <- rna[rna_anno2$ensembl_gene_id,]
rownames(data_to_plot) <- rna_anno2$hgnc_symbol

key_bcell_genes_ordered <- key_bcell_genes[rownames(data_to_plot),]

# Plot
exprs_pattern <- Heatmap(data_to_plot, 
                         top_annotation = HeatmapAnnotation(Subpopulation = rna_metadata$CellType,
                                                            col = list(Subpopulation = bcell_colors),
                                                            show_annotation_name = FALSE), 
                         row_title_gp = gpar(fontsize = 7),
                         cluster_columns = FALSE,
                         column_split = rna_metadata$CellType,
                         col=jdb_palette("wolfgang_extra"),
                         row_split = key_bcell_genes_ordered$class,
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

pdf(paste0(format(Sys.time(), "%Y%m%d"),"_bcell_expression_biomarkers.pdf"), height = 10)
  exprs_pattern
dev.off()
# ------------------------------------------------------------------------------------ #



# Figure 1 - Panel H -----------------------------------------------------------------
# B-cell signature enrichment in mice subpopulations using GSVA
# Heatmap representation of the gene set activity score (GSVA) computed 
# for each human B-cell subpopulation signature (n=8; in rows) in each mouse sample (n=17; in columns). 
# The GSVA score represented in the heatmap is scaled by column to highlight the human B-cell subpopulation
# with higher enrichment in each mouse B-cell subpopulation (in dark blue).

# >>> Required packages
library(ComplexHeatmap)
library(BuenColors)
library(GSVA)

# >>> Required auxiliar functions
source("auxiliar_functions.R")

# >>> Input data
# Gene Expression data from mice (GSE109125)
rna_mice <- readRDS(paste0(data_wd,"/osfstorage-archive/06_mice/rna_norm_data_mice.rds"))
rna_metadata_mice <- readRDS(paste0(data_wd,"/osfstorage-archive/06_mice/rna_metadata_mice.rds"))
# Human B-cell subpopulations biomarkers
bcell_biomarkers <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/differential_expression_analysis/Bcell_biomarkers.rds"))

# >>> Computation:

# 1 - Mouse to human gene names correspondence
orthology <- convertMouseGeneList(rownames(rna_mice))

# 2 - Translate B-cell subpopulations biomarkers to mice gene symbols
bcell_biomarkers_mice <- bcell_biomarkers

for(i in names(bcell_biomarkers)){
  bcell_biomarkers_mice[[i]] <- unique(orthology[orthology$Gene.stable.ID%in%bcell_biomarkers[[i]],"MGI.symbol"])
}

# 3 - GSVA score of B-cell biomarkers
gs <- bcell_biomarkers_mice
names(gs) <- names(bcell_colors)
gsvaPar <- gsvaParam(rna_mice, gs)
gsva.es <- gsva(gsvaPar, verbose=FALSE)

# >>> Plotting
# Plot colors
bcell_colors <- readRDS(paste0(data_wd,"/osfstorage-archive/00_plot_parameters/bcell_colors.rds"))
bcell_colors_mice <- readRDS(paste0(data_wd,"/osfstorage-archive/00_plot_parameters/bcell_colors_mice.rds"))

# Plot
h <- Heatmap(scale(gsva.es),
             top_annotation = HeatmapAnnotation(mice = rna_metadata_mice$Bcell_subtype, 
                                                col = list(mice = bcell_colors_mice),
                                                show_annotation_name = FALSE),
             left_annotation = rowAnnotation(human = factor(names(bcell_colors), levels=names(bcell_colors)), 
                                             col = list(human = bcell_colors),
                                             show_annotation_name = FALSE),
             col = jdb_palette("brewer_purple"),
             column_order = order(rna_metadata_mice$Bcell_subtype),
             column_split = rna_metadata_mice$Bcell_subtype,
             column_gap = unit(0.5, "mm"),
             row_order = names(bcell_colors),
             row_split = factor(names(bcell_colors), levels=names(bcell_colors)),
             row_gap = unit(0.5, "mm"),
             show_column_names = FALSE,
             show_row_names = FALSE,
             column_title = NULL,
             row_title = NULL,
             width = unit(4, "cm"),
             height = unit(4, "cm"),
             row_names_gp = gpar(fontsize = 6),
             name="GSVA score (scaled)"
)

pdf(paste(format(Sys.time(), "%Y%m%d"),"_heatmap_gsva_Bcell_biomarkers_in_mice.pdf",sep=""))
  draw(h)
dev.off()
# ------------------------------------------------------------------------------------ #


