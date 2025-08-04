#***************************************************
## Paper: "Uncovering the regulatory landscape of early human B-cell lymphopoiesis and its implications in the pathogenesis of B-ALL"
## Authors: Planell N et al.
## Date: 2025
## Code: Figure 4
## Input data available at: https://osf.io/gswpy/?view_only=39cfd50d5a7e44e5a2d141869dcb1315
#***************************************************

# Define the directory were input data from OSF has been downloaded
data_wd <- setwd("../GitHub_code/")

# ****** Figure 5 ****** #

# Single-cell analysis -------------------------------------------------

# Figure 4 - Panel A, B -----------------------------------------------------------------
# UMAP representation of single-cell datasets with cell-type annotation

# >>> Required packages
library(Seurat)

# >>> Input data
# Public single-cell dataset - seurat object
scPublic <- readRDS("scPublic.rds")
# In-house single-cell dataset - seurat object
scInHouse <- readRDS("scInHouse.rds")

# >>> Ploting
# Plot colors
bcell_colors_sc <- readRDS(paste0(data_wd,"/osfstorage-archive/00_plot_parameters/bcell_colors_sc.rds"))

# Plot UMAP
# Public single-cell dataset - seurat object
DimPlot(scPublic, reduction = "umap.integrated.harmony_subset", label=TRUE, group.by="cell_type_reannotation", cols = bcell_colors_sc)
  
# In-house single-cell dataset - seurat object
DimPlot(scInHouse, reduction = "umap", label=TRUE, group.by="cell_type_reannotation", cols = bcell_colors_sc)
###-----------------------------------------------------------------------------------###


# Figure 4 - Panel C, D -----------------------------------------------------------------
# Heatmap with regulons score based on AUCell

# >>> Required packages
library(Seurat)
library(ComplexHeatmap)
library(scales)
library(circlize)
library(BuenColors)

# >>> Input data
# Public single-cell dataset - seurat object
scPublic <- readRDS("scPublic.rds")
# In-house single-cell dataset - seurat object
scInHouse <- readRDS("scInHouse.rds")
# Regulons
regulons_group <- readRDS(paste0(data_wd,"/osfstorage-archive/04_gene_regulatory_networks/regulons/TF_regulons_clusters.rds"))
regulons_group_plot <- c(rep("1",length(regulons_group[[1]])),rep("2",length(regulons_group[[2]])),rep("3",length(regulons_group[[3]])),rep("4",length(regulons_group[[4]])),rep("5",length(regulons_group[[5]])))
names(regulons_group_plot) <- c(regulons_group[[1]],regulons_group[[2]],regulons_group[[3]],regulons_group[[4]],regulons_group[[5]])

# >>> Ploting
# Plot colors
regulons_color <- c("1"="#92c5de","2"="#2166ac","3"="#91cf60","4"="#dd1c77","5"="black")

cols <- jdb_palette("brewer_spectra")
col_fun <- colorRamp2(seq(-2, 2, length.out = length(cols)), cols)


# Plot  
regulons_score_to_plot <- grep("AUCell",names(scInHouse@meta.data),value=TRUE)

## Public
data_to_plot <- scPublic@meta.data[,regulons_score_to_plot]
colnames(data_to_plot) <- gsub("AUCell_","",regulons_score_to_plot)
cell_group <- scPublic@meta.data$cell_type_reannotation

h_scaled <- Heatmap(t(scale(data_to_plot)),
                    top_annotation = HeatmapAnnotation("Cell" = cell_group,
                                                       col = list("Cell"= bcell_colors_sc)),
                    left_annotation = rowAnnotation("TF_subgroup" = regulons_group_plot[colnames(data_to_plot)],
                                                    col=list("TF_subgroup" = regulons_color)),
                    col = col_fun,
                    cluster_columns = FALSE,
                    cluster_rows = FALSE,
                    column_split = cell_group,
                    row_split = factor(regulons_group_plot[colnames(data_to_plot)],levels=c("1","2","3","4","5")),
                    column_gap = unit(0.5, "mm"),
                    show_column_names = FALSE,
                    show_row_names = FALSE,
                    row_title = NULL,
                    column_title = NULL,
                    row_names_gp = gpar(fontsize = 6),
                    name="GSVA score (scaled)"
)

pdf(paste0(format(Sys.time(),"%Y%m%d"),"_heatmap_AUCell_regulons_scaled_public.pdf"))
  print(h_scaled)
dev.off()


## In-house
data_to_plot <- scInHouse@meta.data[,regulons_score_to_plot]
colnames(data_to_plot) <- gsub("AUCell_","",regulons_score_to_plot)
cell_group <- scInHouse@meta.data$cell_type_reannotation

h_scaled <- Heatmap(t(scale(data_to_plot)),
                    top_annotation = HeatmapAnnotation("Cell" = cell_group,
                                                       col = list("Cell"= bcell_colors_sc)),
                    left_annotation = rowAnnotation("TF_subgroup" = regulons_group_plot[colnames(data_to_plot)],
                                                    col=list("TF_subgroup" = regulons_color)),
                    col = col_fun,
                    cluster_columns = FALSE,
                    cluster_rows = FALSE,
                    column_split = cell_group,
                    row_split = factor(regulons_group_plot[colnames(data_to_plot)],levels=c("1","2","3","4","5")),
                    column_gap = unit(0.5, "mm"),
                    show_column_names = FALSE,
                    show_row_names = FALSE,
                    row_title = NULL,
                    column_title = NULL,
                    row_names_gp = gpar(fontsize = 6),
                    name="GSVA score (scaled)"
)

pdf(paste0(format(Sys.time(),"%Y%m%d"),"_heatmap_AUCell_regulons_scaled_inhouse.pdf"))
  print(h_scaled)
dev.off()
###-----------------------------------------------------------------------------------###



# Figure 4 - Panel E, F -----------------------------------------------------------------
# Dotplot of target TF regulons gene expression

# >>> Required packages
library(Seurat)

# >>> Input data
# Public single-cell dataset - seurat object
scPublic <- readRDS("scPublic.rds")
# In-house single-cell dataset - seurat object
scInHouse <- readRDS("scInHouse.rds")

# >>> Ploting
# Plot  
features <- rev(c("PPARA","ELK3","NR3C1","MYBL2","FOXM1","E2F7"))

## Public
DotPlot(scPublic, features=features, group.by="cell_type_reannotation") +
  theme(axis.text.x = element_text(angle = 90)) +
  coord_flip()

## In-house
DotPlot(scInHouse, features=features, group.by="cell_type_reannotation") +
  theme(axis.text.x = element_text(angle = 90)) +
  coord_flip()
###-----------------------------------------------------------------------------------###



# Figure 4 - Panel G-I -----------------------------------------------------------------
# Boxplot of target regulons activity (AUCell score)

# >>> Required packages
library(Seurat)
library(ggplot2)
library(reshape2)

# >>> Input data
# Public single-cell dataset - seurat object
scPublic <- readRDS("scPublic.rds")
# In-house single-cell dataset - seurat object
scInHouse <- readRDS("scInHouse.rds")

# >>> Ploting
# Plot colors
bcell_colors_sc <- readRDS(paste0(data_wd,"/osfstorage-archive/00_plot_parameters/bcell_colors_sc.rds"))

# Plot  
features <- c("ELK3","NR3C1","PPARA","MYBL2","FOXM1","E2F7")
tt <- paste0("AUCell_",features)

## Public
dd <- melt(scPublic@meta.data[,tt])
data_to_plot <- data.frame(value=dd$value,
                           gene=factor(gsub("AUCell_","",dd$variable),levels=c("MYBL2","FOXM1","E2F7","ELK3","NR3C1","PPARA")),
                           Cell=rep(scPublic@meta.data$cell_type_reannotation,length(tt)))

p <- ggplot(data_to_plot, aes(x=Cell, y=value, fill=Cell))+
  geom_boxplot(outlier.size = 1) +
  theme_bw() +
  scale_fill_manual(values=bcell_colors_sc) +
  facet_wrap(~gene, scales = "free", ncol=4)+
  ylab("AUCell") +
  xlab("") +
  theme(axis.text.x = element_text(angle=45, hjust=1)) 

pdf(paste0(format(Sys.Date(),"%Y%m%d"),"_AUCell_scpublic_boxplot_targetRegulons.pdf"), height=8, width=14)
print(p)
dev.off()


## In-house
dd <- melt(scInHouse@meta.data[,tt])
data_to_plot <- data.frame(value=dd$value,
                           gene=factor(gsub("AUCell_","",dd$variable),levels=c("MYBL2","FOXM1","E2F7","ELK3","NR3C1","PPARA")),
                           Cell=rep(scInHouse@meta.data$cell_type_reannotation,length(tt)))

p <- ggplot(data_to_plot, aes(x=Cell, y=value, fill=Cell))+
  geom_boxplot(outlier.size = 1) +
  theme_bw() +
  scale_fill_manual(values=bcell_colors_sc) +
  facet_wrap(~gene, scales = "free", ncol=4)+
  ylab("AUCell") +
  xlab("") +
  theme(axis.text.x = element_text(angle=45, hjust=1)) 

pdf(paste0(format(Sys.Date(),"%Y%m%d"),"_AUCell_scinhouse_boxplot_targetRegulons.pdf"), height=8, width=14)
print(p)
dev.off()
###-----------------------------------------------------------------------------------###
