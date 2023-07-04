#***************************************************
## Paper: "Uncovering the cis-regulatory program of early human B-cell commitment and its implications
## in the pathogenesis of acute lymphoblastic leukemia"
## Authors: Planell N et al.
## Date: 2023
## Code: Figure 2
## Input data available at: 
#***************************************************




#***************************************************
## Figure 2 A
#***************************************************
## Heatmap representation of genes included in the GRN (n=7,310 genes). 

library(ggplot2)
library(BuenColors)
library(circlize)
library(ComplexHeatmap)

# Input data:
# Normalized RNA: "rna_norm_data.txt" 
  # available at OSF files: 01_rna_seq_data/rna_norm_data.txt
# Metadata RNA: "rna_metadata.txt" 
  # available at OSF files: 01_rna_seq_data/rna_metadata.txt
# Gene annotation: "rna_gene_annotation.txt"
  # available at OSF files: 01_rna_seq_data/rna_gene_annotation.txt
# List TF-OCR-gene interactions: "tf_ocr_gene_interactions_curated.txt"
  # available at OSF files: 04_gene_regulatory_networks/tf_ocr_gene_interactions_curated.txt
# DEGs clustering:  "rna_deg_clustering.txt"
  # available at OSF files: 01_rna_seq_data/differential_expression_analysis/rna_deg_clustering.txt
# B cell subpopulation colors: "bcell_colors.rds"
# RNA expression colors: "rna_expression_colors.rds"

rna_norm <- read.table("rna_norm_data.txt",header=TRUE, dec=".", sep="\t")
rna_metadata <- read.table("rna_metadata.txt",header=TRUE, dec=".", sep="\t")
rna_annotation <- read.table("rna_gene_annotation.txt",header=TRUE, dec=".", sep="\t")

bcell_colors <- readRDS("bcell_colors.rds")

tf_ocr_gene <- read.table("tf_ocr_gene_interactions_curated.txt",header=TRUE, dec=".", sep="\t")
deg_clust <- read.table("rna_deg_clustering.txt",header=TRUE, dec=".", sep="\t")

data_to_plot <- t(scale(t(rna_norm[unique(tf_ocr_gene$gene_ensembl),]),center=TRUE, scale=TRUE))
data_to_plot2 <- deg_clust[unique(tf_ocr_gene$gene_ensembl),]
row_order <- order(data_to_plot2$pattern_id)

data_to_plot <- data_to_plot[row_order,]
data_to_plot2 <- data_to_plot2[row_order,]

ha = HeatmapAnnotation(stage = rna_metadata$CellType,
                       col = list(stage = bcell_colors))

ha2 = HeatmapAnnotation(stage = labels(bcell_colors),
                        col = list(stage = bcell_colors))

target_genes <- c("PAX5","HOXB6","BCL2","NOTCH1","CDH2","MEIS1","CCND1",
                  "EBF1","TAL1","GATA1","TCF3","CD19","MS4A1","MME","CD79A","CD79B",
                  "CD38","CD34","CD5","IGHM","DNTT","HMGB2","EZH2","HMGB1","MYBL2","BIRC5",
                  "CD83","BTG1")

sel_pos <- vector()
target_genes2 <- vector()
for(i in 1:length(target_genes)){
  ss <- which(rownames(data_to_plot)==rna_annotation$ensembl_gene_id[rna_annotation$hgnc_symbol==target_genes[i]])
  if(length(ss)>0){
    sel_pos <- c(sel_pos,ss)
    target_genes2 <- c(target_genes2,target_genes[i])
  }
}

rha_text <- rowAnnotation(genes=anno_mark(at = sel_pos, labels = target_genes2))

my_colors_rna <- readRDS("rna_expression_colors.rds")

exprs_pattern <- Heatmap(data_to_plot, name = "Expression patterns",
                         top_annotation = ha, column_order=order(rna_metadata$CellType), 
                         col=colorRamp2(seq(-3,3, length.out=length(my_colors_rna)),my_colors_rna),
                         show_column_names = FALSE, show_row_names = FALSE, row_names_gp = gpar(fontsize = 3), 
                         cluster_rows = FALSE, row_split=data_to_plot2$pattern_id,
                         cluster_columns = FALSE, row_gap = unit(0.5, "mm"), use_raster=TRUE, width = unit(4, "cm"))

split <- data.frame(c(1:length(labels(bcell_colors))))
rownames(split) <- labels(bcell_colors)

binary_pattern <- Heatmap(data_to_plot2[,5:12], name = "Peaks vs Cell Type",
                          col = colorRamp2(c(0, 1), c("beige","darkgoldenrod")),
                          top_annotation = ha2,
                          right_annotation = rha_text, 
                          show_column_names = TRUE, show_row_names = FALSE, row_names_gp = gpar(fontsize = 3), 
                          cluster_rows = FALSE, row_split=data_to_plot2$pattern_id,
                          cluster_columns = FALSE, row_gap = unit(0.5, "mm"), column_gap = unit(1, "mm"), use_raster=TRUE, 
                          border="grey", column_split=split, width = unit(3, "cm"))

exprs_pattern + binary_pattern

#***************************************************




#***************************************************
## Figure 2 B
#***************************************************
## Heatmap representation of OCRs included in the GRN (n=16,074 OCRs). 

library(ggplot2)
library(BuenColors)
library(circlize)
library(ComplexHeatmap)

# Input data:
# Normalized ATAC: "atac_norm_data.txt" 
  # available at OSF files: 01_atac_seq_data/atac_norm_data.txt
# Metadata ATAC: "atac_metadata.txt"
  # available at OSF files: 01_atac_seq_data/atac_metadata.txt
# OCR annotation: "atac_OCR_annotation.txt"
  # available at OSF files: 01_atac_seq_data/atac_OCR_annotation.txt
# List TF-OCR-gene interactions: "tf_ocr_gene_interactions_curated.txt"
  # available at OSF files: 04_gene_regulatory_networks/tf_ocr_gene_interactions_curated.txt
# DARs clustering:  "atac_dar_clustering.txt"
  # available at OSF files: 01_atac_seq_data/differential_binding_analysis/atac_dar_clustering.txt
# B cell subpopulation colors: "bcell_colors.rds"
# ATAC expression colors: "atac_colors.rds"

atac_norm <- read.table("atac_norm_data.txt",header=TRUE, dec=".", sep="\t") 
atac_metadata <- read.table("atac_metadata.txt",header=TRUE, dec=".", sep="\t") 
atac_annotation <- read.table("atac_OCR_annotation.txt",header=TRUE, dec=".", sep="\t", quote="")

bcell_colors <- readRDS("bcell_colors.rds")

tf_ocr_gene <- read.table("tf_ocr_gene_interactions_curated.txt",header=TRUE, dec=".", sep="\t")
dar_clust <- read.table("atac_dar_clustering.txt",header=TRUE, dec=".", sep="\t")

data_to_plot <- t(scale(t(atac_norm[unique(tf_ocr_gene$ocr),]),center=TRUE, scale=TRUE))
data_to_plot2 <- dar_clust[unique(tf_ocr_gene$ocr),]
row_order <- order(data_to_plot2$pattern_id)

data_to_plot <- data_to_plot[row_order,]
data_to_plot2 <- data_to_plot2[row_order,]
atac_annotation2 <- atac_annotation[rownames(data_to_plot),]

ha = HeatmapAnnotation(stage = atac_metadata$CellType,
                       col = list(stage = bcell_colors))

ha2 = HeatmapAnnotation(stage = labels(bcell_colors),
                        col = list(stage = bcell_colors))

peak_color <- c("#F6313E","#46A040","#490C65")
names(peak_color) <- c(unique(atac_annotation2$annotation_simplified))	   

rha = rowAnnotation(peak = atac_annotation2$annotation_simplified,
                    col = list(peak = peak_color))

target_genes <- c("PAX5","HOXB6","BCL2","NOTCH1","CDH2","MEIS1","CCND1",
                  "EBF1","TAL1","GATA1","TCF3","CD19","MS4A1","MME","CD79A","CD79B",
                  "CD38","CD34","CD5","IGHM","DNTT","HMGB2","EZH2","HMGB1","MYBL2","BIRC5",
                  "CD83","BTG1")

sel_pos <- vector()
target_genes2 <- vector()
for(i in 1:length(target_genes)){
  ss <- which(rownames(data_to_plot)%in%rownames(atac_annotation2)[which(atac_annotation2$annotation_simplified=="Promoter" & atac_annotation2$SYMBOL==target_genes[i])])
  if(length(ss)>0){
    sel_pos <- c(sel_pos,ss)
    target_genes2 <- c(target_genes2,paste0(target_genes[i],"_",rownames(data_to_plot)[ss]))
  }
}

rha_text <- rowAnnotation(genes=anno_mark(at = sel_pos, labels = target_genes2))

my_colors_atac <- readRDS("atac_color.rds")

exprs_pattern <- Heatmap(data_to_plot, name = "Expression patterns",
                         top_annotation = ha, column_order=order(atac_metadata$CellType),
                         right_annotation = rha,		
                         col=colorRamp2(seq(-3,3, length.out=length(my_colors_atac)),my_colors_atac),
                         show_column_names = FALSE, show_row_names = FALSE, row_names_gp = gpar(fontsize = 3), 
                         cluster_rows = FALSE, row_split=data_to_plot2$pattern_id,
                         cluster_columns = FALSE, row_gap = unit(0.5, "mm"), use_raster=TRUE, width = unit(4, "cm"))

split <- data.frame(c(1:length(labels(bcell_colors))))
rownames(split) <- labels(bcell_colors)

binary_pattern <- Heatmap(data_to_plot2[,5:12], name = "Peaks vs Cell Type",
                          col = colorRamp2(c(0, 1), c("beige","darkgoldenrod")),
                          top_annotation = ha2,
                          right_annotation = rha_text, 
                          show_column_names = TRUE, show_row_names = FALSE, row_names_gp = gpar(fontsize = 3), 
                          cluster_rows = FALSE, row_split=data_to_plot2$pattern_id,
                          cluster_columns = FALSE, row_gap = unit(0.5, "mm"), column_gap = unit(1, "mm"), use_raster=TRUE, 
                          border="grey", column_split=split, width = unit(3, "cm"))

exprs_pattern + binary_pattern

#***************************************************




#***************************************************
## Figure 2 C
#***************************************************
## TFBS activity score shown as Z-score for all TFs within the GRN (n=169 TFs).  

library(ggplot2)
library(BuenColors)
library(circlize)
library(ComplexHeatmap)

# Input data:
# TFBS activity score: "TFBS_scores_from_ocr_tf_matrix_curated.txt"
  # available at OSF files: 03_TFs_OCR_links/chromVAR_tf_cCREs/TFBS_scores_from_ocr_tf_matrix_curated.txt
# Metadata ATAC: "atac_metadata.txt"
  # available at OSF files: 01_atac_seq_data/atac_metadata.txt
# B cell subpopulation colors: "bcell_colors.rds"
# TFBS activity score colors: "tfbs_color.rds"

tfbs_zscore <- read.table("TFBS_scores_from_ocr_tf_matrix_curated.txt",header=TRUE, dec=".", sep="\t") 
atac_metadata <- read.table("atac_metadata.txt",header=TRUE, dec=".", sep="\t") 

bcell_colors <- readRDS("bcell_colors.rds")
atac_metadata$CellType <- factor(atac_metadata$CellType, levels=names(bcell_colors))

data_to_plot <- tfbs_zscore

ha = HeatmapAnnotation(stage = atac_metadata$CellType,
                       col = list(stage = bcell_colors))

target_TF <- c("SPIB","IRF4","SPI1","BCL11A","JUND","FLI1","PAX5","FOXO4","FOXO1",
               "GFI1","IKZF1","TAL1","GATA1","MYB","MLXIP","RUNX3","ELK3","ETV6",
               "ETS2","ERG","E2F3","RUNX1","EBF1","TCF3","TCF4","LEF1")

target_TF_position <- vector()
for(i in 1:length(target_TF)){
  target_TF_position <- c(target_TF_position,grep(target_TF[i],rownames(data_to_plot)))
}

rha = rowAnnotation(foo = anno_mark(at = target_TF_position, 
                                    labels = rownames(data_to_plot)[target_TF_position],
                                    labels_gp = gpar(fontsize = 6)))

my_colors_tfbs <- readRDS("tfbs_color.rds")

exprs_pattern <- Heatmap(data_to_plot, name = "TFBS patterns",
                         top_annotation = ha, column_order=order(atac_metadata$CellType),
                         right_annotation = rha,
                         col=colorRamp2(seq(-1*min(abs(min(data_to_plot)),max(data_to_plot)), min(abs(min(data_to_plot)),max(data_to_plot)), length.out=9),my_colors_tfbs),
                         show_column_names = FALSE, show_row_names = FALSE, 
                         cluster_rows = TRUE)

exprs_pattern

#***************************************************




#***************************************************
## Figure 2 D
#***************************************************
## Multidimensional Scaling representation from the Pearson correlation of cell subpopulation specific GRNs.

###-----------------------------------------------------------------------------------###




#***************************************************
## Figure 2 E
#***************************************************
library(BuenColors)
library(ComplexHeatmap)
library(igraph)

# Input data:
# igraph object for whole GRNs: bcell_grn.rds
  #available at OSF files: 04_gene_regulatory_networks/bcell_grn.rds
# igraph object for cell type specific GRNs: bcell_grn_by_cell_type.rds
  #available at OSF files: 04_gene_regulatory_networks/bcell_grn_by_cell_type.rds
# B cell subpopulation colors: "bcell_colors.rds"
# List TF-OCR-gene interactions: "tf_ocr_gene_interactions_curated.txt"
  # available at OSF files: 04_gene_regulatory_networks/tf_ocr_gene_interactions_curated.txt
# Gene annotation: "rna_gene_annotation.txt"
  # available at OSF files: 01_rna_seq_data/rna_gene_annotation.txt

bcell_colors <- readRDS("bcell_colors.rds")

tf_ocr_gene <- read.table("tf_ocr_gene_interactions_curated.txt",header=TRUE, dec=".", sep="\t")
gene_anno <- read.table("rna_gene_annotation.txt",header=TRUE, dec=".", sep="\t")

# Whole GRN
grn <- readRDS("bcell_grn.rds")
  
# Centrality measures
V(grn)$degree <- igraph::degree(grn)        # Degree centrality
V(grn)$betweenness <- betweenness(grn)      # Vertex betweenness centrality

# GRN specific per cell type
list_grn <- readRDS("bcell_grn_by_cell_type.rds")

for(i in 1:length(list_grn)){
  # Centrality measures
  V(list_grn[[i]])$degree <- igraph::degree(list_grn[[i]])      # Degree centrality
  V(list_grn[[i]])$betweenness <- betweenness(list_grn[[i]])    # Vertex betweenness centrality
}
names(list_grn) <- cell_type


# Data frame with centrality measures

centrality <- data.frame(row.names   = V(grn)$name,						
                         degree      = V(grn)$degree,
                         betweenness = V(grn)$betweenness)

for(i in 1:length(list_grn)){
  centrality2 <- data.frame(row.names   = V(list_grn[[i]])$name,						
                            degree      = V(list_grn[[i]])$degree,
                            betweenness = V(list_grn[[i]])$betweenness)
  colnames(centrality2) <- paste0(colnames(centrality2),"_",names(list_grn)[i])
  centrality <- merge(centrality,centrality2, by.x=0, by.y=0, all=TRUE)
  rownames(centrality) <- centrality$Row.names
  centrality$Row.names <- NULL
}

centrality_tf <- centrality[unique(tf_ocr_gene$tf_ensembl),]
rownames(centrality_tf) <- gene_anno[rownames(centrality_tf),"hgnc_symbol"]


# Heatmap representation
ha = HeatmapAnnotation(stage = cell_type,
                       col = list(stage = bcell_colors))

# Degree representation of whole GRN
grn_properties_degree <- Heatmap(centrality_tf[,"degree"], 
                                 row_names_gp = gpar(fontsize = 2), width = unit(1, "cm"),
                                 col=jdb_palette("flame_light"),name="degree",
                                 na_col = "white")

# Heatmap used to define the clustering 
mm <- centrality_tf[,grep("degree_",colnames(centrality_tf),value=TRUE)]
mm[is.na(mm)] <- 0
grn_properties_degree_cell_type0 <- Heatmap(mm, 
                                            row_names_gp = gpar(fontsize = 2), width = unit(4, "cm"),
                                            top_annotation=ha,
                                            col=jdb_palette("flame_light"),name="degree0")

# Degree representation of cell type specific GRN
grn_properties_degree_cell_type <- Heatmap(centrality_tf[,grep("degree_",colnames(centrality_tf),value=TRUE)], 
                                           row_names_gp = gpar(fontsize = 2), width = unit(4, "cm"),
                                           col=jdb_palette("flame_light"),name="degree2",
                                           top_annotation=ha,
                                           na_col = "white", cluster_rows=FALSE, cluster_columns=FALSE)


# Betweeness representation of whole GRN
grn_properties_betweenness <- Heatmap(centrality_tf[,c("betweenness")], 
                                      row_names_gp = gpar(fontsize = 2), width = unit(1, "cm"),
                                      col=jdb_palette("flame_artic"), name="betweenness",
                                      na_col = "white")


# Betweeness representation of cell type specific GRN
grn_properties_betweenness_cell_type <- Heatmap(centrality_tf[,grep("betweenness_",colnames(centrality_tf),value=TRUE)], 
                                                row_names_gp = gpar(fontsize = 2), width = unit(4, "cm"),
                                                col=jdb_palette("flame_artic"), name="betweenness2",
                                                top_annotation=ha,
                                                na_col = "white", cluster_rows=FALSE, cluster_columns=FALSE)

# Combined plot
grn_properties_degree_cell_type0+grn_properties_degree+grn_properties_betweenness+grn_properties_degree_cell_type+grn_properties_betweenness_cell_type

###-----------------------------------------------------------------------------------###




#***************************************************
## Figure 2 F
#***************************************************

# After generate Figure 2 E 

# Identify the top TF per cell type based on betweenness
sel <- unique(c(order(centrality_tf$betweenness_HSC,decreasing=TRUE)[1:10],
                order(centrality_tf$betweenness_CLP,decreasing=TRUE)[1:10],
                order(centrality_tf$betweenness_proB,decreasing=TRUE)[1:10],
                order(centrality_tf$betweenness_preB,decreasing=TRUE)[1:10],
                order(centrality_tf$betweenness_ImmatureB,decreasing=TRUE)[1:10],
                order(centrality_tf$betweenness_Transitional_B,decreasing=TRUE)[1:10],
                order(centrality_tf$betweenness_Naive_CD5pos,decreasing=TRUE)[1:10],
                order(centrality_tf$betweenness_Naive_CD5neg,decreasing=TRUE)[1:10]))

centrality_tf <- centrality_tf[sel,]

# Degree representation of whole GRN
grn_properties_degree <- Heatmap(centrality_tf[,"degree"], 
                                 row_names_gp = gpar(fontsize = 2), width = unit(1, "cm"),
                                 col=jdb_palette("flame_light"),name="degree",
                                 na_col = "white")

# Heatmap used to define the clustering 
mm <- centrality_tf[,grep("degree_",colnames(centrality_tf),value=TRUE)]
mm[is.na(mm)] <- 0
grn_properties_degree_cell_type0 <- Heatmap(mm, 
                                            row_names_gp = gpar(fontsize = 2), width = unit(4, "cm"),
                                            top_annotation=ha,
                                            col=jdb_palette("flame_light"),name="degree0")

# Degree representation of cell type specific GRN
grn_properties_degree_cell_type <- Heatmap(centrality_tf[,grep("degree_",colnames(centrality_tf),value=TRUE)], 
                                           row_names_gp = gpar(fontsize = 2), width = unit(4, "cm"),
                                           col=jdb_palette("flame_light"),name="degree2",
                                           top_annotation=ha,
                                           na_col = "white", cluster_rows=FALSE, cluster_columns=FALSE)


# Betweeness representation of whole GRN
grn_properties_betweenness <- Heatmap(centrality_tf[,c("betweenness")], 
                                      row_names_gp = gpar(fontsize = 2), width = unit(1, "cm"),
                                      col=jdb_palette("flame_artic"), name="betweenness",
                                      na_col = "white")


# Betweeness representation of cell type specific GRN
grn_properties_betweenness_cell_type <- Heatmap(centrality_tf[,grep("betweenness_",colnames(centrality_tf),value=TRUE)], 
                                                row_names_gp = gpar(fontsize = 6), width = unit(4, "cm"),
                                                col=jdb_palette("flame_artic"), name="betweenness2",
                                                top_annotation=ha,
                                                na_col = "white", cluster_rows=FALSE, cluster_columns=FALSE)

# Combined plot
grn_properties_degree_cell_type0+grn_properties_degree+grn_properties_betweenness+grn_properties_degree_cell_type+grn_properties_betweenness_cell_type

###-----------------------------------------------------------------------------------###




#***************************************************
## Figure 2 G
#***************************************************
library(ggplot2)
library(ggridges)

# ELK3 exploration

# Input data:
# Normalized RNA: "rna_norm_data.txt" 
  # available at OSF files: 01_rna_seq_data/rna_norm_data.txt
# Metadata RNA: "rna_metadata.txt" 
  # available at OSF files: 01_rna_seq_data/rna_metadata.txt
# Gene annotation: "rna_gene_annotation.txt"
  # available at OSF files: 01_rna_seq_data/rna_gene_annotation.txt
# RNA expression colors: "rna_expression_colors.rds"


# RNA - Ridgeplot ELK3

rna_norm <- read.table("rna_norm_data.txt",header=TRUE, dec=".", sep="\t")
rna_metadata <- read.table("rna_metadata.txt",header=TRUE, dec=".", sep="\t")
rna_annotation <- read.table("rna_gene_annotation.txt",header=TRUE, dec=".", sep="\t")
rna_color <- readRDS("rna_expression_colors.rds")

group <- factor(as.character(rna_metadata$CellType),
                levels=c("Naive_CD5pos","Naive_CD5neg","Transitional_B","ImmatureB","preB","proB","CLP","HSC"))

data_to_plot <- data.frame(cell_type=group, score=as.numeric(rna_norm[rna_annotation$ensembl_gene_id[rna_annotation$hgnc_symbol=="ELK3"],]))
ggplot(data_to_plot, aes(x = score, y = cell_type, fill= ..x..)) +
        geom_density_ridges_gradient(scale=2) +
        scale_fill_gradientn(colors = rna_color) +
        labs(title = "ELK3 - RNA expression (Z-score)") +
        theme_ridges()


# TFBS - Ridgeplot ELK3

tfbs_score <- read.table("TFBS_scores_from_ocr_tf_matrix_curated.txt",header=TRUE, dec=".", sep="\t")
atac_metadata <- read.table("atac_metadata.txt",header=TRUE, dec=".", sep="\t") 
tfbs_color <- readRDS("tfbs_color.rds")

group <- factor(as.character(atac_metadata$CellType),
                levels=c("Naive_CD5pos","Naive_CD5neg","Transitional_B","ImmatureB","preB","proB","CLP","HSC"))


data_to_plot <- data.frame(cell_type=group, score=as.numeric(tfbs_score["ELK3",]))
ggplot(data_to_plot, aes(x = score, y = cell_type, fill= ..x..)) +
  geom_density_ridges_gradient(scale=2) +
  scale_fill_gradientn(colors = tfbs_color) +
  labs(title = "ELK3 - TFBS (Z-score)") +
  theme_ridges()

###-----------------------------------------------------------------------------------###
