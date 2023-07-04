#***************************************************
## Paper: "Uncovering the cis-regulatory program of early human B-cell commitment and its implications
## in the pathogenesis of acute lymphoblastic leukemia"
## Authors: Planell N et al.
## Date: 2023
## Code: Figure 3
## Input data available at: 
#***************************************************




#***************************************************
## Figure 3 A
#***************************************************
## GRN representation in Cytoscape
#***************************************************




#***************************************************
## Figure 3 B
#***************************************************
## GRN representation in Cytoscape
#***************************************************




#***************************************************
## Figure 3 C
#***************************************************
## EBF1 upstream regulation

library(ggplot2)
library(BuenColors)
library(circlize)
library(ComplexHeatmap)

# Input data:
# List TF-OCR-gene interactions: "tf_ocr_gene_interactions_curated.txt"
  # available at OSF files: 04_gene_regulatory_networks/tf_ocr_gene_interactions_curated.txt
# Normalized RNA: "rna_norm_data.txt" 
  # available at OSF files: 01_rna_seq_data/rna_norm_data.txt
# Metadata RNA: "rna_metadata.txt" 
  # available at OSF files: 01_rna_seq_data/rna_metadata.txt
# Gene annotation: "rna_gene_annotation.txt"
  # available at OSF files: 01_rna_seq_data/rna_gene_annotation.txt
# B cell subpopulation colors: "bcell_colors.rds"
# RNA expression colors: "rna_expression_colors.rds"

tf_ocr_gene <- read.table("tf_ocr_gene_interactions_curated.txt",header=TRUE, dec=",", sep="\t")
rna_norm <- read.table("rna_norm_data.txt",header=TRUE, dec=".", sep="\t")
rna_metadata <- read.table("rna_metadata.txt",header=TRUE, dec=".", sep="\t")
rna_annotation <- read.table("rna_gene_annotation.txt",header=TRUE, dec=".", sep="\t")

bcell_color <- readRDS("bcell_colors.rds")

# Select the upstream regulators of EBF1
tf_regulators <- unique(tf_ocr_gene[tf_ocr_gene$gene_symbol=="EBF1",c("tf_ensembl","gene_tf_rho")])

rna_metadata$CellType <- factor(rna_metadata$CellType, levels=names(bcell_colors))


# Heatmap

#Norma data
data_to_plot <- rna_norm[tf_regulators$tf_ensembl,]
rownames(data_to_plot) <- rna_annotation[rownames(data_to_plot),"hgnc_symbol"]

# RNA color 
palette <- colorRampPalette(colors=c("#E5CCFF", "#4C0099","#FFFF00"))
cols <- palette(100)
my_rna_col <- colorRamp2(seq(-10,10, length.out=length(cols)),cols)

#rho color for annotation
col_fun1 = colorRamp2(c(-0.9,-0.8,-0.7,-0.6, 0, 0.6, 0.7, 0.8, 0.9), jdb_palette("solar_flare",9,"discrete"))


target_genes <- c("PAX5", "TCF3", "FOXO1", "SPIB", "RUNX1", "CEBPA", "LEF1", "MEF2A", "IRF4", "TCF4", "POU4F1", "FOXO6", "MLXIP")
sel_pos <- vector("numeric",length=length(target_genes))
for(i in 1:length(target_genes)){
  sel_pos[i] <- which(rownames(data_to_plot)==target_genes[i])
}

ha = HeatmapAnnotation(stage = rna_metadata$CellType, 
                       col = list(stage = bcell_color),
                       rna = anno_barplot(as.numeric(rna_norm[rna_annotation$ensembl_gene_id[rna_annotation$hgnc_symbol=="EBF1"],])))


rha = rowAnnotation(rho=tf_regulators$gene_tf_rho,
                    col = list(rho = col_fun1),
                    genes=anno_mark(at = sel_pos, labels = target_genes))



Heatmap(data_to_plot, top_annotation = ha, right_annotation=rha, 
                  show_row_names = TRUE, show_column_names =FALSE, 
                  column_names_gp = gpar(fontsize = 5), row_names_gp = gpar(fontsize = 4),
                  col=my_rna_col, 
                  heatmap_legend_param = list(title = "Norm RNA"),
                  column_order = order(rna_metadata$CellType),
                  width = unit(5, "cm"), height = unit(10, "cm"),
                  row_names_side = "left")

#***************************************************




#***************************************************
## Figure 3 D
#***************************************************
## ELK3 upstream regulation

library(ggplot2)
library(BuenColors)
library(circlize)
library(ComplexHeatmap)

# Input data:
# List TF-OCR-gene interactions: "tf_ocr_gene_interactions_curated.txt"
  # available at OSF files: 04_gene_regulatory_networks/tf_ocr_gene_interactions_curated.txt
# Normalized RNA: "rna_norm_data.txt" 
  # available at OSF files: 01_rna_seq_data/rna_norm_data.txt
# Metadata RNA: "rna_metadata.txt" 
  # available at OSF files: 01_rna_seq_data/rna_metadata.txt
# Gene annotation: "rna_gene_annotation.txt"
  # available at OSF files: 01_rna_seq_data/rna_gene_annotation.txt
# B cell subpopulation colors: "bcell_colors.rds"
# RNA expression colors: "rna_expression_colors.rds"

tf_ocr_gene <- read.table("tf_ocr_gene_interactions_curated.txt",header=TRUE, dec=",", sep="\t")
rna_norm <- read.table("rna_norm_data.txt",header=TRUE, dec=".", sep="\t")
rna_metadata <- read.table("rna_metadata.txt",header=TRUE, dec=".", sep="\t")
rna_annotation <- read.table("rna_gene_annotation.txt",header=TRUE, dec=".", sep="\t")

bcell_color <- readRDS("bcell_colors.rds")

# Select the upstream regulators of ELK3
tf_regulators <- unique(tf_ocr_gene[tf_ocr_gene$gene_symbol=="ELK3",c("tf_ensembl","gene_tf_rho")])

rna_metadata$CellType <- factor(rna_metadata$CellType, levels=names(bcell_colors))


# Heatmap

#Norma data
data_to_plot <- rna_norm[tf_regulators$tf_ensembl,]
rownames(data_to_plot) <- rna_annotation[rownames(data_to_plot),"hgnc_symbol"]

# RNA color 
palette <- colorRampPalette(colors=c("#E5CCFF", "#4C0099","#FFFF00"))
cols <- palette(100)
my_rna_col <- colorRamp2(seq(-10,10, length.out=length(cols)),cols)

#rho color for annotation
col_fun1 = colorRamp2(c(-0.9,-0.8,-0.7,-0.6, 0, 0.6, 0.7, 0.8, 0.9), jdb_palette("solar_flare",9,"discrete"))


target_genes <- c("ETV6","SPIB","POU2F2","ELK3","E2F3","STAT5A","STAT5B","ETS2","ERG")
sel_pos <- vector("numeric",length=length(target_genes))
for(i in 1:length(target_genes)){
  sel_pos[i] <- which(rownames(data_to_plot)==target_genes[i])
}

ha = HeatmapAnnotation(stage = rna_metadata$CellType, 
                       col = list(stage = bcell_color),
                       rna = anno_barplot(as.numeric(rna_norm[rna_annotation$ensembl_gene_id[rna_annotation$hgnc_symbol=="ELK3"],])))


rha = rowAnnotation(rho=tf_regulators$gene_tf_rho,
                    col = list(rho = col_fun1),
                    genes=anno_mark(at = sel_pos, labels = target_genes))



Heatmap(data_to_plot, top_annotation = ha, right_annotation=rha, 
        show_row_names = TRUE, show_column_names =FALSE, 
        column_names_gp = gpar(fontsize = 5), row_names_gp = gpar(fontsize = 4),
        col=my_rna_col, 
        heatmap_legend_param = list(title = "Norm RNA"),
        column_order = order(rna_metadata$CellType),
        width = unit(5, "cm"), height = unit(10, "cm"),
        row_names_side = "left")

#***************************************************




#***************************************************
## Figure 3 E
#***************************************************
library(ggplot2)
library(ggridges)

# MLXIP exploration

# Input data:
# Normalized RNA: "rna_norm_data.txt" 
  # available at OSF files: 01_rna_seq_data/rna_norm_data.txt
# Metadata RNA: "rna_metadata.txt" 
  # available at OSF files: 01_rna_seq_data/rna_metadata.txt
# Gene annotation: "rna_gene_annotation.txt"
  # available at OSF files: 01_rna_seq_data/rna_gene_annotation.txt
# RNA expression colors: "rna_expression_colors.rds"


# RNA - Ridgeplot MLXIP

rna_norm <- read.table("rna_norm_data.txt",header=TRUE, dec=".", sep="\t")
rna_metadata <- read.table("rna_metadata.txt",header=TRUE, dec=".", sep="\t")
rna_annotation <- read.table("rna_gene_annotation.txt",header=TRUE, dec=".", sep="\t")
rna_color <- readRDS("rna_expression_colors.rds")

group <- factor(as.character(rna_metadata$CellType),
                levels=c("Naive_CD5pos","Naive_CD5neg","Transitional_B","ImmatureB","preB","proB","CLP","HSC"))

data_to_plot <- data.frame(cell_type=group, score=as.numeric(rna_norm[rna_annotation$ensembl_gene_id[rna_annotation$hgnc_symbol=="MLXIP"],]))
ggplot(data_to_plot, aes(x = score, y = cell_type, fill= ..x..)) +
  geom_density_ridges_gradient(scale=2) +
  scale_fill_gradientn(colors = rna_color) +
  labs(title = "MLXIP - RNA expression (Z-score)") +
  theme_ridges()


# TFBS - Ridgeplot MLXIP

tfbs_score <- read.table("TFBS_scores_from_ocr_tf_matrix_curated.txt",header=TRUE, dec=".", sep="\t")
atac_metadata <- read.table("atac_metadata.txt",header=TRUE, dec=".", sep="\t") 
tfbs_color <- readRDS("tfbs_color.rds")

group <- factor(as.character(atac_metadata$CellType),
                levels=c("Naive_CD5pos","Naive_CD5neg","Transitional_B","ImmatureB","preB","proB","CLP","HSC"))


data_to_plot <- data.frame(cell_type=group, score=as.numeric(tfbs_score["MLXIP",]))
ggplot(data_to_plot, aes(x = score, y = cell_type, fill= ..x..)) +
  geom_density_ridges_gradient(scale=2) +
  scale_fill_gradientn(colors = tfbs_color) +
  labs(title = "MLXIP - TFBS (Z-score)") +
  theme_ridges()

###-----------------------------------------------------------------------------------###
