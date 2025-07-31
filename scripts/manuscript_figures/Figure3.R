#***************************************************
## Paper: "Uncovering the regulatory landscape of early human B-cell lymphopoiesis and its implications in the pathogenesis of B-ALL"
## Authors: Planell N et al.
## Date: 2025
## Code: Figure 3
## Input data available at: https://osf.io/gswpy/?view_only=39cfd50d5a7e44e5a2d141869dcb1315
#***************************************************

# Define the directory were input data from OSF has been downloaded
data_wd <- setwd("../GitHub_code/")

# ****** Figure 3 ****** #

# Figure 3 - Panel A, B -----------------------------------------------------------------
# Box plot showing the distribution of ELK3 gene expression and ELK3 TFBS activity score by cell subpopulation.

# >>> Required packages
library(ggplot2)
library(ggsignif)

# >>> Input data
# Gene Expression data
rna <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/rna_norm_data.rds"))
rna_metadata <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/rna_metadata.rds"))
rna_anno <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/rna_gene_annotation.rds"))
# TFBS data
tfbs <- read.table(paste0(data_wd,"/osfstorage-archive/03_TFs_OCR_links/chromVAR_tf/TFBS_scores_from_ocr_tf_matrix.txt"))
atac_metadata <- readRDS(paste0(data_wd,"/osfstorage-archive/01_atac_seq_data/atac_metadata.rds"))

# >>> Plotting
# Plot colors
bcell_colors <- readRDS(paste0(data_wd,"/osfstorage-archive/00_plot_parameters/bcell_colors.rds"))

# Plot
data_to_plot <- data.frame(value=c(as.numeric(rna[rna_anno$ensembl_gene_id[rna_anno$hgnc_symbol=="ELK3"],]),as.numeric(tfbs["ELK3",])), 
                           gene=c(rep("ELK3 gene",nrow(rna_metadata)), rep("ELK3 TFBS", nrow(atac_metadata))),
                           cell_type=c(rna_metadata$CellType,atac_metadata$CellType))

pdf(paste0(format(Sys.time(), "%Y%m%d"),"_ELK3_rna_tfbs_stats.pdf"), width = 10, height = 5)
ggplot(data_to_plot, aes(x=cell_type, y=value, fill=cell_type)) + 
  geom_boxplot() +
  scale_fill_manual(values=bcell_colors)+
  ylim(c(-5,10)) +
  facet_wrap(~gene, nrow=1, scales="free_y") +
  geom_signif(comparisons = list(c("HSC", "CLP"),
                                 c("CLP", "proB"),
                                 c("proB", "preB"),
                                 c("preB", "ImmatureB"),
                                 c("ImmatureB", "Transitional_B"),
                                 c("Transitional_B", "Naive_CD5pos"),
                                 c("Transitional_B", "Naive_CD5neg")), 
              map_signif_level=TRUE,
              y_position = c(7.5, 8, 9, 7,5,6.5,8),
              tip_length = 0) +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) 
dev.off()
# ------------------------------------------------------------------------------------ #



# Figure 3 - Panel C -----------------------------------------------------------------
# Heatmap representation of gene expression for TF regulators of ELK3. 

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
target_gene <- "ELK3"
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
target_genes <- c("SPIB","POU2F2","JUNB","IRF4","ELK3","ERG","STAT5A","ETS2","ETV6","ZNF467","E2F3")
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



# Figure 3 - Panel D -----------------------------------------------------------------
# Heatmap representation of gene expression for TF regulated by ELK3. 

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
#Genes regulated by ELK3 according to our network.
target_gene <- "ELK3"
list_target_genes_ELK3 <- unique(tf_ocr_gene[which(tf_ocr_gene$tf_symbol==target_gene),"gene_ensembl"])
list_target_tfs_ELK3 <- list_target_genes_ELK3[list_target_genes_ELK3%in%unique(rna_anno$ensembl_gene_id[rna_anno$hgnc_symbol%in%colnames(ocr_tf_matrix)])]
target_grn <- unique(tf_ocr_gene[which(tf_ocr_gene$tf_symbol==target_gene & tf_ocr_gene$gene_ensembl%in%list_target_tfs_ELK3),c("gene_ensembl","tf_ensembl","gene_tf_rho")])


# Data to plot
M <- rna[target_grn$gene_ensembl,]
rownames(M) <- rna_anno[rownames(M),"hgnc_symbol"]

# >>> Plotting
# Plot colors
bcell_colors <- readRDS(paste0(data_wd,"/osfstorage-archive/00_plot_parameters/bcell_colors.rds"))
col_fun1 = colorRamp2(c(-0.9,-0.8,-0.7,-0.6, 0, 0.6, 0.7, 0.8, 0.9), jdb_palette("solar_flare",9,"discrete"))

# Plot
# Row annotation
target_genes <- c("IKZF1", "SPI1", "POU2F2", "BCL11A", "STAT5A", "IRF4", "SPIB", "PAX5", "MEIS1")
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



# Figure 3 - Panel E, F -----------------------------------------------------------------
# Representation of the main enriched GO terms derived from the pathway analysis of the ELK3-regulon-activated and ELK3-regulon-repressed genes for each csGRN. 

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
#Genes regulated by ELK3 according to our network.
target_gene <- "ELK3"
list_target_genes_ELK3 <- unique(tf_ocr_gene[which(tf_ocr_gene$tf_symbol==target_gene),])

# Split by positive or negative regulation
list_target_genes_ELK3_pos <- unique(list_target_genes_ELK3[list_target_genes_ELK3$gene_tf_rho>0,"gene_ensembl"])
list_target_genes_ELK3_neg <- unique(list_target_genes_ELK3[list_target_genes_ELK3$gene_tf_rho<0,"gene_ensembl"])

# ELK3-regulon ORA analysis from csGRN
list_target_genes_activator_ELK3_by_cell_type <- vector(mode="list", length=length(bcell_csGRN))
names(list_target_genes_activator_ELK3_by_cell_type) <- names(bcell_csGRN)

list_target_genes_repressor_ELK3_by_cell_type <- vector(mode="list", length=length(bcell_csGRN))
names(list_target_genes_repressor_ELK3_by_cell_type) <- names(bcell_csGRN)

for(i in names(bcell_csGRN)){
  aa <- as_edgelist(bcell_csGRN[[i]])
  bb <- unique(aa[aa[,1]==rna_anno$ensembl_gene_id[rna_anno$hgnc_symbol==target_gene],2])
  if(length(bb)>0){
    list_target_genes_activator_ELK3_by_cell_type[[i]] <- bb[bb%in%list_target_genes_ELK3_pos]
    list_target_genes_repressor_ELK3_by_cell_type[[i]] <- bb[bb%in%list_target_genes_ELK3_neg]
  }
}

str(list_target_genes_activator_ELK3_by_cell_type)
str(list_target_genes_repressor_ELK3_by_cell_type)

# activator ELK3-regulon per csGRN
ck_go_005_activator <- compareCluster(geneCluster = list_target_genes_activator_ELK3_by_cell_type, 
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

# repressor ELK3-regulon per csGRN
ck_go_005_repressor <- compareCluster(geneCluster = list_target_genes_repressor_ELK3_by_cell_type, 
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

# repressor ELK3-regulon per proB cells, the only cell type with enriched pathways.
ora_GO_ELK3_neg_proB <- enrichGO(gene          = list_target_genes_repressor_ELK3_by_cell_type[["proB"]],
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

# >>> Plotting
pdf(paste0(format(Sys.time(), "%Y%m%d"),"_ora_GO_activator_ELK3_regulon_from_csGRNs_dotplot_targeted.pdf"))
  dotplot(ck_go_005_activator, size="Count", showCategory =c("unsaturated fatty acid metabolic process",
                                                     "cellular oxidant detoxification",
                                                     "small molecule biosynthetic process",
                                                     "positive regulation of phosphatidylinositol 3-kinase/protein kinase B signal transduction",
                                                     "positive regulation of MAPK cascade",
                                                     "positive regulation of cell cycle",
                                                     "DNA replication",
                                                     "ERK1 and ERK2 cascade",
                                                     "nucleoside monophosphate biosynthetic process",
                                                     "positive regulation of cell cycle process",
                                                     "cell cycle DNA replication",
                                                     "nuclear division",
                                                     "chromosome segregation",
                                                     "regulation of cell cycle phase transition",
                                                     "cell cycle G2/M phase transition")) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf(paste0(format(Sys.time(), "%Y%m%d"),"_ora_GO_ELK3_regulon_from_csGRNs_dotplot_proB_targeted_neg.pdf"), width = 7)
  dotplot(ora_GO_ELK3_neg_proB, size="Count", showCategory =c("immune response-regulating cell surface receptor signaling pathway",
                                                            "B cell activation",
                                                            "B cell differentiation",
                                                            "mononuclear cell differentiation",
                                                            "regulation of lymphocyte activation",
                                                            "lymphocyte differentiation",
                                                            "B cell receptor signaling pathway",
                                                            "Fc receptor signaling pathway",
                                                            "B cell proliferation",
                                                            "mature B cell differentiation")) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
# ------------------------------------------------------------------------------------ #
