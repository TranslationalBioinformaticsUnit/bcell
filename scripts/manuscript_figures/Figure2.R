#***************************************************
## Paper: "Uncovering the regulatory landscape of early human B-cell lymphopoiesis and its implications in the pathogenesis of B-ALL"
## Authors: Planell N et al.
## Date: 2025
## Code: Figure 2
## Input data available at: https://osf.io/gswpy/?view_only=39cfd50d5a7e44e5a2d141869dcb1315
#***************************************************

# Define the directory were input data from OSF has been downloaded
data_wd <- setwd("../GitHub_code/")

# ****** Figure 2 ****** #

# Figure 2 - Panel A/B -----------------------------------------------------------------
#  Regulon activity, based on ATAC-Seq accessibility score based on GSVA and TF expression.
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
source("auxiliar_functions.R")

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

# >>> Computing
# Get the list of TFs in the B-cell GRN
tf_regulons <- unique(tf_ocr_gene[,c("tf_symbol","tf_ensembl")])

# Compute the accessibility of the atac for each regulon spliting by positive/negative correlation
gsva_tf_regulons_pos_mean <- matrix(NA, nrow=nrow(tf_regulons), ncol=length(levels(rna_metadata$CellType)),
                                    dimnames = list(tf_regulons$tf_ensembl, levels(rna_metadata$CellType)))
gsva_tf_regulons_neg_mean <- matrix(NA, nrow=nrow(tf_regulons), ncol=length(levels(rna_metadata$CellType)),
                                    dimnames = list(tf_regulons$tf_ensembl, levels(rna_metadata$CellType)))
gsva_tf_regulons_pos <- matrix(NA, nrow=nrow(tf_regulons), ncol=ncol(atac),
                               dimnames = list(tf_regulons$tf_ensembl, colnames(atac)))
gsva_tf_regulons_neg <- matrix(NA, nrow=nrow(tf_regulons), ncol=ncol(atac),
                               dimnames = list(tf_regulons$tf_ensembl, colnames(atac)))
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
  gs <- list(pos=unique(target_tf_regulon_pos$ocr),
             neg=unique(target_tf_regulon_neg$ocr))
  gsvaPar <- gsvaParam(as.matrix(atac), gs)
  gsva.es <- gsva(gsvaPar, verbose=FALSE)
  gsva_mean_by_celltype <- apply(gsva.es,1,FUN=function(x){
    mm <- by(x,atac_metadata$CellType,mean)
    return(as.numeric(mm))
  })
  rownames(gsva_mean_by_celltype) <- levels(atac_metadata$CellType)
  gsva_tf_regulons_pos_mean[i,] <- gsva_mean_by_celltype[,"pos"]
  gsva_tf_regulons_neg_mean[i,] <- gsva_mean_by_celltype[,"neg"]
  gsva_tf_regulons_pos[i,] <- gsva.es["pos",]
  gsva_tf_regulons_neg[i,] <- gsva.es["neg",]
}

# Defined thresholds: 0.2
## Binarize the atac_seq matrix

gsva_tf_regulons_pos_binary <- matrix(0, nrow=nrow(tf_regulons), ncol=ncol(atac),
                                      dimnames = list(tf_regulons$tf_ensembl, colnames(atac)))

gsva_tf_regulons_pos_binary[gsva_tf_regulons_pos>0.2] <- 1
rownames(gsva_tf_regulons_pos_binary) <- tf_regulons$tf_symbol

gsva_tf_regulons_neg_binary <- matrix(0, nrow=nrow(tf_regulons), ncol=ncol(atac),
                                      dimnames = list(tf_regulons$tf_ensembl, colnames(atac)))

gsva_tf_regulons_neg_binary[gsva_tf_regulons_neg>0.2] <- 1
rownames(gsva_tf_regulons_neg_binary) <- tf_regulons$tf_symbol

# Calculates Regulon specificity score (RSS)
rrs_pos <- calculate_rss(metadata=atac_metadata, binary_regulons=gsva_tf_regulons_pos_binary, cell_type_column="CellType")
rrs_neg <- calculate_rss(metadata=atac_metadata, binary_regulons=gsva_tf_regulons_neg_binary, cell_type_column="CellType")

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
set.seed(123)
rrs_matrix <- matrix(nrow=nrow(data_to_plot2), ncol=ncol(data_to_plot2), 
                     dimnames = list(rownames(data_to_plot2), colnames(data_to_plot2)))

for(i in colnames(rrs_matrix)){
  tmp <- rrs_pos[rrs_pos$cell_type==i,]
  rownames(tmp) <- tmp$regulon
  rrs_matrix[,i] <- tmp[rownames(rrs_matrix),"RSS"]
}

pp <- Heatmap(data_to_plot2, 
              name= "Normalized rna",
              cluster_columns = FALSE,
              top_annotation = HeatmapAnnotation(cell = as.character(colnames(data_to_plot2)),
                                                 col = list(cell = bcell_colors)),
              col=jdb_palette("solar_extra"),
              row_km=5,
              row_names_gp = gpar(fontsize = 6),
              show_column_names = FALSE,
              right_annotation = rowAnnotation(foo = anno_empty(border = FALSE, 
                                                                width = unit(4, "mm"))),
              layer_fun = function(j, i, x, y, w, h, fill) {
                ind_mat = restore_matrix(j, i, x, y)
                v = pindex(rrs_matrix, i, j)
                vpos0 = v <= 0
                vpos01 = v>0 & v <= 0.1
                vpos02 = v>0.1 & v <= 0.2
                vpos03 = v>0.2 & v <= 0.3
                vpos04 = v>0.3 
                grid.points(x[ind_mat[vpos0]],y[ind_mat[vpos0]], pch = 19, size = unit(0.1, "mm"))
                grid.points(x[ind_mat[vpos01]],y[ind_mat[vpos01]], pch = 19, size = unit(0.3, "mm"))
                grid.points(x[ind_mat[vpos02]],y[ind_mat[vpos02]], pch = 19, size = unit(0.8, "mm"))
                grid.points(x[ind_mat[vpos03]],y[ind_mat[vpos03]], pch = 19, size = unit(1, "mm"))
                grid.points(x[ind_mat[vpos04]],y[ind_mat[vpos04]], pch = 19, size = unit(1.5, "mm"))
              }
)

ht_pos = draw(pp)
rgroups <- row_order(ht_pos)
rgroups_named <- lapply(rgroups,FUN=function(x){return(rownames(data_to_plot2)[x])})

# Include the cluster info in the formatted data to plot (data1)
for(i in unique(data1$Var1)){
  group <- lapply(rgroups_named,FUN=function(x){i%in%x})
  data1$group[data1$Var1==i] <- names(group)[group==TRUE]
}
data1$Var1 <- factor(data1$Var1, levels=unlist(rgroups_named))

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

ggsave(plot = ht_pos, width = 4, height = 14, filename = paste0("RSS_tf_exprss_pos.pdf")) 

# RSS tf-gene negatively correlated (putative repressor TFs)
ht_neg <- ggplot(data=data1, mapping=aes(x=Var2, y=Var1)) + 
  geom_raster(aes(fill=value)) +
  scale_fill_gradientn(colors=jdb_palette("brewer_celsius")) +
  geom_point(mapping=aes(size=rrs_neg), color="black") +
  scale_radius(range=c(min.size=0.1, max.size=2)) +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.title.y=element_blank(), 
        axis.text.x=element_text(angle=90, hjust=1, size=6),
        axis.text.y=element_text(size=6))

ggsave(plot = ht_neg, width = 4, height = 14, filename = paste0("RSS_tf_exprss_neg.pdf")) 

# Select TFs with coordinated activator or repressor activity
sel_pos <- vector()
sel_neg <- vector()
for(i in tf_regulons$tf_symbol){
  dd <- data1[data1$Var1==i,]
  a <- cor.test(dd$rrs_pos,dd$value, method="spearman")
  if(a$estimate>0){
    sel_pos <- c(sel_pos,i)
  }
  a <- cor.test(dd$rrs_neg,dd$value, method="spearman")
  if(a$estimate>0){
    sel_neg <- c(sel_neg,i)
  }
}

length(sel_neg); length(sel_pos)

# RSS tf-gene negatively correlated (putative repressor TFs) from TFexpression and RSS cooridnation

P <- ggplot(data=data1[data1$Var1%in%sel_neg,], mapping=aes(x=Var2, y=Var1)) + 
  geom_raster(aes(fill=value)) +
  scale_fill_gradientn(colors=jdb_palette("brewer_celsius")) +
  geom_point(mapping=aes(size=rrs_neg), color="black") +
  scale_radius(range=c(min.size=0.1, max.size=2)) +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.title.y=element_blank(), 
        axis.text.x=element_text(angle=90, hjust=1, size=6),
        axis.text.y=element_text(size=6))

ggsave(plot = P, width = 4, height = 5, filename = paste0("RSS_tf_exprss_neg_sig.pdf")) 
# ------------------------------------------------------------------------------------ #



# Figure 2 - Panel C -----------------------------------------------------------------
# Pathway analysis of regulons.
# Circus plot summarizing the main biological pathways enriched for each regulon.

# >>> Required packages
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(biomaRt)
library(circlize)

# >>> Input data
# Gene Expression data
rna <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/rna_norm_data.rds"))
# GRN data
tf_ocr_gene <- read.table(paste0(data_wd,"/osfstorage-archive/04_gene_regulatory_networks/tf_ocr_gene_interactions_curated.txt"),
                          sep="\t",
                          dec=",",
                          header=TRUE)
# Regulons clusters
regulon_clusters <- readRDS(paste0(data_wd,"/osfstorage-archive/04_gene_regulatory_networks/regulons/TF_regulons_clusters.rds"))


# >>> Computing
### ORA analysis of full list of regulons
# Insted of run this section you can upload the following file:
# "/osfstorage-archive/04_gene_regulatory_networks/regulons/list_regulon_ORA_GO.rds"
# "/osfstorage-archive/04_gene_regulatory_networks/regulons/list_regulon_ORA_KEGG.rds"
# corresponding to list_regulon_ORA_GO and list_regulon_ORA_KEGG, respectively.

# Retrieve the ensembl info
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "dec2021.archive.ensembl.org")
entrez_IDS_universe <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", 'entrezgene_id', 'entrezgene_accession'), values = rownames(rna), mart = human)

# List of B-cell differentiation regulons
target_tfs <- unique(tf_ocr_gene$tf_symbol)  #169

list_regulon_ORA_GO <- vector(mode="list", length=length(target_tfs))
names(list_regulon_ORA_GO) <- target_tfs

list_regulon_ORA_KEGG <- vector(mode="list", length=length(target_tfs))
names(list_regulon_ORA_KEGG) <- target_tfs

for(i in target_tfs){
  cat(">> ",i," ....\n")
  #Genes regulated by target TF according to our network.
  list_target_genes <- unique(tf_ocr_gene[which(tf_ocr_gene$tf_symbol==i),"gene_ensembl"])
  length(list_target_genes)
  
  ## ORA based on GO terms with all the genes linked to ELK3
  list_regulon_ORA_GO[[i]] <- enrichGO(gene          = list_target_genes,
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
  
  cat("ORA GO done\n")  
  ## ORA based on KEGG terms with all the genes linked to ELK3
  # Filter the ensembl ID
  entrez_IDS_list_target_genes <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id", 'entrezgene_id', 'entrezgene_accession'), values = list_target_genes, mart = human)
  
  list_regulon_ORA_KEGG[[i]] <- enrichKEGG(gene         = entrez_IDS_list_target_genes$entrezgene_id,
                                           organism     = 'hsa',
                                           universe      = as.character(entrez_IDS_universe$entrezgene_id),
                                           minGSSize = 20,
                                           pvalueCutoff = 0.05)
  cat("ORA KEGG done\n")  
}

# >>> Plotting
# Select categories:
target_categories <- c("Immune system","Environmental adaptation","Cell growth and death","Replication and repair","Global and overview maps",
                       "Development and regeneration","Metabolism of other amino acids","Signal transduction","Signaling molecules and interaction","Transport and catabolism",
                       "Sensory system","Circulatory system","Nucleotide metabolism","Xenobiotics biodegradation and metabolism","Cellular community - eukaryotes",
                       "Carbohydrate metabolism","Metabolism of cofactors and vitamins","Cell motility","Lipid metabolism","Chromosome","Amino acid metabolism","Translation")

# Full list of pathways:
overall_enrich_description <- lapply(list_regulon_ORA_KEGG, FUN=function(x){return(unique(x@result[which(x@result$p.adjust<0.05 & x@result$subcategory%in%target_categories),"Description"]))})
names(overall_enrich_description) <- names(list_regulon_ORA_KEGG)

tfs_enrich <- vector()
pathways_enrich <- vector()

target_pathways <- c("VEGF signaling pathway","Tight junction", "Inositol phosphate metabolism","Apoptosis","ErbB signaling pathway","Hedgehog signaling pathway",
                     "Adherens junction","Phosphatidylinositol signaling system","PI3K-Akt signaling pathway",
                     "Chemokine signaling pathway","Phospholipase D signaling pathway",
                     "Fc gamma R-mediated phagocytosis","Fc epsilon RI signaling pathway",
                     "Gap junction","JAK-STAT signaling pathway",
                     "NF-kappa B signaling pathway","MAPK signaling pathway",
                     "Wnt signaling pathway","HIF-1 signaling pathway","C-type lectin receptor signaling pathway",
                     "cGMP-PKG signaling pathway","cAMP signaling pathway",
                     "Calcium signaling pathway","Leukocyte transendothelial migration","NOD-like receptor signaling pathway",
                     "Rap1 signaling pathway", "Necroptosis",
                     "FoxO signaling pathway","Cell adhesion molecules","Phagosome",
                     "Hematopoietic cell lineage","B cell receptor signaling pathway","Cellular senescence",
                     "Antigen processing and presentation","DNA replication",
                     "Cell cycle","p53 signaling pathway",
                     "Pyrimidine metabolism","Notch signaling pathway",
                     "Hippo signaling pathway","Nucleotide metabolism",
                     "Ribosome","Glutathione metabolism",
                     "Fatty acid metabolism","Fatty acid elongation","Base excision repair","Mismatch repair",
                     "Polycomb repressive complex","Purine metabolism",
                     "Lysine degradation","Glycerolipid metabolism","Cysteine and methionine metabolism")

for(i in 1:length(overall_enrich_description)){
  if(length(overall_enrich_description[[i]])>0){
    tfs_enrich <- c(tfs_enrich,rep(names(overall_enrich_description)[i],sum(overall_enrich_description[[i]]%in%target_pathways)))
    pathways_enrich <- c(pathways_enrich,overall_enrich_description[[i]][overall_enrich_description[[i]]%in%target_pathways])
  }
}

data_to_plot <- data.frame(tfs=tfs_enrich, pathway=pathways_enrich)

for(i in unique(data_to_plot$tfs)){
  group <- lapply(regulon_clusters,FUN=function(x){i%in%x})
  data_to_plot$cell_type[data_to_plot$tfs==i] <- names(group)[group==TRUE]
}

data_to_plot$cell_type <- factor(data_to_plot$cell_type, levels=names(regulon_clusters))

data_to_plot$pathway <- factor(data_to_plot$pathway, 
                               levels=rev(target_pathways))

data_to_plot$tfs <- factor(data_to_plot$tfs, levels=c(rev(regulon_clusters[[1]]),
                                                      rev(regulon_clusters[[2]]),
                                                      rev(regulon_clusters[[3]]),
                                                      rev(regulon_clusters[[4]]),
                                                      rev(regulon_clusters[[5]])))

mat <- table(data_to_plot$tfs, data_to_plot$pathway)

# Plotting colors
grid.col <- setNames(c(rep("#00FFFF",length(regulon_clusters[[1]])),
                       rep("#1E90FF",length(regulon_clusters[[2]])),
                       rep("#32CD32",length(regulon_clusters[[3]])),
                       rep("#DC143C",length(regulon_clusters[[4]])),
                       rep("#252525",length(regulon_clusters[[5]])),
                       rep("#fa9fb5",10),rep("#756bb1",5),rep("#fec44f",9),rep("#a6bddb",13), rep("#1c9099",15)), union(rownames(mat), colnames(mat)))

# Plot
pdf(paste0(format(Sys.time(), "%Y%m%d"),"_circus_plot_regulons_complete4.pdf")) #, height = 20, width=20)
  chordDiagram(mat, annotationTrack = "grid", preAllocateTracks = 1, grid.col = grid.col)
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.name = get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex=0.3)
    circos.axis(h = "top", labels.cex = 0.5, major.tick.percentage = 0.2, sector.index = sector.name, track.index = 2)
  }, bg.border = NA)
dev.off()

###-----------------------------------------------------------------------------------###



# Figure 2 - Panel D -----------------------------------------------------------------
# Multidimensional Scaling representation from the Pearson correlation of cell subpopulation specific GRNs. 

# >>> Required packages
library(igraph)
library(ComplexHeatmap)
library(ggplot2)

# >>> Required auxiliar functions
source("auxiliar_functions.R")

# >>> Input data
# GRN data
bcell_GRN <- readRDS(paste0(data_wd,"/osfstorage-archive/04_gene_regulatory_networks/bcell_grn.rds"))
bcell_csGRN <- readRDS(paste0(data_wd,"/osfstorage-archive/04_gene_regulatory_networks/bcell_grn_by_cell_type.rds"))

# >>> Computing and plotting
# Plot colors
bcell_colors <- readRDS(paste0(data_wd,"/osfstorage-archive/00_plot_parameters/bcell_colors.rds"))

# plot
pearson_mds(GRN=bcell_GRN, csGRN=bcell_csGRN, name="B-cells", pca_colors=bcell_colors)

# ------------------------------------------------------------------------------------ #



# Figure 2 - Panel E -----------------------------------------------------------------
# Betweenness and degree GRN by cell type plot.
# Heatmap representation of the degree (left) and betweenness (right) 
# centrality measures computed within each csGRN. 

# >>> Required packages
library(igraph)
library(ComplexHeatmap)
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
bcell_GRN <- readRDS(paste0(data_wd,"/osfstorage-archive/04_gene_regulatory_networks/bcell_grn.rds"))
bcell_csGRN <- readRDS(paste0(data_wd,"/osfstorage-archive/04_gene_regulatory_networks/bcell_grn_by_cell_type.rds"))
# Regulons clusters
regulon_clusters <- readRDS(paste0(data_wd,"/osfstorage-archive/04_gene_regulatory_networks/regulons/TF_regulons_clusters.rds"))
                        
# >>> Computing

# Centrality measures:

# Whole GRN
grn <- upgrade_graph(bcell_GRN)
V(grn)$degree <- igraph::degree(grn)        # Degree centrality
V(grn)$betweenness <- betweenness(grn)      # Vertex betweenness centrality

# GRN specific per cell type (csGRN)
for(i in 1:length(bcell_csGRN)){
  # Centrality measures
  V(bcell_csGRN[[i]])$degree <- igraph::degree(bcell_csGRN[[i]])      # Degree centrality
  V(bcell_csGRN[[i]])$betweenness <- betweenness(bcell_csGRN[[i]])    # Vertex betweenness centrality
}

# Data frame with centrality measures
centrality <- data.frame(row.names   = V(grn)$name,						
                         degree      = V(grn)$degree,
                         betweenness = V(grn)$betweenness)

for(i in 1:length(bcell_csGRN)){
  centrality2 <- data.frame(row.names   = V(bcell_csGRN[[i]])$name,						
                            degree      = V(bcell_csGRN[[i]])$degree,
                            betweenness = V(bcell_csGRN[[i]])$betweenness)
  colnames(centrality2) <- paste0(colnames(centrality2),"_",names(bcell_csGRN)[i])
  centrality <- merge(centrality,centrality2, by.x=0, by.y=0, all=TRUE)
  rownames(centrality) <- centrality$Row.names
  centrality$Row.names <- NULL
}

centrality_tf <- centrality[unique(tf_ocr_gene$tf_ensembl),]
rownames(centrality_tf) <- rna_anno[rownames(centrality_tf),"hgnc_symbol"]

# >>> Plotting

# Plot colors
bcell_colors <- readRDS(paste0(data_wd,"/osfstorage-archive/00_plot_parameters/bcell_colors.rds"))

# Plot
# Get regulons clustering
regulon_groups <- vector()
for(i in 1:length(regulon_clusters)){
  regulon_groups <- c(regulon_groups,rep(names(regulon_clusters)[i], length(regulon_clusters[[i]])))
}
regulon_groups <- factor(regulon_groups,levels=rev(names(regulon_clusters)))


# Get top10 tf degree per cell type
data_to_explore <- centrality_tf[unlist(regulon_clusters),grep("degree_",colnames(centrality_tf),value=TRUE)]
top10_degree <- vector()
for(i in 1:ncol(data_to_explore)){
  top10_degree <- c(top10_degree,rownames(data_to_explore)[order(data_to_explore[,i], decreasing = TRUE)[1:10]])
}

data_to_explore <- centrality_tf[unlist(regulon_clusters),grep("betweenness_",colnames(centrality_tf),value=TRUE)]
top10_betweenness <- vector()
for(i in 1:ncol(data_to_explore)){
  a <- rownames(data_to_explore)[order(data_to_explore[,i], decreasing = TRUE)[1:10]]
  top10_betweenness <- c(top10_betweenness,a)
}

top10_to_highlight <- unique(c(top10_degree,top10_betweenness))

# Degree representation of cell type specific GRN
grn_properties_degree_cell_type2 <- Heatmap(centrality_tf[rev(unlist(regulon_clusters)),grep("degree_",colnames(centrality_tf),value=TRUE)], 
                                            row_names_gp = gpar(fontsize = 5), 
                                            width = unit(4, "cm"),
                                            border = "black",
                                            left_annotation = rowAnnotation(foo = anno_empty(border = FALSE, 
                                                                                             width = unit(4, "mm"))),
                                            col=jdb_palette("flame_light"),
                                            name="degree",
                                            row_names_side = "left",
                                            top_annotation=HeatmapAnnotation(stage = names(bcell_colors),
                                                                             col = list(stage = bcell_colors)),
                                            row_split = rev(regulon_groups),
                                            na_col = "white", 
                                            cluster_rows=FALSE, 
                                            cluster_columns=FALSE, 
                                            show_column_names = FALSE)

# Betweeness representation of cell type specific GRN
grn_properties_betweenness_cell_type2 <- Heatmap(centrality_tf[rev(unlist(regulon_clusters)),grep("betweenness_",colnames(centrality_tf),value=TRUE)], 
                                                 row_names_gp = gpar(fontsize = 5), 
                                                 width = unit(4, "cm"),
                                                 border = "black",
                                                 col=jdb_palette("flame_artic"), name="betweenness",
                                                 top_annotation=HeatmapAnnotation(stage = names(bcell_colors),
                                                                                  col = list(stage = bcell_colors)),
                                                 na_col = "white",
                                                 cluster_rows=FALSE, 
                                                 cluster_columns=FALSE, 
                                                 show_column_names = FALSE,
                                                 show_row_names = FALSE)

pdf(paste(format(Sys.time(), "%Y%m%d"),"_grn_cell_type_properties_weighted_clusters.pdf",sep=""), width=15, height=15)
grn_properties_degree_cell_type2+grn_properties_betweenness_cell_type2
cluster_color <- c()
for(i in 1:5) {
  decorate_annotation("foo", slice = i, {
    grid.rect(x = 0, width = unit(2, "mm"), gp = gpar(fill = 6-i, col = NA), just = "left")
  })
}
dev.off()

# ------------------------------------------------------------------------------------ #
