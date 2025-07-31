#***************************************************
## Paper: "Uncovering the regulatory landscape of early human B-cell lymphopoiesis and its implications in the pathogenesis of B-ALL"
## Authors: Planell N et al.
## Date: 2025
## Code: Figure 5
## Input data available at: https://osf.io/gswpy/?view_only=39cfd50d5a7e44e5a2d141869dcb1315
#***************************************************

# Define the directory were input data from OSF has been downloaded
data_wd <- setwd("../GitHub_code/")

# ****** Figure 5 ****** #

# B-ALL Public datasets analysis -------------------------------------------------

###-----------------------------------------------------------------------------------###
# ****** Li J. multi-dataset (EGAD00001002704, EGAD00001002692, EGAD00001002151) ****** #
###  Li J. multi-dataset ----------

# >>> Input data
# B-ALL Gene Expression data
ball_rna_Li <- read.csv("Li_dataset/raw_counts.csv", row.names = 1)
ball_metadata_Li <- read.csv("Li_dataset/metadata.csv")

# Subset specific mutations
ball_metadata_Li <- ball_metadata_Li[ball_metadata_Li$subgroup%in%c("G1","G2","G3","G4","G5","G6","G7","G8"),]
ball_rna_Li <- ball_rna_Li[,ball_metadata_Li$samples_ID]

ball_metadata_Li$mutation <- recode_factor(ball_metadata_Li$subgroup, 
                                      G1="MEF2D fusions",
                                      G2="TCF3-PBX1",
                                      G3="ETV6-RUNX1/like",
                                      G4="DUX4 fusions",
                                      G5="ZNF384 fusions",
                                      G6="BCR-ABL1/ph-like",
                                      G7="Hyperdiploidy",
                                      G8="KTM2A fusions")

# >>>  Data normalization
require(edgeR)
# Create DEGlist object with count data and defining the grouping of the variables
dge <- DGEList(counts=ball_rna_Li, group=ball_metadata_Li$subgroup)

# Filter out low expressed genes
filtering_criteria_edgeR <- filterByExpr(dge, group =ball_metadata_Li$subgroup,  min.count = 5, min.total.count = 10)
dge.filt <- dge[filtering_criteria_edgeR,]
dim(dge.filt)

# Define the model matrix
dge.filt <- calcNormFactors(dge.filt) # computes the correction factor based on TMM (default method)
group <- as.factor(ball_metadata_Li$subgroup)
design <- model.matrix( ~ 0 + group)
colnames(design) <- gsub("group","", colnames(design))

# Voom transformation including log2 CPM normalization
v <- voom(dge.filt, design = design, plot = TRUE)
ball_norm_rna_Li <- v$E

# Correct batch effect
library(sva)
ball_norm_rna_Li <- ComBat(dat = ball_norm_rna_Li, batch = ball_metadata_Li$batch, mod = NULL, par.prior = TRUE, prior.plots = FALSE)


# >>> Exploratory analysis
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


# >>> Differential expression analysis:

DEA_output <- list()
DEA_BALL_biomarkers <- list()
min_contrasts <- 5
for(target_group in levels(group)){
  vv <- paste0(target_group,"-",colnames(design)[!colnames(design)%in%target_group])
  # To get the changes for each comparison we need to define a contrast matrix
  fit <- lmFit(ball_norm_rna_Li, design)
  contr <- makeContrasts(vv[1],
                         vv[2],
                         vv[3],
                         vv[4],
                         vv[5],
                         vv[6],
                         vv[7],
                         levels = colnames(coef(fit)))
  tmp <- contrasts.fit(fit, contr)
  fit <- eBayes(tmp)
  
  # Results for first contrast:
  DEA_output[[target_group]] <- list()
  for(i in vv){
    tmp_output <- topTable(fit, coef=i,number=nrow(ball_norm_rna_Li), sort.by="none")
    DEA_output[[target_group]][[i]] <- tmp_output
  }
  sig_genes <- lapply(DEA_output[[target_group]],FUN=function(x){return(rownames(x)[x$adj.P.Val<0.01 & x$logFC>0.58])})
  tt <- table(unlist(sig_genes))
  DEA_BALL_biomarkers[[target_group]] <- names(tt)[tt>min_contrasts]
}  
str(DEA_BALL_biomarkers)  
names(DEA_BALL_biomarkers) <- recode(names(DEA_BALL_biomarkers), 
                                 G1="MEF2D fusions",
                                 G2="TCF3-PBX1",
                                 G3="ETV6-RUNX1/like",
                                 G4="DUX4 fusions",
                                 G5="ZNF384 fusions",
                                 G6="BCR-ABL1/ph-like",
                                 G7="Hyperdiploidy",
                                 G8="KTM2A fusions")

saveRDS(DEA_BALL_biomarkers,file="ALL_subgroup_biomarkers.rds")


# >>> Computing GSVA
# GRN data
tf_ocr_gene <- read.table(paste0(data_wd,"/osfstorage-archive/04_gene_regulatory_networks/tf_ocr_gene_interactions_curated.txt"),
                          sep="\t",
                          dec=",",
                          header=TRUE)

### GSVA
# Get the TF regulons gene sets for the positively regulated genes.
df_regulons_pos <- unique(tf_ocr_gene[tf_ocr_gene$likelihood=="activator",c("tf_symbol","gene_ensembl")])

df_regulons_list_pos <- list()
for(i in unique(df_regulons_pos$tf_symbol)){
  df_regulons_list_pos[[i]] <- df_regulons_pos$gene_ensembl[df_regulons_pos$tf_symbol==i]
}

gs <- df_regulons_list_pos
gsvaPar <- gsvaParam(ball_norm_rna_Li, gs, kcdf="Gaussian")
gsva.BALL_Li <- gsva(gsvaPar, verbose=FALSE)


# >>> Save results
save(ball_norm_rna_Li,ball_metadata_Li,DEA_output,gsva.BALL_Li,file="BALL_Li.RData")
###-----------------------------------------------------------------------------------###


###-----------------------------------------------------------------------------------###
# ****** Rainer dataset (GSE73578) ****** #
###  Rainer dataset ----------
# Microarray-based whole-genome expression profiling of lymphoblasts from 46 patients

library(Biobase)
library(GEOquery)

# >>> Input data
# B-ALL Gene Expression data
# Get the data from GEO database
gse <- getGEO("GSE73578", GSEMatrix = TRUE)
show(gse)
gse <- gse[[1]]

ball_norm_rna_Rainer <- exprs(gse)
ball_metadata_Rainer <- pData(gse)

# Subset specific mutations
ball_metadata_Rainer <- ball_metadata_Rainer[ball_metadata_Rainer$"all subtype by clustering:ch1"%in%c("BCR/ABL","E2A/PBX1","ETV6/RUNX1","Hyperdiploid"),]
ball_norm_rna_Rainer <- ball_norm_rna_Rainer[,rownames(ball_metadata_Rainer)]

ball_metadata_Rainer$mutation <- recode_factor(ball_metadata_Rainer$"all subtype by clustering:ch1",
                                           "E2A/PBX1"="TCF3-PBX1",
                                           "ETV6/RUNX1"="ETV6-RUNX1/like",
                                           "BCR/ABL"="BCR-ABL1/ph-like",
                                           "Hyperdiploid"="Hyperdiploidy")

ball_metadata_Rainer$timepoint <- factor(recode(ball_metadata_Rainer$"treatment time point:ch1",
                                               "6"="6/8",
                                               "8"="6/8"), levels=c("0","6/8","24"))

# Translate probe ID to gene names:
## Exclude genes with multiple annotations
gene_name_Rainer <- fData(gse)[,"Gene Symbol"]
to_exclude <- grep("/",gene_name_Rainer)
length(to_exclude)
ball_norm_rna_Rainer <- ball_norm_rna_Rainer[-to_exclude,]
rownames(ball_norm_rna_Rainer) <- gene_name_Rainer[-to_exclude]


# >>> Exploratory analysis
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
gsvaPar <- gsvaParam(ball_norm_rna_Rainer, gs_symbol, kcdf="Gaussian")
gsva.BALL_Rainer <- gsva(gsvaPar, verbose=FALSE)


# >>> Save results
save(ball_norm_rna_Rainer,ball_metadata_Rainer,gsva.BALL_Rainer,file="BALL_Rainer.RData")
###-----------------------------------------------------------------------------------###


###-----------------------------------------------------------------------------------###
# ****** Iacobucci dataset (GSE13204) ****** #
###  Iacobucci dataset ----------
# Microarray gene expression profiling of 144 adult B-ALL patients and B-ALL negative for known molecular rearrangements.

# >>> Input data
# B-ALL Gene Expression data

load("Iacobucci2012.RData", verbose = TRUE)
ball_norm_rna_Iacobucci <- dataset@assayData$exprs
ball_metadata_Iacobucci <- dataset@phenoData@data

# Translate probe ID to gene names:
geneNames_ord <- geneNames[rownames(ball_norm_rna_Iacobucci)]
rownames(ball_norm_rna_Iacobucci) <- geneNames_ord

ball_metadata_Iacobucci$mutation <- recode_factor(ball_metadata_Iacobucci$GeneticBackground,
                                               "WildType"="Others",
                                               "BcrAbl"="BCR-ABL1/ph-like")

# >>> Exploratory analysis
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
gsvaPar <- gsvaParam(ball_norm_rna_Iacobucci, gs_symbol, kcdf="Gaussian")
gsva.BALL_Iacobucci <- gsva(gsvaPar, verbose=FALSE)


# >>> Save results
save(ball_norm_rna_Iacobucci,ball_metadata_Iacobucci,gsva.BALL_Iacobucci,file="BALL_Iacobucci.RData")
###-----------------------------------------------------------------------------------###



# Figure 5 - Panel A -----------------------------------------------------------------
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

pdf(paste(format(Sys.time(), "%Y%m%d"),"_pca_BALL_Li.pdf",sep=""), width=6, height = 6)
print(p_subgroup)
print(p_subgroup+theme(legend.position = "none"))
dev.off()
###-----------------------------------------------------------------------------------###


# Figure 5 - Panel B -----------------------------------------------------------------

# >>> Input data
# B-ALL Li Gene Expression data
load("BALL_Li.RData", verbose = TRUE)

# >>> Ploting
# Plot colors
ball_colors <- jdb_palette("corona", length(levels(ball_metadata_Li$mutation)))
names(ball_colors) <- levels(ball_metadata_Li$mutation)

regulons_color <- c("1"="#92c5de","2"="#2166ac","3"="#91cf60","4"="#dd1c77","5"="black")
age_color <- c("Paediatric"="#DAA520","Adult"="#2F4F4F")


# Plot  
# Get the regulons classification
bcell_biomarkers <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/differential_expression_analysis/bcell_biomarkers.rds"))

regulons_group <- readRDS(paste0(data_wd,"/osfstorage-archive/04_gene_regulatory_networks/regulons/TF_regulons_clusters.rds"))
regulons_group_plot <- c(rep("1",length(regulons_group[[1]])),rep("2",length(regulons_group[[2]])),rep("3",length(regulons_group[[3]])),rep("4",length(regulons_group[[4]])),rep("5",length(regulons_group[[5]])))
names(regulons_group_plot) <- c(regulons_group[[1]],regulons_group[[2]],regulons_group[[3]],regulons_group[[4]],regulons_group[[5]])

# Define heatmap color
cols <- jdb_palette("brewer_spectra")
col_fun <- colorRamp2(seq(-2, 2, length.out = length(cols)), cols)

# Heatmap
h_scaled_BALL_Li <- Heatmap(t(scale(t(gsva.BALL_Li))),
                            top_annotation = HeatmapAnnotation("B-ALL" = ball_metadata_Li$mutation,
                                                               "Age" = ball_metadata_Li$age,
                                                               col=list("B-ALL"= ball_colors,
                                                                        "Age" = age_color)),
                            left_annotation = rowAnnotation("TF_subgroup" = regulons_group_plot[rownames(gsva.BALL_Li)],
                                                            col=list("TF_subgroup" = regulons_color)),
                            col = col_fun,
                            column_split = ball_metadata_Li$mutation,
                            row_split = factor(regulons_group_plot[rownames(gsva.BALL_Li)],levels=c("1","2","3","4","5")),
                            column_gap = unit(0.5, "mm"),
                            row_gap = unit(0.5, "mm"),
                            show_column_names = FALSE,
                            show_row_names = FALSE,
                            row_title = NULL,
                            column_title = NULL,
                            row_names_gp = gpar(fontsize = 6),
                            name="GSVA score (scaled)"
)

h_scaled_BALL_Li_clust <- Heatmap(t(scale(t(gsva.BALL_Li))),
                                  top_annotation = HeatmapAnnotation("B-ALL" = ball_metadata_Li$mutation,
                                                                     "Age" = ball_metadata_Li$age,
                                                                     col=list("B-ALL"= ball_colors,
                                                                              "Age" = age_color)),
                                  left_annotation = rowAnnotation("TF_subgroup" = regulons_group_plot[rownames(gsva.BALL_Li)],
                                                                  col=list("TF_subgroup" = regulons_color)),
                                  col = col_fun,
                                  #column_split = group,
                                  row_split = factor(regulons_group_plot[rownames(gsva.BALL_Li)],levels=c("1","2","3","4","5")),
                                  column_gap = unit(0.5, "mm"),
                                  row_gap = unit(0.5, "mm"),
                                  show_column_names = FALSE,
                                  show_row_names = FALSE,
                                  row_title = NULL,
                                  column_title = NULL,
                                  row_names_gp = gpar(fontsize = 6),
                                  name="GSVA score (scaled)"
)

pdf(paste0(format(Sys.time(), "%Y%m%d"),"_heatmap_GSVA_regulons_BALL_Li.pdf"), width = 15, heigh=5)
h_scaled_BALL_Li
h_scaled_BALL_Li_clust
dev.off()
###-----------------------------------------------------------------------------------###


# Figure 5 - Panel C -----------------------------------------------------------------
# Over Representation Analysis (ORA) of B-ALL subtypes signatures over the B-cell subpopulations signatures. 

# >>> Required packages
library(dplyr)
library(viridis)
library(enrichplot)
library(clusterProfiler)
library(ggplot2)

# >>> Required auxiliar functions
source(paste0(data_wd,"/auxiliar_functions.R"))

# >>> Input data
# Human B-cell subpopulations biomarkers
bcell_biomarkers <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/differential_expression_analysis/Bcell_biomarkers.rds"))

# Human B-ALL subtype biomarkers
ball_biomarkers <- readRDS(paste0(data_wd,"/osfstorage-archive/05_B-ALL/ALL_subgroup_biomarkers.rds"))

# >>> Computing & Plotting
#For each mutation, is there an enrichment of a cell types signatures?
a0 <- ORA(list_query = ball_biomarkers, list_terms = bcell_biomarkers, x_lab="Leukemia subtype", y_lab="B-cell subpopulation")

pdf(paste0(format(Sys.time(), "%Y%m%d"),"_ORA_ALL_subtypes_cellType.pdf"), width = 6)
a0[[2]]
dev.off()
###-----------------------------------------------------------------------------------###


# Figure 5 - Panel D -----------------------------------------------------------------
# ORA for RUNX1-ETV6 high affinity binding sites (n=711) over the consensus OCRs 
# for each B-cell subpopulation.

# >>> Required packages
library(GSVA)
library(ComplexHeatmap)
library(BuenColors)
library(GenomicRanges)

# >>> Required auxiliar functions
source("auxiliar_functions.R")

# >>> Input data
# Human B-cell consensus peak per cell type
cell_type <- c("HSC","CLP","proB","preB","ImmatureB","Transitional_B","Naive_CD5pos","Naive_CD5neg")
cell_type_gr <- vector(mode="list", length=length(cell_type))
for(i in 1:length(cell_type)){
  dd <- read.table(paste0(data_wd,"/osfstorage-archive/01_atac_seq_data/consensus_OCRs/consensus_",cell_type[i],".bed"), header=FALSE, sep="\t")
  colnames(dd)[1:3] <- c("Chr","Start","End")
  cell_type_gr[[i]] <- makeGRangesFromDataFrame(dd)
}
names(cell_type_gr) <- cell_type
# Human B-cell consensus peaks
atac <- readRDS(paste0(data_wd,"/osfstorage-archive/01_atac_seq_data/atac_norm_data.rds"))
# Read ChiP ETV6-RUNX1 public data
chip_etv6_runx1 <- read.table(paste0(data_wd,"/osfstorage-archive/05_B-ALL/ETV6_RUNX1_ChiP.txt"), sep="\t", header=TRUE)
chip_etv6_runx1 <- matrix(unlist(strsplit(unique(chip_etv6_runx1$Peak),"_")),ncol=3, byrow = TRUE)
chip_etv6_runx1 <- data.frame(Chr=chip_etv6_runx1[,1],
                              Start=as.numeric(chip_etv6_runx1[,2]),
                              End=as.numeric(chip_etv6_runx1[,3]))
chip_etv6_runx1_gr <- makeGRangesFromDataFrame(chip_etv6_runx1)

bcell_colors <- readRDS(paste0(data_wd,"/osfstorage-archive/00_plot_parameters/bcell_colors.rds"))

# >>> Computing
enrich_output <- data.frame(matrix(nrow = 8, ncol = 10)) 
colnames(enrich_output) <- c("cell_type", "query_TF", "pval", "est",
                             "f11", "f01", "f10", "f00", "Lquery", "Lsubj")

for(i in 1:length(cell_type_gr)){
  cat(">>>>>",names(cell_type_gr)[i],"\n")
  output <- my_enrichPeakOverlap(chip_etv6_runx1_gr,cell_type_gr[[i]],nrow(atac))
  enrich_output$cell_type[i] <- names(cell_type_gr)[i]
  enrich_output$query_TF[i] <- "RUNX1-ETV6_ChiP"
  enrich_output$pval[i] <- output["pval"]
  enrich_output$est[i] <- output["estimate.odds ratio"]
  enrich_output$f11[i] <- output["f11"]
  enrich_output$f01[i] <- output["f01"]
  enrich_output$f10[i] <- output["f10"]
  enrich_output$f00[i] <- output["f00"]
  enrich_output$Lquery[i] <- output["Lquery"]
  enrich_output$Lsubj[i] <- output["Lsubj"]
}

data_to_plot <- enrich_output
data_to_plot$log10_pval <- log10(data_to_plot$pval)*-1
data_to_plot$cell_type <- factor(data_to_plot$cell_type, levels=cell_type)

# >>> Plotting
pdf(paste0(format(Sys.time(), "%Y%m%d"),"_ORA_ALL_subtypes_cellType_ATAC_ETV6RUNX1.pdf"), width = 6)
ggplot(data_to_plot, aes(cell_type, query_TF, colour= est, size=log10_pval)) + 
  geom_point() +
  scale_colour_gradientn(colors=jdb_palette("solar_extra")) +
  ylab("OCRs") + 
  xlab("Cell Type") +
  ggtitle("-log10(pval) overlap ETV6-RUNX1 ChiP-Seq by OCRs per cell type") +
  theme_minimal()
dev.off()
###-----------------------------------------------------------------------------------###




# Figure 5 - Panel E -----------------------------------------------------------------
# ORA of B-ALL subtypes signatures over the TF regulons defined by bulk approach (n=169).

# >>> Required packages
library(dplyr)
library(viridis)
library(enrichplot)
library(clusterProfiler)
library(ggplot2)

# >>> Required auxiliar functions
source("auxiliar_functions.R")

# >>> Input data
# Human B-cell regulons
# GRN data
tf_ocr_gene <- read.table(paste0(data_wd,"/osfstorage-archive/04_gene_regulatory_networks/tf_ocr_gene_interactions_curated.txt"),
                          sep="\t",
                          dec=",",
                          header=TRUE)

# Human B-ALL subtype biomarkers
ball_biomarkers <- readRDS(paste0(data_wd,"/osfstorage-archive/05_B-ALL/ALL_subgroup_biomarkers.rds"))
names(ball_biomarkers) <- recode(names(ball_biomarkers), 
                                 G1="MEF2D fusions",
                                 G2="TCF3-PBX1",
                                 G3="ETV6-RUNX1/like",
                                 G4="DUX4 fusions",
                                 G5="ZNF384 fusions",
                                 G6="BCR-ABL1/ph-like",
                                 G7="Hyperdiploidy",
                                 G8="KTM2A fusions")

# >>> Computing & Plotting
#There are an enrichment of the regulons within each ALL subtype?
regulons <- unique(tf_ocr_gene$tf_symbol)

list_regulons <- vector(mode="list", length=length(regulons))
names(list_regulons) <- regulons

for(i in regulons){
  list_regulons[[i]] <- unique(tf_ocr_gene[tf_ocr_gene$tf_symbol==i,"gene_ensembl"])
}

regulons_group <- readRDS(paste0(data_wd,"/osfstorage-archive/04_gene_regulatory_networks/regulons/TF_regulons_clusters.rds"))
regulons_group_plot <- c(rep("1",length(regulons_group[[1]])),rep("2",length(regulons_group[[2]])),rep("3",length(regulons_group[[3]])),rep("4",length(regulons_group[[4]])),rep("5",length(regulons_group[[5]])))
names(regulons_group_plot) <- c(regulons_group[[1]],regulons_group[[2]],regulons_group[[3]],regulons_group[[4]],regulons_group[[5]])

regulons_group_plot_ord <- c(rev(regulons_group[[1]]),
                             rev(regulons_group[[2]]),
                             rev(regulons_group[[3]]),
                             rev(regulons_group[[4]]),
                             rev(regulons_group[[5]]))


a1 <- ORA2_rev2(list_query = ball_biomarkers, 
                list_terms = list_regulons, 
                x_lab="Leukemia subtype", 
                TF_order=regulons_group_plot_ord,
                y_lab="Regulons",
                y_elements_to_highlight=c("MEF2D","PBX1","RUNX1","ETV6","TCF3"),
                col_rect_highlight="gray",
                order="both")

pdf(paste0(format(Sys.time(), "%Y%m%d"),"_ORA_regulons_ALL_subtypes.pdf"),height = 17, width=5)
print(a1[[2]])
print(a1[[3]])
dev.off()
###-----------------------------------------------------------------------------------###


# Figure 5 - Panel F-I -----------------------------------------------------------------
# >>> Required packages
library(ggplot2)
library(BuenColors)
library(tidyverse)

# >>> Input data
# B-ALL Gene Expression public datasets
load("BALL_Li.RData", verbose = TRUE)
load("BALL_Rainer.RData", verbose = TRUE)
load("BALL_Iacobucci.RData", verbose = TRUE)

# >>> Ploting
# Plot colors
ball_colors <- jdb_palette("corona", length(levels(ball_metadata_Li$mutation)))
names(ball_colors) <- levels(ball_metadata_Li$mutation)

# Plot  
# Get the regulons classification
regulons_group <- readRDS(paste0(data_wd,"/osfstorage-archive/04_gene_regulatory_networks/regulons/regulons_clusters.rds"))
regulons_group_plot <- c(rep("1",length(regulons_group[[1]])),rep("2",length(regulons_group[[2]])),rep("3",length(regulons_group[[3]])),rep("4",length(regulons_group[[4]])),rep("5",length(regulons_group[[5]])))
names(regulons_group_plot) <- c(regulons_group[[1]],regulons_group[[2]],regulons_group[[3]],regulons_group[[4]],regulons_group[[5]])

# Boxplots

study <- c("Li","Rainer","Iacobucci")
gsva_data_list <- list("Li"=gsva.BALL_Li,"Rainer"=gsva.BALL_Rainer,"Iacobucci"=gsva.BALL_Iacobucci)
metadata_list <- list("Li"=ball_metadata_Li,"Rainer"=ball_metadata_Rainer,"Iacobucci"=ball_metadata_Iacobucci)
target_gsva <- c("TCF3","MEF2D","FOXD3","SPI1","STAT5B","FOSB","OTX1","IRX2")

pdf(paste0(format(Sys.time(), "%Y%m%d"),"_boxplots_GSVA_regulons_BALL_paper.pdf"), width=10, height=5)
for(i in study){
  dd <- gsva_data_list[[i]][rownames(gsva_data_list[[i]])%in%target_gsva,]
  data_to_plot <- data.frame(melt(dd),ALL_subgroup=rep(metadata_list[[i]]$mutation, each=nrow(dd)))
  ball_colors <- jdb_palette("corona", length(levels(ball_metadata_Li$mutation)))
  names(ball_colors) <- levels(ball_metadata_Li$mutation)
  
  if(i=="Iacobucci"){
    ball_colors <- c("BCR-ABL1/ph-like"="#8c564b","Others"="#ADD8E6")}
  
  p <- ggplot(data_to_plot,aes(x=ALL_subgroup,y=value, fill=ALL_subgroup)) +
    geom_boxplot(outlier.size = 1) +
    facet_wrap(~Var1, ncol=4) +
    scale_fill_manual(values=ball_colors) +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1)) 
  print(p)
}
dev.off()
###-----------------------------------------------------------------------------------###








