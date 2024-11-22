#***************************************************
## Paper: "Uncovering the regulatory landscape of early human B-cell lymphopoiesis and its 
## implications in the pathogenesis of B-cell acute lymphoblastic leukemia"
## Authors: Planell N et al.
## Date: 2024
## Code: Figure 5
## Input data available at: https://osf.io/gswpy/?view_only=39cfd50d5a7e44e5a2d141869dcb1315
#***************************************************

# Define the directory were input data from OSF has been downloaded
data_wd <- setwd("../GitHub_code/")

# ****** Figure 5 ****** #

# Figure 5 - Panel A -----------------------------------------------------------------
# Schematic representation of B-ALL public data explored and PCA representation

# >>> Required packages
library(ggplot2)
library(BuenColors)
library(tidyverse)

# >>> Input data
# B-ALL Gene Expression data
ball_rna <- readRDS(paste0(data_wd,"/osfstorage-archive/05_B-ALL/ball_rna_norm_data.rds"))
ball_metadata <- readRDS(paste0(data_wd,"/osfstorage-archive/05_B-ALL/ball_metadata.rds"))

ball_metadata$subgroup <- recode_factor(ball_metadata$subgroup, 
                                        G1="MEF2D fusions",
                                        G2="TCF3-PBX1",
                                        G3="ETV6-RUNX1/like",
                                        G4="DUX4 fusions",
                                        G5="ZNF384 fusions",
                                        G6="BCR-ABL1/ph-like",
                                        G7="Hyperdiploidy",
                                        G8="KTM2A fusions")

# >>> Computing
# Multidimensional scaling plots (MDS) of Spearman's correlation matrix
res.cor <- 1-cor(ball_rna, method = "spearman")
fit <- cmdscale(res.cor,eig=TRUE, k=2) # k is the number of dim

data_to_plot <- data.frame(PC1=fit$points[,1],PC2=fit$points[,2], ALL_group=ball_metadata$subgroup)

# >>> Plotting
# Plot colors
ball_colors <- jdb_palette("corona", length(levels(ball_metadata$subgroup)))
names(ball_colors) <- levels(ball_metadata$subgroup)

# Plot
p <- ggplot(data_to_plot, aes(x = PC1, y = PC2, color = ALL_group)) +
  geom_point(size=3) +
  scale_color_manual(values=ball_colors) +
  labs(title="ALL subgroups", x = "Coordinate 1", y = "Coordinate 2") +
  theme_classic()

pdf(paste(format(Sys.time(), "%Y%m%d"),"_pca_BALL.pdf",sep=""), width=6, height = 6)
  print(p)
  print(p+theme(legend.position = "none"))
dev.off()
# ------------------------------------------------------------------------------------ #



## Figure 5 - 2023/07/30

library(ggplot2)
library(BuenColors)
library(circlize)
library(viridis)
library(RColorBrewer)
library(enrichplot)
library(clusterProfiler)
library(tidyverse)

setwd("C:/Users/nplanellpic/Documents/npp/project/ongoing/bcell/")


# Figure 5 - Panel B -----------------------------------------------------------------
# Over Representation Analysis (ORA) of B-ALL subtypes signatures over the B-cell subpopulations signatures. 

# >>> Required packages
library(dplyr)
library(viridis)
library(enrichplot)
library(clusterProfiler)
library(ggplot2)

# >>> Required auxiliar functions
source("auxiliar_functions.R")

# >>> Input data
# Human B-cell subpopulations biomarkers
bcell_biomarkers <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/differential_expression_analysis/Bcell_biomarkers.rds"))

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
#For each mutation, is there an enrichment of a cell types signatures?
a0 <- ORA(list_query = ball_biomarkers, list_terms = bcell_biomarkers, x_lab="Leukemia subtype", y_lab="B-cell subpopulation")

pdf(paste0(format(Sys.time(), "%Y%m%d"),"_ORA_ALL_subtypes_cellType.pdf"), width = 6)
  a0[[2]]
dev.off()
###-----------------------------------------------------------------------------------###



# Figure 5 - Panel C -----------------------------------------------------------------
# Heatmap representation of the gene set activity score (GSVA) computed for 
# each human B-cell subpopulation signature (n=8; in rows) in each B-ALL sample (n=356; in columns). 

# >>> Required packages
library(GSVA)
library(ComplexHeatmap)
library(BuenColors)

# >>> Input data
# Human B-cell subpopulations biomarkers
bcell_biomarkers <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/differential_expression_analysis/bcell_biomarkers.rds"))
# B-ALL Gene Expression data
ball_rna <- readRDS(paste0(data_wd,"/osfstorage-archive/05_B-ALL/ball_rna_norm_data.rds"))
ball_metadata <- readRDS(paste0(data_wd,"/osfstorage-archive/05_B-ALL/ball_metadata.rds"))

# >>> Computing
# GSVA score of target genes
gs <- bcell_biomarkers
gsvaPar <- gsvaParam(ball_rna, gs)
gsva.es <- gsva(gsvaPar, verbose=FALSE)

ball_metadata$subgroup <- recode_factor(ball_metadata$subgroup, 
                                   G1="MEF2D fusions",
                                   G2="TCF3-PBX1",
                                   G3="ETV6-RUNX1/like",
                                   G4="DUX4 fusions",
                                   G5="ZNF384 fusions",
                                   G6="BCR-ABL1/ph-like",
                                   G7="Hyperdiploidy",
                                   G8="KTM2A fusions")

data_to_plot <- data.frame(melt(gsva.es),ALL_subgroup=rep(ball_metadata$subgroup, each=nrow(gsva.es)))

# >>> Plotting
# Plot colors
bcell_colors <- readRDS(paste0(data_wd,"/osfstorage-archive/00_plot_parameters/bcell_colors.rds"))
ball_colors <- jdb_palette("corona", length(levels(ball_metadata$subgroup)))
names(ball_colors) <- levels(ball_metadata$subgroup)

# Plot  
h <- Heatmap(scale(gsva.es),
             top_annotation = HeatmapAnnotation("B-ALL" = ball_metadata$subgroup, 
                                                col = list("B-ALL" = ball_colors),
                                                show_annotation_name = FALSE),
             left_annotation = rowAnnotation(human = factor(names(bcell_colors), levels=names(bcell_colors)), 
                                             col = list(human = bcell_colors),
                                             show_annotation_name = FALSE),
             col = jdb_palette("brewer_purple"),
             column_order = order(ball_metadata$subgroup),
             column_split = ball_metadata$subgroup,
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


pdf(paste(format(Sys.time(), "%Y%m%d"),"_heatmap_gsva_Bcell_biomarkers_in_BALL.pdf",sep=""))
  draw(h)
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
ggplot(data_to_plot, aes(cell_type, query_TF, colour= est, size=log10_pval)) + 
  geom_point() +
  scale_colour_gradientn(colors=jdb_palette("solar_extra")) +
  ylab("TF OCRs") + 
  xlab("Cell Type") +
  ggtitle("-log10(pval) overlap ETV6-RUNX1 ChiP-Seq by OCRs per cell type") +
  theme_minimal()
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

a1 <- ORA2(list_query = ball_biomarkers, 
           list_terms = list_regulons, 
           x_lab="Leukemia subtype", 
           y_lab="Regulons",
           y_elements_to_highlight=c("MEF2D","PBX1","RUNX1","ETV6","TCF3"),
           col_rect_highlight="gray",
           order="both")

png(paste0(format(Sys.time(), "%Y%m%d"),"_ORA_regulons_ALL_subtypes.png"),height = 1300, width=350)
  a1[[2]]
  a1[[3]]
dev.off()

###-----------------------------------------------------------------------------------###



# Figure 5 - Panel F -----------------------------------------------------------------
# ORA of B-ALL subtypes signatures over the single-cell B-cell subpopulations signatures.  

# >>> Required packages
library(dplyr)
library(viridis)
library(enrichplot)
library(clusterProfiler)
library(ggplot2)
library(SeqGSEA)

# >>> Required auxiliar functions
source("auxiliar_functions.R")

# >>> Input data
# Human B-cell single-cell subpopulations biomarkers
scbcell_biomarkers <- readRDS(paste0(data_wd,"/osfstorage-archive/07_single-cell/sc_rnaseq_markers_per_cell_type.rds"))
sc_cell_type <- c("HSC","MPP","LMPP","CLP","cycling_proB","proB","cycling_preBI","preBI")
list_scbcell_biomarkers <- vector(mode="list", length=length(sc_cell_type))
names(list_scbcell_biomarkers) <- sc_cell_type

for(i in sc_cell_type){
  gg <- unique(scbcell_biomarkers[scbcell_biomarkers$cluster==i,"gene"])
  #Translate to ENSEMBL ID
  gg_ensmbl_id <-   convertSymbol2Ensembl(gg)
  list_scbcell_biomarkers[[i]] <- unique(gg_ensmbl_id$ensembl_gene_id)
}

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
#For each mutation, is there an enrichment of a cell types signatures?
a0 <- ORA(list_query = ball_biomarkers, list_terms = list_scbcell_biomarkers, x_lab="Leukemia subtype", y_lab="B-cell subpopulations (Single-cell)")

pdf(paste0(format(Sys.time(), "%Y%m%d"),"_ORA_ALL_subtypes_cellType_singleCell.pdf"), width = 6)
  a0[[2]]
dev.off()

###-----------------------------------------------------------------------------------###


# Figure 5 - Panel G -----------------------------------------------------------------
# Disease exploration - ORA B-ALL subtypes in single-cell regulons

# >>> Required packages
library(dplyr)
library(viridis)
library(enrichplot)
library(clusterProfiler)
library(ggplot2)
library(SeqGSEA)

# >>> Required auxiliar functions
source("auxiliar_functions.R")

# >>> Input data
# Human B-cell regulons
# GRN data single-cell
sc_regulons <- read.csv(paste0(data_wd,"/osfstorage-archive/07_single-cell/eRegulon_metadata_filtered.csv"))
regulons <- unique(sc_regulons$Region_signature_name)
list_regulons <- vector(mode="list", length=length(regulons))
names(list_regulons) <- regulons

for(i in regulons){
  gg <- unique(sc_regulons[sc_regulons$Region_signature_name==i,"Gene"])
  #Translate to ENSEMBL ID
  gg_ensmbl_id <-   convertSymbol2Ensembl(gg)
  list_regulons[[i]] <- unique(gg_ensmbl_id$ensembl_gene_id)
}

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

# >>> Computing
#There are an enrichment of the regulons within each ALL subtype?
# Significant regulons
regulons_names <- c("TCF4_+_(1365r)","RUNX1_+_(6147r)","PAX5_+_(854r)","LEF1_+_(39r)","EBF1_+_(1610r)",
                    "ZNF367_extended_+_(252r)","MYBL2_+_(8r)","CEBPA_+_(54r)","NFE2_+_(1004r)","THRB_+_(951r)",
                    "RUNX2_+_(4210r)","PLAGL1_extended_+_(653r)","PBX1_extended_+_(316r)","NFIA_+_(438r)","MEIS1_+_(2276r)",
                    "MECOM_+_(153r)","LYL1_+_(335r)","HMGA2_+_(567r)","HLF_+_(162r)","FOS_+_(2445r)",
                    "FOSB_+_(165r)","RUNX1_-_(38r)","PAX5_-_(167r)","NR3C1_-_(59r)","LEF1_-_(1210r)","EBF1_-_(732r)","BACH2_-_(2433r)",
                    "ZNF385D_-_(190r)","THRB_-_(73r)","PRDM16_extended_-_(37r)","NFIA_-_(65r)",
                    "MEIS1_-_(297r)","KLF9_-_(81r)","HMGA2_-_(70r)")

list_regulons_ordered <- list_regulons[regulons_names]

a1 <- ORA2(list_query = ball_biomarkers, 
           list_terms = list_regulons_ordered, 
           x_lab="Leukemia subtype", 
           y_lab="Regulons",
           y_elements_to_highlight=c(grep("EBF1",regulons_names,value=TRUE),grep("PAX5",regulons_names,value=TRUE),grep("ETV6",regulons_names,value=TRUE),grep("RUNX2",regulons_names,value=TRUE),grep("MEIS1",regulons_names,value=TRUE)),
           col_rect_highlight="gray",
           order="both")

pdf(paste0(format(Sys.time(), "%Y%m%d"),"_ORA_significant_sc_regulons_ALL_subtypes.pdf"), height=10, width=5.5)
  a1[[2]]
  a1[[3]]
dev.off()

###-----------------------------------------------------------------------------------###


# Figure 5 - Panel H -----------------------------------------------------------------
# Box plot showing the distribution of ELK3 gene expression by B-ALL subtypes. 

# >>> Required packages
library(ggplot2)
library(BuenColors)
library(tidyverse)
library(ggsignif)

# >>> Input data
# Gene Expression data
rna_anno <- readRDS(paste0(data_wd,"/osfstorage-archive/01_rna_seq_data/rna_gene_annotation.rds"))

# B-ALL Gene Expression data
ball_rna <- readRDS(paste0(data_wd,"/osfstorage-archive/05_B-ALL/ball_rna_norm_data.rds"))
ball_metadata <- readRDS(paste0(data_wd,"/osfstorage-archive/05_B-ALL/ball_metadata.rds"))

ball_metadata$subgroup <- recode_factor(ball_metadata$subgroup, 
                                        G1="MEF2D fusions",
                                        G2="TCF3-PBX1",
                                        G3="ETV6-RUNX1/like",
                                        G4="DUX4 fusions",
                                        G5="ZNF384 fusions",
                                        G6="BCR-ABL1/ph-like",
                                        G7="Hyperdiploidy",
                                        G8="KTM2A fusions")

# >>> Plotting
data_to_plot <- data.frame(value=ball_rna[rna_anno$ensembl_gene_id[rna_anno$hgnc_symbol=="ELK3"],],
                           gene=rep("ELK3",nrow(ball_metadata)),
                           ALL_subgroup=ball_metadata$subgroup)

p <- ggplot(data_to_plot, aes(x=ALL_subgroup, y=value, fill=ALL_subgroup))+
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values=jdb_palette("corona")) +
  geom_signif(comparisons = list(c("TCF3-PBX1","ETV6-RUNX1/like"),
                                 c("MEF2D fusions","ETV6-RUNX1/like"),
                                 c("ETV6-RUNX1/like","DUX4 fusions"),
                                 c("ETV6-RUNX1/like","ZNF384 fusions"),
                                 c("ETV6-RUNX1/like","BCR-ABL1/ph-like"),
                                 c("ETV6-RUNX1/like","Hyperdiploidy"),
                                 c("ETV6-RUNX1/like","KTM2A fusions")), 
              map_signif_level=T,
              textsize=3,
              size=0.5,
              tip_length = 0,
              step_increase = 0.1) +
  facet_wrap(~gene)+
  ylab("Gene expression") +
  xlab("") +
  theme(axis.text.x = element_text(angle=45, hjust=1)) 


pdf(paste(format(Sys.Date(),"%Y%m%d"),"_ELK3_expression_in_ALL.pdf",sep=""), height = 5)
  print(p)
dev.off()
###-----------------------------------------------------------------------------------###



# Figure 5 - Panel I -----------------------------------------------------------------
# Box plot showing the distribution of the gene set activity score (GSVA) computed on the ELK3-regulon-repressed gene set by B-ALL subtypes.

# >>> Required packages
library(ggplot2)
library(BuenColors)
library(tidyverse)
library(ggsignif)
library(GSVA)

# >>> Input data
# B-ALL Gene Expression data
ball_rna <- readRDS(paste0(data_wd,"/osfstorage-archive/05_B-ALL/ball_rna_norm_data.rds"))
ball_metadata <- readRDS(paste0(data_wd,"/osfstorage-archive/05_B-ALL/ball_metadata.rds"))

ball_metadata$subgroup <- recode_factor(ball_metadata$subgroup, 
                                        G1="MEF2D fusions",
                                        G2="TCF3-PBX1",
                                        G3="ETV6-RUNX1/like",
                                        G4="DUX4 fusions",
                                        G5="ZNF384 fusions",
                                        G6="BCR-ABL1/ph-like",
                                        G7="Hyperdiploidy",
                                        G8="KTM2A fusions")
# GRN data
tf_ocr_gene <- read.table(paste0(data_wd,"/osfstorage-archive/04_gene_regulatory_networks/tf_ocr_gene_interactions_curated.txt"),
                          sep="\t",
                          dec=",",
                          header=TRUE)

# >>> Computing
target_gene <- "ELK3"
elk3_regulon <- unique(tf_ocr_gene[tf_ocr_gene$tf_symbol==target_gene,"gene_ensembl"])

list_target_genes_ELK3_pos <- unique(tf_ocr_gene[which(tf_ocr_gene$tf_symbol==target_gene & tf_ocr_gene$gene_tf_rho>0),"gene_ensembl"])
list_target_genes_ELK3_neg <- unique(tf_ocr_gene[which(tf_ocr_gene$tf_symbol==target_gene & tf_ocr_gene$gene_tf_rho<0),"gene_ensembl"])

# GSVA score of target genes
gs <- list(elk3=list_target_genes_ELK3_pos)
gsvaPar <- gsvaParam(ball_rna, gs)
gsva.es_pos <- gsva(gsvaPar, verbose=FALSE)

gs <- list(elk3=list_target_genes_ELK3_neg)
gsvaPar <- gsvaParam(ball_rna, gs)
gsva.es_neg <- gsva(gsvaPar, verbose=FALSE)

data_to_plot <- data.frame(value_pos=gsva.es_pos[1,],
                           value_neg=gsva.es_neg[1,],
                           gene=rep("ELK3 regulon",nrow(ball_metadata)),
                           ALL_subgroup=ball_metadata$subgroup)

# >>> Plotting
p1 <- ggplot(data_to_plot, aes(x=ALL_subgroup, y=value_pos, fill=ALL_subgroup))+
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values = jdb_palette("corona")) +
  geom_signif(comparisons =list(c("TCF3-PBX1","ETV6-RUNX1/like"),
                                c("MEF2D fusions","ETV6-RUNX1/like"),
                                c("ETV6-RUNX1/like","DUX4 fusions"),
                                c("ETV6-RUNX1/like","ZNF384 fusions"),
                                c("ETV6-RUNX1/like","BCR-ABL1/ph-like"),
                                c("ETV6-RUNX1/like","Hyperdiploidy"),
                                c("ETV6-RUNX1/like","KTM2A fusions")), 
              map_signif_level=T,
              textsize=3,
              size=0.5,
              tip_length = 0,
              step_increase = 0.1) +
  facet_wrap(~gene, nrow = 2)+
  ylab("GSVA score") +
  xlab("") +
  theme(axis.text.x = element_text(angle=45, hjust=1)) 

p2 <- ggplot(data_to_plot, aes(x=ALL_subgroup, y=value_neg, fill=ALL_subgroup))+
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values = jdb_palette("corona")) +
  geom_signif(comparisons =list(c("TCF3-PBX1","ETV6-RUNX1/like"),
                                c("MEF2D fusions","ETV6-RUNX1/like"),
                                c("ETV6-RUNX1/like","DUX4 fusions"),
                                c("ETV6-RUNX1/like","ZNF384 fusions"),
                                c("ETV6-RUNX1/like","BCR-ABL1/ph-like"),
                                c("ETV6-RUNX1/like","Hyperdiploidy"),
                                c("ETV6-RUNX1/like","KTM2A fusions")), 
              map_signif_level=T,
              textsize=3,
              size=0.5,
              tip_length = 0,
              step_increase = 0.1) +
  facet_wrap(~gene, nrow = 2)+
  ylab("GSVA score") +
  xlab("") +
  theme(axis.text.x = element_text(angle=45, hjust=1)) 

pdf(paste(format(Sys.time(), "%Y%m%d"),"_boxplot_BALL_ELK3_regulons_gsva_score.pdf",sep=""), height = 5)
  print(p1)
  print(p2)
dev.off()
###-----------------------------------------------------------------------------------###
