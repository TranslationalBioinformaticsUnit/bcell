#***************************************************
## Paper: "Uncovering the cis-regulatory program of early human B-cell commitment and its implications
## in the pathogenesis of acute lymphoblastic leukemia"
## Authors: Planell N et al.
## Date: 2023
## Code: Figure 4
## Input data available at: 
#***************************************************




#***************************************************
## Figure 4 A
#***************************************************
## Schematic diagram of samples
#***************************************************




#***************************************************
## Figure 4 B
#***************************************************
## ORA of BCP ALL subtypes biomarkers over the biomarkers of different B cell types

library(clusterProfiler)
library(ggplot2)
library(circlize)
library(viridis)

# Input data:
# Disease subgroup biomarkers: Supplementary Table 4 (ALL_subgroup_biomarkers) or "ALL_subgroup_biomarkers.rds"
  # available at OSF files: 01_rna_seq_data/ALL_subgroup_biomarkers.rds
# B cell subpopulation biomarkers: "Bcell_subpopulations_biomarkers.rds"
  # available at OSF files: 01_rna_seq_data/rna_norm_data.txt

# Auxiliar functions
source("auxiliar_functions.R")

ALL_subgroup_biomarkers <- readRDS("ALL_subgroup_biomarkers.rds")
Bcell_subpopulations_biomarkers <- readRDS("Bcell_subpopulations_biomarkers.rds")

  
#For each mutation we assessed if there is an enrichment of a cell types signatures.
output <- ORA(list_query = ALL_subgroup_biomarkers, list_terms = Bcell_subpopulations_biomarkers)

#***************************************************




#***************************************************
## Figure 4 C
#***************************************************
## ORA for BCP ALL subtypes G1 (MEF2D fusion), G2 (TCF3-PBX1), and G3 (ETV6-RUNX1) 
## over the set of genes regulated by MEFD2, TCF3 and PBX1 or ETV6 and RUNX1, respectively, 
## within each cell type.

library(clusterProfiler)
library(ggplot2)
library(circlize)
library(viridis)
library(igraph)

# Input data:
# Disease subgroup biomarkers: Supplementary Table 4 (ALL_subgroup_biomarkers) or "ALL_subgroup_biomarkers.rds"
  # available at OSF files: 01_rna_seq_data/ALL_subgroup_biomarkers.rds
# B cell subpopulation biomarkers: "Bcell_subpopulations_biomarkers.rds"
  # available at OSF files: 01_rna_seq_data/rna_norm_data.txt
# GRN for each cell type: "bcell_grn_by_cell_type.rds"
  # available at OSF files: 04_gene_regulatory_networks/bcell_grn_by_cell_type.rds
# Gene annotation: "rna_gene_annotation.txt"
  # available at OSF files: 01_rna_seq_data/rna_gene_annotation.txt


# Auxiliar functions
source("auxiliar_functions.R")

ALL_subgroup_biomarkers <- readRDS("ALL_subgroup_biomarkers.rds")
grns_by_cell_type <- readRDS("bcell_grn_by_cell_type.rds")
rna_annotation <- read.table("rna_gene_annotation.txt",header=TRUE, dec=".", sep="\t")

# MEF2D (G1)---------------------------------
cell_type_tf_mef2d <- lapply(grns_by_cell_type, FUN=function(x){
  aa <- as_edgelist(x, names = TRUE)
  return(unique(aa[aa[,1]==rna_annotation$ensembl_gene_id[rna_annotation$hgnc_symbol=="MEF2D"],2]))
})

#For MEF2D mutation, is there an enrichment of a cell types signatures?
ora_mef2d <- ORA(list_query = ALL_subgroup_biomarkers[1], list_terms = cell_type_tf_mef2d)


#PBX1 - TCF3 (G2)---------------------------------
cell_type_tf_pbx1_tcf3 <- lapply(grns_by_cell_type, FUN=function(x){
  aa <- as_edgelist(x, names = TRUE)
  return(unique(aa[aa[,1]%in%c(rna_annotation$ensembl_gene_id[rna_annotation$hgnc_symbol=="PBX1"], rna_annotation$ensembl_gene_id[rna_annotation$hgnc_symbol=="TCF3"]),2]))
})

#For PBX1-TCF3 mutation, is there an enrichment of a cell types signatures?
ora_pbx1_tcf3 <- ORA(list_query = ALL_subgroup_biomarkers[2], list_terms = cell_type_tf_pbx1_tcf3)
# -----------------------------------------------


#RUNX1 - ETV6 (G3)---------------------------------
cell_type_tf_runx1_etv6 <- lapply(grns_by_cell_type, FUN=function(x){
  aa <- as_edgelist(x, names = TRUE)
  return(unique(aa[aa[,1]%in%c(rna_annotation$ensembl_gene_id[rna_annotation$hgnc_symbol=="ETV6"], rna_annotation$ensembl_gene_id[rna_annotation$hgnc_symbol=="RUNX1"]),2]))
})

#For RUNX1-ETV6 mutation, is there an enrichment of a cell types signatures?
ora_runx1_etv6 <- ORA(list_query = ALL_subgroup_biomarkers[3], list_terms = cell_type_tf_runx1_etv6)

#***************************************************




#***************************************************
## Figure 4 D
#***************************************************
## ORA for RUNX1-ETV6 high affinity binding sites (n=711) over the consensus OCRs for each cell type.

library(clusterProfiler)
library(ggplot2)
library(circlize)
library(BuenColors)
library(igraph)
library(GenomicRanges)

# Input data:
# RUNX1-ETV6 high affinity binding sites (n=711): "ETV6_RUNX1_Fig4d.txt"
  # available at OSF files: 05_B-ALL/ETV6_RUNX1_Fig4d.txt
# Consensus OCRs per cell type: "consensus_XXX.bed" 
  # available at OSF files: 01_atac_seq_data/consensus_OCRs/consensus_XXX.bed
# B cell subpopulation colors: "bcell_colors.rds"


# Auxiliar functions
source("auxiliar_functions.R")

## ETV6-RUNX1 (union) vs 711 

# Read public data
chip_etv6_runx1 <- read.table("ETV6_RUNX1_Fig4d.txt", sep="\t", header=TRUE)
chip_etv6_runx1 <- matrix(unlist(strsplit(unique(chip_etv6_runx1$Peak),"_")),ncol=3, byrow = TRUE)
chip_etv6_runx1 <- data.frame(Chr=chip_etv6_runx1[,1],
                              Start=as.numeric(chip_etv6_runx1[,2]),
                              End=as.numeric(chip_etv6_runx1[,3]))
chip_etv6_runx1_gr <- makeGRangesFromDataFrame(chip_etv6_runx1) #711


#Consensus OCRs per cell type
bcell_colors <- readRDS("bcell_colors.rds")
cell_type <- names(bcell_colors)

cell_type_gr <- vector(mode="list", length=length(cell_type))

for(i in 1:length(cell_type)){
  dd <- read.table(paste0("consensus_",cell_type[i],".bed"), header=FALSE, sep="\t")
  colnames(dd)[1:3] <- c("Chr","Start","End")
  cell_type_gr[[i]] <- makeGRangesFromDataFrame(dd)
}
names(cell_type_gr) <- cell_type


#by cell type
enrich_output <- data.frame(matrix(nrow = 8, ncol = 10)) 
colnames(enrich_output) <- c("cell_type", "query_TF", "pval", "est",
                             "f11", "f01", "f10", "f00", "Lquery", "Lsubj")

for(i in 1:length(cell_type_gr)){
  cat(">>>>>",names(cell_type_gr)[i],"\n")
  output <- my_enrichPeakOverlap(chip_etv6_runx1_gr,cell_type_gr[[i]],nrow(ocr_tf))
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

#odds ratio
ggplot(data_to_plot, aes(cell_type, query_TF, colour= est, size=log10_pval)) + 
  geom_point() +
  scale_colour_gradientn(colors=jdb_palette("solar_extra")) +
  ylab("TF OCRs") + 
  xlab("Cell Type") +
  ggtitle("-log10(pval) overlap ETV6-RUNX1 ChiP-Seq by OCRs per cell type") +
  theme_minimal()
###-----------------------------------------------------------------------------------###




#***************************************************
## Figure 4 E
#***************************************************
## Violin plot of ELK3 gene expression by BCP ALL subtypes

library(ggplot2)

# Input data:
# B-ALL normalized rna data: "BALL_rna_norm_data.txt"
  # available at OSF files: 05_B-ALL/BALL_rna_norm_data.txt
# B-ALL metadata: "BALL_metadata.txt"
  # available at OSF files: 05_B-ALL/BALL_metadata.txt
# Gene annotation: "rna_gene_annotation.txt"
  # available at OSF files: 01_rna_seq_data/rna_gene_annotation.txt

# Auxiliar functions
source("auxiliar_functions.R")

rna_annotation <- read.table("rna_gene_annotation.txt",header=TRUE, dec=".", sep="\t")
BALL_metadata <- read.table("BALL_metadata.txt", sep="\t", dec=".", header=TRUE)
BALL_rna_norm_data <- read.table("BALL_rna_norm_data.txt", sep="\t", dec=".", header=TRUE)

data_to_plot <- data.frame(gene=as.numeric(BALL_rna_norm_data[rna_annotation$ensembl_gene_id[rna_annotation$hgnc_symbol=="ELK3"],]),group=BALL_metadata$subgroup)

print(ggplot(data_to_plot, aes(x=group, y=gene, fill=group)) + 
        geom_violin(trim=FALSE) +
        scale_fill_manual(values=jdb_palette("corona")) +
        stat_summary(fun.data=data_summary) +
        theme_classic())

###-----------------------------------------------------------------------------------###




#***************************************************
## Figure 4 F
#***************************************************
## Representation of ATAC-Seq signal for rs17481869 genomic region for the 8 cell subpopulations.

library(ggplot2)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

# Bigwig image
# Input data:
# Bigwig information: "bcell_bwplots.RData"
# Consensus OCRs: "consensus_OCRs_all_cell_subpopulations.bed"
  # available at OSF files: 01_atac_seq_data/consensus_OCRs/consensus_OCRs_all_cell_subpopulations.bed
# H3K4me1 regions: "H3K4me1_regions.bed"
  # available at OSF files: 02_cCREs_OCR_gene_links/public_epigenetic_data_for_validation/histone_ChiPseq/h3k4me1/H3K4me1_regions.bed
# H3K4me3 regions: "H3K4me3_regions.bed"
  # available at OSF files: 02_cCREs_OCR_gene_links/public_epigenetic_data_for_validation/histone_ChiPseq/h3k4me3/H3K4me3_regions.bed
# H3K27ac regions: "H3K27ac_regions.bed"
  # available at OSF files: 02_cCREs_OCR_gene_links/public_epigenetic_data_for_validation/histone_ChiPseq/h3k27ac/H3K27ac_regions.bed

# Auxiliar functions
source("auxiliar_functions.R")


load("bcell_bwplots.RData", verbose = TRUE)
consensus_OCRs_granges <- import("consensus_OCRs_all_cell_subpopulations.bed", format = "BED")
h3k4me1 <- import("H3K4me1_regions.bed", format = "BED")
h3k4me3 <- import("H3K4me3_regions.bed", format = "BED")
h3k27ac <- import("H3K27ac_regions.bed", format = "BED")

# Bigwig plot
#rs17481869 (chr2:145,365,198-145,371,227) 2:145366886
snp <- makeGRangesFromDataFrame(df=data.frame(chr="chr2",start=145366886,end=145366886))

gene <- "rs17481869"
bw_plot3(accDT_by_group, peaks=consensus_OCRs_granges, gene=gene, ylims=c(0,50), start=145365198, end=145371227, chr="chr2",
         h3k4me1, h3k4me3, h3k27ac, snp)

###-----------------------------------------------------------------------------------###




#***************************************************
## Figure 4 G
#***************************************************
## Representation of ATAC-Seq signal for the loci within the ELK3 gene (rs2268508, rs17736737, and rs4762284).

library(ggplot2)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

# Bigwig image
# Input data:
# Bigwig information: "bcell_bwplots.RData"
# Consensus OCRs: "consensus_OCRs_all_cell_subpopulations.bed"
  # available at OSF files: 01_atac_seq_data/consensus_OCRs/consensus_OCRs_all_cell_subpopulations.bed
# H3K4me1 regions: "H3K4me1_regions.bed"
  # available at OSF files: 02_cCREs_OCR_gene_links/public_epigenetic_data_for_validation/histone_ChiPseq/h3k4me1/H3K4me1_regions.bed
# H3K4me3 regions: "H3K4me3_regions.bed"
  # available at OSF files: 02_cCREs_OCR_gene_links/public_epigenetic_data_for_validation/histone_ChiPseq/h3k4me3/H3K4me3_regions.bed
# H3K27ac regions: "H3K27ac_regions.bed"
  # available at OSF files: 02_cCREs_OCR_gene_links/public_epigenetic_data_for_validation/histone_ChiPseq/h3k27ac/H3K27ac_regions.bed

# Auxiliar functions
source("auxiliar_functions.R")


load("bcell_bwplots.RData", verbose = TRUE)
consensus_OCRs_granges <- import("consensus_OCRs_all_cell_subpopulations.bed", format = "BED")
h3k4me1 <- import("H3K4me1_regions.bed", format = "BED")
h3k4me3 <- import("H3K4me3_regions.bed", format = "BED")
h3k27ac <- import("H3K27ac_regions.bed", format = "BED")

# Bigwig plot
#rs2268508 (chr12:96208965), rs17736737 (chr12:96215566), rs4762284 (12:96218984) (chr12:96,206,290-96,221,319)
snp <- makeGRangesFromDataFrame(df=data.frame(chr=c("chr12","chr12","chr12"),
                                              start=c(96208965, 96215566, 96218984),
                                              end=c(96208965, 96215566, 96218984)))

gene <- "rs2268508"
bw_plot3(accDT_by_group, peaks=consensus_OCRs_granges, gene=gene, ylims=c(0,600), start=96206290, end=96221319, chr="chr12",
         h3k4me1, h3k4me3, h3k27ac, snp)

###-----------------------------------------------------------------------------------###




#***************************************************
## Figure 4 H
#***************************************************
## Representation of ATAC-Seq signal for the rs2239633 locus, located in the CEBPE promoter region.

library(ggplot2)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

# Bigwig image
# Input data:
# Bigwig information: "bcell_bwplots.RData"
# Consensus OCRs: "consensus_OCRs_all_cell_subpopulations.bed"
  # available at OSF files: 01_atac_seq_data/consensus_OCRs/consensus_OCRs_all_cell_subpopulations.bed
# H3K4me1 regions: "H3K4me1_regions.bed"
  # available at OSF files: 02_cCREs_OCR_gene_links/public_epigenetic_data_for_validation/histone_ChiPseq/h3k4me1/H3K4me1_regions.bed
# H3K4me3 regions: "H3K4me3_regions.bed"
  # available at OSF files: 02_cCREs_OCR_gene_links/public_epigenetic_data_for_validation/histone_ChiPseq/h3k4me3/H3K4me3_regions.bed
# H3K27ac regions: "H3K27ac_regions.bed"
  # available at OSF files: 02_cCREs_OCR_gene_links/public_epigenetic_data_for_validation/histone_ChiPseq/h3k27ac/H3K27ac_regions.bed

# Auxiliar functions
source("auxiliar_functions.R")


load("bcell_bwplots.RData", verbose = TRUE)
consensus_OCRs_granges <- import("consensus_OCRs_all_cell_subpopulations.bed", format = "BED")
h3k4me1 <- import("H3K4me1_regions.bed", format = "BED")
h3k4me3 <- import("H3K4me3_regions.bed", format = "BED")
h3k27ac <- import("H3K27ac_regions.bed", format = "BED")

# Bigwig plot
#rs2239633 (14:23119848): chr14:23,118,522-23,122,139 / chr14:23,117,135-23,122,561
snp <- makeGRangesFromDataFrame(df=data.frame(chr="chr14",start=23119848,end=23119848))

gene <- "rs2239633" 
bw_plot3(accDT_by_group, peaks=consensus_OCRs_granges, gene=gene, ylims=c(0,200), start=23117135, end=23122561, chr="chr14",
         h3k4me1, h3k4me3, h3k27ac, snp)

###-----------------------------------------------------------------------------------###




#***************************************************
## Figure 4 I
#***************************************************
## Violin plot of CEBPE gene expression by BCP ALL subtypes

library(ggplot2)

# Input data:
# B-ALL normalized rna data: "BALL_rna_norm_data.txt"
  # available at OSF files: 05_B-ALL/BALL_rna_norm_data.txt
# B-ALL metadata: "BALL_metadata.txt"
  # available at OSF files: 05_B-ALL/BALL_metadata.txt
# Gene annotation: "rna_gene_annotation.txt"
  # available at OSF files: 01_rna_seq_data/rna_gene_annotation.txt

# Auxiliar functions
source("auxiliar_functions.R")

rna_annotation <- read.table("rna_gene_annotation.txt",header=TRUE, dec=".", sep="\t")
BALL_metadata <- read.table("BALL_metadata.txt", sep="\t", dec=".", header=TRUE)
BALL_rna_norm_data <- read.table("BALL_rna_norm_data.txt", sep="\t", dec=".", header=TRUE)

data_to_plot <- data.frame(gene=as.numeric(BALL_rna_norm_data[rna_annotation$ensembl_gene_id[rna_annotation$hgnc_symbol=="CEBPE"],]),group=BALL_metadata$subgroup)

print(ggplot(data_to_plot, aes(x=group, y=gene, fill=group)) + 
        geom_violin(trim=FALSE) +
        scale_fill_manual(values=jdb_palette("corona")) +
        stat_summary(fun.data=data_summary) +
        theme_classic())

###-----------------------------------------------------------------------------------###




#***************************************************
## Figure 4 J
#***************************************************
## Representation of ATAC-Seq signal for the rs35837782 locus.

library(ggplot2)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)

# Bigwig image
# Input data:
# Bigwig information: "bcell_bwplots.RData"
# Consensus OCRs: "consensus_OCRs_all_cell_subpopulations.bed"
  # available at OSF files: 01_atac_seq_data/consensus_OCRs/consensus_OCRs_all_cell_subpopulations.bed
# H3K4me1 regions: "H3K4me1_regions.bed"
  # available at OSF files: 02_cCREs_OCR_gene_links/public_epigenetic_data_for_validation/histone_ChiPseq/h3k4me1/H3K4me1_regions.bed
# H3K4me3 regions: "H3K4me3_regions.bed"
  # available at OSF files: 02_cCREs_OCR_gene_links/public_epigenetic_data_for_validation/histone_ChiPseq/h3k4me3/H3K4me3_regions.bed
# H3K27ac regions: "H3K27ac_regions.bed"
  # available at OSF files: 02_cCREs_OCR_gene_links/public_epigenetic_data_for_validation/histone_ChiPseq/h3k27ac/H3K27ac_regions.bed

# Auxiliar functions
source("auxiliar_functions.R")


load("bcell_bwplots.RData", verbose = TRUE)
consensus_OCRs_granges <- import("consensus_OCRs_all_cell_subpopulations.bed", format = "BED")
h3k4me1 <- import("H3K4me1_regions.bed", format = "BED")
h3k4me3 <- import("H3K4me3_regions.bed", format = "BED")
h3k27ac <- import("H3K27ac_regions.bed", format = "BED")

# Bigwig plot
#rs35837782 (10:124604740): chr10:124,602,438-124,608,467
snp <- makeGRangesFromDataFrame(df=data.frame(chr="chr10",start=124604740,end=124604740))

gene <- "rs35837782" 
bw_plot3(accDT_by_group, peaks=consensus_OCRs_granges, gene=gene, ylims=c(0,400), start=124602438, end=124608467, chr="chr10",
         h3k4me1, h3k4me3, h3k27ac, snp)

###-----------------------------------------------------------------------------------###




#***************************************************
## Figure 4 K
#***************************************************
## Density plot of gene expression by cell subpopulations for LHPP, FAM53B-AS1 and CHST15 gene.

library(ggplot2)
library(ggridges)

# Input data:
# Normalized RNA: "rna_norm_data.txt" 
  # available at OSF files: 01_rna_seq_data/rna_norm_data.txt
# Metadata RNA: "rna_metadata.txt" 
  # available at OSF files: 01_rna_seq_data/rna_metadata.txt
# Gene annotation: "rna_gene_annotation.txt"
  # available at OSF files: 01_rna_seq_data/rna_gene_annotation.txt
# RNA expression colors: "rna_expression_colors.rds"


# RNA - Ridgeplot 

rna_norm <- read.table("rna_norm_data.txt",header=TRUE, dec=".", sep="\t")
rna_metadata <- read.table("rna_metadata.txt",header=TRUE, dec=".", sep="\t")
rna_annotation <- read.table("rna_gene_annotation.txt",header=TRUE, dec=".", sep="\t")
rna_color <- readRDS("rna_expression_colors.rds")

group <- factor(as.character(rna_metadata$CellType),
                levels=c("Naive_CD5pos","Naive_CD5neg","Transitional_B","ImmatureB","preB","proB","CLP","HSC"))

#LHPP
data_to_plot <- data.frame(cell_type=group, score=as.numeric(rna_norm[rna_annotation$ensembl_gene_id[rna_annotation$hgnc_symbol=="LHPP"],]))
ggplot(data_to_plot, aes(x = score, y = cell_type, fill= ..x..)) +
  geom_density_ridges_gradient(scale=2) +
  scale_fill_gradientn(colors = rna_color) +
  labs(title = "LHPP - RNA expression (Z-score)") +
  theme_ridges()

#FAM53B-AS1
data_to_plot <- data.frame(cell_type=group, score=as.numeric(rna_norm[rna_annotation$ensembl_gene_id=="ENSG00000278831",]))
ggplot(data_to_plot, aes(x = score, y = cell_type, fill= ..x..)) +
  geom_density_ridges_gradient(scale=2) +
  scale_fill_gradientn(colors = rna_color) +
  labs(title = "FAM53B-AS1 - RNA expression (Z-score)") +
  theme_ridges()

#CHST15
data_to_plot <- data.frame(cell_type=group, score=as.numeric(rna_norm[rna_annotation$ensembl_gene_id[rna_annotation$hgnc_symbol=="CHST15"],]))
ggplot(data_to_plot, aes(x = score, y = cell_type, fill= ..x..)) +
  geom_density_ridges_gradient(scale=2) +
  scale_fill_gradientn(colors = rna_color) +
  labs(title = "CHST15 - RNA expression (Z-score)") +
  theme_ridges()

###-----------------------------------------------------------------------------------###




#***************************************************
## Figure 4 L
#***************************************************
## Violin plot of LHPP and CHST15 gene expression by BCP ALL subtypes

library(ggplot2)
library(BuenColors)

# Input data:
# B-ALL normalized rna data: "BALL_rna_norm_data.txt"
  # available at OSF files: 05_B-ALL/BALL_rna_norm_data.txt
# B-ALL metadata: "BALL_metadata.txt"
  # available at OSF files: 05_B-ALL/BALL_metadata.txt
# Gene annotation: "rna_gene_annotation.txt"
  # available at OSF files: 01_rna_seq_data/rna_gene_annotation.txt

# Auxiliar functions
source("auxiliar_functions.R")

rna_annotation <- read.table("rna_gene_annotation.txt",header=TRUE, dec=".", sep="\t")
BALL_metadata <- read.table("BALL_metadata.txt", sep="\t", dec=".", header=TRUE)
BALL_rna_norm_data <- read.table("BALL_rna_norm_data.txt", sep="\t", dec=".", header=TRUE)

#LHPP
data_to_plot <- data.frame(gene=as.numeric(BALL_rna_norm_data[rna_annotation$ensembl_gene_id[rna_annotation$hgnc_symbol=="LHPP"],]),group=BALL_metadata$subgroup)

print(ggplot(data_to_plot, aes(x=group, y=gene, fill=group)) + 
        geom_violin(trim=FALSE) +
        scale_fill_manual(values=jdb_palette("corona")) +
        stat_summary(fun.data=data_summary) +
        theme_classic())

#CHST15
data_to_plot <- data.frame(gene=as.numeric(BALL_rna_norm_data[rna_annotation$ensembl_gene_id[rna_annotation$hgnc_symbol=="CHST15"],]),group=BALL_metadata$subgroup)

print(ggplot(data_to_plot, aes(x=group, y=gene, fill=group)) + 
        geom_violin(trim=FALSE) +
        scale_fill_manual(values=jdb_palette("corona")) +
        stat_summary(fun.data=data_summary) +
        theme_classic())

###-----------------------------------------------------------------------------------###
