#***************************************************
## Paper: "Uncovering the cis-regulatory program of early human B-cell commitment and its implications
## in the pathogenesis of acute lymphoblastic leukemia"
## Authors: Planell N et al.
## Date: 2023
## Code: Figure 1
## Input data available at: 
#***************************************************




#***************************************************
## Figure 1 A
#***************************************************
## Schematic picture
#***************************************************




#***************************************************
## Figure 1 B
#***************************************************
## Principal Component Analysis (PCA) for RNA-Seq and ATAC-Seq
library(ggplot2)

# PCA of RNA-Seq
# Input data:
  # Normalized RNA: "rna_norm_data.txt" 
      # available at OSF files: 01_rna_seq_data/rna_norm_data.txt
  # Metadata RNA: "rna_metadata.txt" 
      # available at OSF files: 01_rna_seq_data/rna_metadata.txt
  # B cell subpopulation colors: "bcell_colors.rds"

rna_norm <- read.table("rna_norm_data.txt",header=TRUE, dec=".", sep="\t")
rna_metadata <- read.table("rna_metadata.txt",header=TRUE, dec=".", sep="\t")
bcell_colors <- readRDS("bcell_colors.rds")

# PCA computation
pca_rna <- prcomp(t(rna_norm))
var_prcomp <- pca_rna$sdev^2

# PCs variance
pcvar <- data.frame(var=var_prcomp/sum(var_prcomp), 
                    pc=c(1:length(var_prcomp)))

# Plot
data_to_plot <- data.frame(PC1=pca_rna$x[,1],PC2=pca_rna$x[,2], CellType=rna_metadata$CellType)

p <- ggplot(data_to_plot, aes(x = PC1, y = PC2, color = CellType)) +
  geom_point(size = 3, show.legend=FALSE) + coord_fixed(ratio=1.35) + scale_color_manual(values = bcell_colors) + directlabels::geom_dl(aes(label = CellType), method = list("smart.grid", cex=1)) + 
  labs(title="PCA: RNA - Log2(Normalized and Filtered Count table (VOOM))",x = paste("PC1 (",round(pcvar[1,1]*100,0),"%)", sep=""), y = paste("PC2 (",round(pcvar[2,1]*100,0),"%)", sep=""))+
  theme_classic()


# PCA of ATAC-Seq
# Input data:
# Normalized ATAC: "atac_norm_data.txt" 
  # available at OSF files: 01_atac_seq_data/atac_norm_data.txt
# Metadata ATAC: "atac_metadata.txt"
  # available at OSF files: 01_atac_seq_data/atac_metadata.txt
# B cell subpopulation colors: "bcell_colors.rds"

atac_norm <- read.table("atac_norm_data.txt",header=TRUE, dec=".", sep="\t") 
atac_metadata <- read.table("atac_metadata.txt",header=TRUE, dec=".", sep="\t") 
bcell_colors <- readRDS("bcell_colors.rds")

# PCA computation
pca_atac <- prcomp(t(atac_norm))
var_prcomp <- pca_atac$sdev^2

# PCs variance
pcvar <- data.frame(var=var_prcomp/sum(var_prcomp), 
                    pc=c(1:length(var_prcomp)))

# Plot
data_to_plot <- data.frame(PC1=pca_atac$x[,1],PC2=pca_atac$x[,2], CellType=atac_metadata$CellType)

p <- ggplot(data_to_plot, aes(x = PC1, y = PC2, color = CellType)) +
  geom_point(size = 3, show.legend=FALSE) + coord_fixed(ratio=1.35) + scale_color_manual(values = bcell_colors) + directlabels::geom_dl(aes(label = CellType), method = list("smart.grid", cex=1)) + 
  labs(title="PCA: ATAC - Normalized and Filtered Count table (TMM)",x = paste("PC1 (",round(pcvar[1,1]*100,0),"%)", sep=""), y = paste("PC2 (",round(pcvar[2,1]*100,0),"%)", sep=""))+
  theme_classic()

#***************************************************




#***************************************************
## Figure 1 C
#***************************************************
## Bigwig representation for PAX5 region: 36,980,000-37,042,948
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
gene <- "PAX5"
bw_plot2(accDT_by_group, peaks=consensus_OCRs_granges, gene=gene, ylims=c(0,800), start=36980000, end=37042948, chr="chr9",
         h3k4me1, h3k4me3, h3k27ac)


# Gene expression barplot
# Input data:
# Normalized RNA: "rna_norm_data.txt" 
  # available at OSF files: 01_rna_seq_data/rna_norm_data.txt
# Metadata RNA: "rna_metadata.txt" 
  # available at OSF files: 01_rna_seq_data/rna_metadata.txt
# Gene annotation: "rna_gene_annotation.txt"
  # available at OSF files: 01_rna_seq_data/rna_gene_annotation.txt
# B cell subpopulation colors: "bcell_colors.rds"

rna_norm <- read.table("rna_norm_data.txt",header=TRUE, dec=".", sep="\t")
rna_metadata <- read.table("rna_metadata.txt",header=TRUE, dec=".", sep="\t")
rna_annotation <- read.table("rna_gene_annotation.txt",header=TRUE, dec=".", sep="\t")
bcell_colors <- readRDS("bcell_colors.rds")

rna_metadata$CellType <- factor(rna_metadata$CellType, levels=unique(rna_metadata$CellType))

target_gene <- "PAX5"
data_to_plot <- data.frame(value=as.numeric(by(as.numeric(rna_norm[rna_annotation$ensembl_gene_id[rna_annotation$hgnc_symbol==target_gene],]), rna_metadata$CellType, median)),
                           class=factor(levels(rna_metadata$CellType),levels = levels(rna_metadata$CellType)))
p <- ggplot(data_to_plot, aes(x=class, y=value, fill=class)) + 
  geom_bar(stat="identity")+
  scale_fill_manual(values=bcell_colors) +
  theme_bw()

#***************************************************




#***************************************************
## Figure 1 D, F, G
#***************************************************
## Relation between gene expression and chromatin accessibility

library(ggplot2)

# Bigwig image
# Input data:
# Normalized RNA: "rna_norm_data.txt" 
  # available at OSF files: 01_rna_seq_data/rna_norm_data.txt
# Metadata RNA: "rna_metadata.txt" 
  # available at OSF files: 01_rna_seq_data/rna_metadata.txt
# Gene annotation: "rna_gene_annotation.txt"
  # available at OSF files: 01_rna_seq_data/rna_gene_annotation.txt
# Normalized ATAC: "atac_norm_data.txt" 
  # available at OSF files: 01_atac_seq_data/atac_norm_data.txt
# Metadata ATAC: "atac_metadata.txt"
  # available at OSF files: 01_atac_seq_data/atac_metadata.txt
# OCR annotation: "atac_OCR_annotation.txt"
  # available at OSF files: 01_atac_seq_data/atac_OCR_annotation.txt
# B cell subpopulation colors: "bcell_colors.rds"
# Common individuals: "paired_samples.rds"


rna_norm <- read.table("rna_norm_data.txt",header=TRUE, dec=".", sep="\t")
rna_metadata <- read.table("rna_metadata.txt",header=TRUE, dec=".", sep="\t")
rna_annotation <- read.table("rna_gene_annotation.txt",header=TRUE, dec=".", sep="\t")

atac_norm <- read.table("atac_norm_data.txt",header=TRUE, dec=".", sep="\t") 
atac_metadata <- read.table("atac_metadata.txt",header=TRUE, dec=".", sep="\t") 
atac_annotation <- read.table("atac_OCR_annotation.txt",header=TRUE, dec=".", sep="\t", quote="")


bcell_colors <- readRDS("bcell_colors.rds")
paired_samples <- readRDS("paired_samples.rds")


#Get paired data
rna2 <- rna_norm[,paired_samples]
rna_metadata2 <- rna_metadata[paired_samples,]
atac2 <- atac_norm[,paired_samples]


#PAX5 Promoter peak: chr9:37,034,185-37,035,344 (panel D)

data_to_plot <- data.frame(rna=as.numeric(rna2[rna_annotation$ensembl_gene_id[rna_annotation$hgnc_symbol=="PAX5"],]),
                           atac=as.numeric(atac2["chr9_37034185_37035344",]),
                           cell_type=rna_metadata2$CellType)


ggplot(data_to_plot, aes(x=rna,y=atac,colour=cell_type)) + 
  geom_point(size=3) + 
  scale_color_manual(values=bcell_colors) +
  theme_bw()

#PAX5 Enhancer peak: chr9:36,993,757-36,994,690 (panel F)

data_to_plot <- data.frame(rna=as.numeric(rna2[rna_annotation$ensembl_gene_id[rna_annotation$hgnc_symbol=="PAX5"],]),
                           atac=as.numeric(atac2["chr9_36993757_36994690",]),
                           cell_type=rna_metadata2$CellType)


ggplot(data_to_plot, aes(x=rna,y=atac,colour=cell_type)) + 
  geom_point(size=3) + 
  scale_color_manual(values=bcell_colors) +
  theme_bw()


#PAX5 Enhancer peak: chr9:chr9:37,566,825-37,567,195 (panel G)

data_to_plot <- data.frame(rna=as.numeric(rna2[rna_annotation$ensembl_gene_id[rna_annotation$hgnc_symbol=="PAX5"],]),
                           atac=as.numeric(atac2["chr9_37566825_37567195",]),
                           cell_type=rna_metadata2$CellType)


ggplot(data_to_plot, aes(x=rna,y=atac,colour=cell_type)) + 
  geom_point(size=3) + 
  scale_color_manual(values=bcell_colors) +
  theme_bw()


#PAX5 Enhancer peak: chr9:37,029,060-37,029,688 (panel G)

data_to_plot <- data.frame(rna=as.numeric(rna2[rna_annotation$ensembl_gene_id[rna_annotation$hgnc_symbol=="PAX5"],]),
                           atac=as.numeric(atac2["chr9_37029060_37029688",]),
                           cell_type=rna_metadata2$CellType)


ggplot(data_to_plot, aes(x=rna,y=atac,colour=cell_type)) + 
  geom_point(size=3) + 
  scale_color_manual(values=bcell_colors) +
  theme_bw()


#PAX5 Enhancer peak: chr9:37,065,987-37,066,673 (panel G)

data_to_plot <- data.frame(rna=as.numeric(rna2[rna_annotation$ensembl_gene_id[rna_annotation$hgnc_symbol=="PAX5"],]),
                           atac=as.numeric(atac2["chr9_37065987_37066673",]),
                           cell_type=rna_metadata2$CellType)


ggplot(data_to_plot, aes(x=rna,y=atac,colour=cell_type)) + 
  geom_point(size=3) + 
  scale_color_manual(values=bcell_colors) +
  theme_bw()

#***************************************************




#***************************************************
## Figure 1 E
#***************************************************
## Heatmap representation of chromatin accessibility for CREs identified for PAX5

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
# Normalized ATAC: "atac_norm_data.txt" 
  # available at OSF files: 01_atac_seq_data/atac_norm_data.txt
# Metadata ATAC: "atac_metadata.txt"
  # available at OSF files: 01_atac_seq_data/atac_metadata.txt
# OCR annotation: "atac_OCR_annotation.txt"
  # available at OSF files: 01_atac_seq_data/atac_OCR_annotation.txt
# List of cCREs: "candidate_CREs_PCHiC_links_validation.txt"
  # available at OSF files: 02_cCREs_OCR_gene_links/public_epigenetic_data_for_validation/candidate_CREs_PCHiC_links_validation.txt
# B cell subpopulation colors: "bcell_colors.rds"
# Common individuals: "paired_samples.rds"


rna_norm <- read.table("rna_norm_data.txt",header=TRUE, dec=".", sep="\t")
rna_metadata <- read.table("rna_metadata.txt",header=TRUE, dec=".", sep="\t")
rna_annotation <- read.table("rna_gene_annotation.txt",header=TRUE, dec=".", sep="\t")

atac_norm <- read.table("atac_norm_data.txt",header=TRUE, dec=".", sep="\t") 
atac_metadata <- read.table("atac_metadata.txt",header=TRUE, dec=".", sep="\t") 
atac_annotation <- read.table("atac_OCR_annotation.txt",header=TRUE, dec=".", sep="\t", quote="")

paired_samples <- readRDS("paired_samples.rds")

ocr_gene <- read.table("candidate_CREs_PCHiC_links_validation.txt",header=TRUE, dec=".", sep="\t")
ocr_gene_pax5 <- ocr_gene[ocr_gene$hgnc_symbol=="PAX5",]


rownames(ocr_gene_pax5) <- ocr_gene_pax5$peak
ocr_gene_pax5$pchic <- rep("no", nrow(ocr_gene_pax5))
ocr_gene_pax5$pchic[which(ocr_gene_pax5$overlap_PR_Biola_ourgene_in_baitName=="TRUE" | ocr_gene_pax5$overlap_ER_Biola_ourgene_in_baitName=="TRUE")] <- "yes"

ocr_gene_pax5$h3k4me3 <- atac_annotation[ocr_gene_pax5$peak,"h3k4me3"]
ocr_gene_pax5$h3k4me1 <- atac_annotation[ocr_gene_pax5$peak,"h3k4me1"]
ocr_gene_pax5$h3k27ac <- atac_annotation[ocr_gene_pax5$peak,"h3k27ac"]

#Graphical parameters
#Peak type color for annotation
peak_colors <- as.character(jdb_palette("Darjeeling2",2))
names(peak_colors) <- unique(ocr_gene_pax5$peak_type)
#B-cell color for annotation
bcell_colors <- readRDS("bcell_colors.rds")
#TSS dist color for annotation
col_fun = colorRamp2(c(-1000000,-500000,-334000,-167000, 0, 167000, 334000, 500000, 1000000), jdb_palette("ocean_green",9,"discrete"))
#rho color for annotation
col_fun1 = colorRamp2(c(-0.9,-0.8,-0.7,-0.6, 0, 0.6, 0.7, 0.8, 0.9), jdb_palette("solar_flare",9,"discrete"))
#log10pval color for annotation
col_fun2 = colorRamp2(c(0, 2.5, 5), c("white", "gray", "black"))
#Histone marks color for annotation
hist_colors = c("aliceblue","#B980A7")
names(hist_colors) <- c("no","yes")
#PCHiC marks color for annotation
pchic_colors = c("aliceblue","#F4711F")
names(pchic_colors) <- c("no","yes")

#Plot of ATAC + RNA of CREs from a target gene
## Column annotation
ha = HeatmapAnnotation(stage = rna_metadata[paired_samples,"CellType"], 
                       col = list(stage = bcell_colors),
                       rna = anno_barplot(as.numeric(rna_norm[rna_annotation$ensembl_gene_id[rna_annotation$hgnc_symbol=="PAX5"],paired_samples])))

## Row annotation
rha = rowAnnotation(peak = ocr_gene_pax5$peak_type, 
                    TSS_dist = ocr_gene_pax5$dist,
                    rho = ocr_gene_pax5$spearman_rho,
                    log10_pval = (log10(ocr_gene_pax5$spearman_adj.pval+0.00001)*-1),
                    h3k4me3 = ocr_gene_pax5$h3k4me3,
                    h3k4me1 = ocr_gene_pax5$h3k4me1,
                    h3k27ac = ocr_gene_pax5$h3k27ac,
                    pchic = ocr_gene_pax5$pchic,
                    col = list(peak = peak_colors,
                               TSS_dist = col_fun,
                               rho = col_fun1,
                               log10_pval = col_fun2,
                               h3k4me3 = hist_colors,
                               h3k4me1 = hist_colors,
                               h3k27ac = hist_colors,
                               pchic = pchic_colors))

## OCR profile
f1 = jdb_palette("samba_color")
M <- log2(atac_norm[ocr_gene_pax5$peak, paired_samples]+1)

p <- Heatmap(M, top_annotation = ha, right_annotation=rha, 
                  show_row_names = FALSE, show_column_names =FALSE, 
                  column_names_gp = gpar(fontsize = 6), col=f1, 
                  heatmap_legend_param = list(title = "OCRs measure"),
                  column_order = order(rna_metadata[paired_samples,"CellType"]))

#***************************************************




#***************************************************
## Figure 1 H
#***************************************************
## Distribution of chromatin accessibility for PAX5 cCREs

library(ggplot2)

# Bigwig image
# Input data:
# Normalized ATAC: "atac_norm_data.txt" 
  # available at OSF files: 01_atac_seq_data/atac_norm_data.txt
# Metadata ATAC: "atac_metadata.txt"
  # available at OSF files: 01_atac_seq_data/atac_metadata.txt
# List of cCREs: "candidate_CREs_PCHiC_links_validation.txt"
  # available at OSF files: 02_cCREs_OCR_gene_links/public_epigenetic_data_for_validation/candidate_CREs_PCHiC_links_validation.txt


atac_norm <- read.table("atac_norm_data.txt",header=TRUE, dec=".", sep="\t") 
atac_metadata <- read.table("atac_metadata.txt",header=TRUE, dec=".", sep="\t") 

ocr_gene <- read.table("candidate_CREs_PCHiC_links_validation.txt",header=TRUE, dec=".", sep="\t")
ocr_gene_pax5 <- ocr_gene[ocr_gene$hgnc_symbol=="PAX5",]
  
list_OCRs <- unique(ocr_gene_pax5$peak[which(ocr_gene_pax5$overlap_PR_Biola_ourgene_in_baitName==TRUE | ocr_gene_pax5$overlap_ER_Biola_ourgene_in_baitName==TRUE)])
naive_accessibility_PCHiC <- unlist(atac_norm[list_OCRs,grep("Naive",atac_metadata$CellType)])


list_OCRs <- unique(ocr_gene_pax5$peak[which(ocr_gene_pax5$overlap_PR_Biola_ourgene_in_baitName==FALSE & ocr_gene_pax5$overlap_ER_Biola_ourgene_in_baitName==FALSE)])
naive_accessibility_no_PCHiC <- unlist(atac_norm[list_OCRs,grep("Naive",atac_metadata$CellType)])



data_to_plot <- data.frame(value=c(naive_accessibility_PCHiC,naive_accessibility_no_PCHiC), 
                           type=c(rep("validated",length(naive_accessibility_PCHiC)),rep("unvalidated",length(naive_accessibility_no_PCHiC))))

p <- ggplot(data_to_plot, aes(x=value, fill=type)) +
  geom_density(alpha=0.4) +
  scale_fill_grey() +
  theme_bw()

#***************************************************
