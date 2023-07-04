#***************************************************
## Paper: "Uncovering the cis-regulatory program of early human B-cell commitment and its implications
## in the pathogenesis of acute lymphoblastic leukemia"
## Authors: Planell N et al.
## Date: 2023
## Code: Figure 5
## Input data available at: 
#***************************************************




#***************************************************
## Figure 5 A
#***************************************************
##	Dot plot showing the promoter activation. 

require(ComplexHeatmap)
require(ggplot2)
require(RColorBrewer)
require(circlize)

# Input data:
# Normalized RNA: "rna_norm_data.txt" 
  # available at OSF files: 01_rna_seq_data/rna_norm_data.txt
# Metadata RNA: "rna_metadata.txt" 
  # available at OSF files: 01_rna_seq_data/rna_metadata.txt
# Gene annotation: "rna_gene_annotation.txt"
  # available at OSF files: 01_rna_seq_data/rna_gene_annotation.txt
# Differential expression analysis: "rna_dea_transitions.txt"
  # available at OSF files: 01_rna_seq_data/na_dea_transitions.txt
# Normalized ATAC: "atac_norm_data.txt" 
  # available at OSF files: 01_atac_seq_data/atac_norm_data.txt
# Metadata ATAC: "atac_metadata.txt"
  # available at OSF files: 01_atac_seq_data/atac_metadata.txt
# OCR annotation: "atac_OCR_annotation.txt"
  # available at OSF files: 01_atac_seq_data/atac_OCR_annotation.txt
# Differential binding analysis: "atac_dba_transitions.txt"
  # available at OSF files: 01_atac_seq_data/atac_dba_transitions.txt
# B cell subpopulation colors: "bcell_colors.rds"


rna_norm <- read.table("rna_norm_data.txt",header=TRUE, dec=".", sep="\t")
rna_metadata <- read.table("rna_metadata.txt",header=TRUE, dec=".", sep="\t")
rna_annotation <- read.table("rna_gene_annotation.txt",header=TRUE, dec=".", sep="\t")
rna_dea <- read.table("rna_dea_transitions.txt",header=TRUE, dec=".", sep="\t")
rna_clust <- read.table("rna_deg_clustering.txt",header=TRUE, dec=".", sep="\t")
  
transitions <- grep("logFC",colnames(rna_dea),value=TRUE)
transitions_sig <- grep("sig_",colnames(rna_dea),value=TRUE)
rna_dea[,c(transitions,transitions_sig)] <- rna_dea[,c(transitions,transitions_sig)]*-1


atac_norm <- read.table("atac_norm_data.txt",header=TRUE, dec=".", sep="\t") 
atac_metadata <- read.table("atac_metadata.txt",header=TRUE, dec=".", sep="\t") 
atac_annotation <- read.table("atac_OCR_annotation.txt",header=TRUE, dec=".", sep="\t", quote="")
atac_dba <- read.table("atac_dba_transitions.txt",header=TRUE, dec=".", sep="\t")


bcell_colors <- readRDS("bcell_colors.rds")


bcell_tf_ocr_gene <- read.table("tf_ocr_gene_interactions_curated.txt",header=TRUE, dec=".", sep="\t")

bcell_tf_ocr_gene2 <- unique(bcell_tf_ocr_gene[,c("ocr","gene_ensembl")]) #28,091 relations
CREs_promoter <- bcell_tf_ocr_gene2[bcell_tf_ocr_gene2$ocr%in%rownames(atac_annotation)[atac_annotation$annotation_simplified=="Promoter"],] #6,310

# Identify increasing and decreasing genes over differentiation
summary(rna_clust$pattern_id)

clust_up <- vector()
clust_down <- vector()
for(i in 1:max(rna_clust$pattern_id)){
  sel_genes <- rownames(rna_clust)[rna_clust$pattern_id==i]
  aa <- rna_dea[sel_genes,grep("sig",colnames(rna_dea))]
  direction <- median(rowSums(aa))
  if(direction>0){clust_up <- c(clust_up,i)}
  if(direction<0){clust_down <- c(clust_down,i)}
}

rna_up <- rna_clust[rna_clust$pattern_id%in%clust_up,]
rna_down <- rna_clust[rna_clust$pattern_id%in%clust_down,]


# TRANSITION WITH TOP SIGNIFICANT CHANGE

# Promoter Increasing
sel_rna_transition <- vector()
sel_atac_transition <- vector()
transition_gene <- vector()
transition_OCR <- vector()

for(i in 1:nrow(rna_up)){
  if(rownames(rna_up)[i]%in%CREs_promoter$gene_ensembl){
    sel_peaks <- CREs_promoter$ocr[which(CREs_promoter$gene_ensembl==rownames(rna_up)[i])]
    for(j in 1:length(sel_peaks)){
      ss <- which(atac_dba[sel_peaks[j],transitions_sig]==1)
      ss_rna <- which(rna_dea[rownames(rna_up)[i],transitions_sig]==1)
      if(length(ss)==0){next}
      if(length(ss_rna)==0){next}
      if(length(ss)>1){ ss <- ss[which(atac_dba[sel_peaks[j],transitions[ss]]==max(atac_dba[sel_peaks[j],transitions[ss]]))]}
      if(length(ss_rna)>1){ss_rna <- ss_rna[which(rna_dea[rownames(rna_up)[i],transitions[ss_rna]]==max(rna_dea[rownames(rna_up)[i],transitions[ss_rna]]))]}
      sel_rna_transition <- c(sel_rna_transition, ss_rna) #the transition where a the maximum sig increase of RNA occurs
      sel_atac_transition <- c(sel_atac_transition, ss) #the transition where a the maximum sig increaes of promoter accessibility occur
      transition_gene <- c(transition_gene,rownames(rna_up)[i])
      transition_OCR <- c(transition_OCR,sel_peaks[j])
    }
  }
}

regulation_up <- data.frame(gene=transition_gene,OCR=transition_OCR,rna=sel_rna_transition,atac=sel_atac_transition)

#Plot
regulation_up_plot <- data.frame(rna=factor(rep(gsub("logFC_","",transitions)[1:max(sel_rna_transition)],max(sel_atac_transition)), levels=gsub("logFC_","",transitions)[1:max(sel_rna_transition)]),
                                 atac=factor(rep(gsub("logFC_","",transitions)[1:max(sel_atac_transition)], each=max(sel_rna_transition)), levels=gsub("logFC_","",transitions)[1:max(sel_atac_transition)]),
                                 values=as.vector(table(regulation_up$rna, regulation_up$atac)))


p_up <- ggplot(regulation_up_plot, aes(x=atac, y=rna, size = values, color=values)) +
  geom_point(alpha=0.7) +
  scale_size(range = c(.1, 18), name="Genes") +
  scale_y_discrete(limits = rev) +
  scale_color_gradient(low="blue", high="red") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


# Promoter Decreasing
sel_rna_transition <- vector()
sel_atac_transition <- vector()
transition_gene <- vector()
transition_OCR <- vector()


for(i in 1:nrow(rna_down)){
  if(rownames(rna_down)[i]%in%CREs_promoter$gene_ensembl){
    sel_peaks <- CREs_promoter$ocr[which(CREs_promoter$gene_ensembl==rownames(rna_down)[i])]
    for(j in 1:length(sel_peaks)){
      ss <- which(atac_dba[sel_peaks[j],transitions_sig]==-1)
      ss_rna <- which(rna_dea[rownames(rna_down)[i],transitions_sig]==-1)
      if(length(ss)==0){next}
      if(length(ss_rna)==0){next}
      if(length(ss)>1){ss <- ss[which(atac_dba[sel_peaks[j],transitions[ss]]==min(atac_dba[sel_peaks[j],transitions[ss]]))]}
      if(length(ss_rna)>1){ss_rna <- ss_rna[which(rna_dea[rownames(rna_down)[i],transitions[ss_rna]]==min(rna_dea[rownames(rna_down)[i],transitions[ss_rna]]))]}
      sel_rna_transition <- c(sel_rna_transition, ss_rna) #the first transition where a sig idecreases of RNA occurs
      sel_atac_transition <- c(sel_atac_transition, ss) #the first transition where a sig decreaes of promoter accessibility occur
      transition_gene <- c(transition_gene,rownames(rna_down)[i])
      transition_OCR <- c(transition_OCR,sel_peaks[j])
    }
  }
}


regulation_down <- data.frame(gene=transition_gene,OCR=transition_OCR,rna=sel_rna_transition,atac=sel_atac_transition)

#Plot
regulation_down_plot <- data.frame(rna=factor(rep(gsub("logFC_","",transitions)[1:max(sel_rna_transition)],max(sel_atac_transition)), levels=gsub("logFC_","",transitions)[1:max(sel_rna_transition)]),
                                   atac=factor(rep(gsub("logFC_","",transitions)[1:max(sel_atac_transition)], each=max(sel_rna_transition)), levels=gsub("logFC_","",transitions)[1:max(sel_atac_transition)]),
                                   values=as.vector(table(regulation_down$rna, regulation_down$atac)))


p_down <- ggplot(regulation_down_plot, aes(x=atac, y=rna, size = values, color=values)) +
  geom_point(alpha=0.7) +
  scale_size(range = c(.1, 18), name="Genes") +
  scale_y_discrete(limits = rev) +
  scale_color_gradient(low="blue", high="red") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#***************************************************




#***************************************************
## Figure 5 B
#***************************************************
##	Dot plot showing the enhancer activation. 

require(ComplexHeatmap)
require(ggplot2)
require(RColorBrewer)
require(circlize)

# Input data:
# Normalized RNA: "rna_norm_data.txt" 
  # available at OSF files: 01_rna_seq_data/rna_norm_data.txt
# Metadata RNA: "rna_metadata.txt" 
  # available at OSF files: 01_rna_seq_data/rna_metadata.txt
# Gene annotation: "rna_gene_annotation.txt"
  # available at OSF files: 01_rna_seq_data/rna_gene_annotation.txt
# Differential expression analysis: "rna_dea_transitions.txt"
  # available at OSF files: 01_rna_seq_data/na_dea_transitions.txt
# Normalized ATAC: "atac_norm_data.txt" 
  # available at OSF files: 01_atac_seq_data/atac_norm_data.txt
# Metadata ATAC: "atac_metadata.txt"
  # available at OSF files: 01_atac_seq_data/atac_metadata.txt
# OCR annotation: "atac_OCR_annotation.txt"
  # available at OSF files: 01_atac_seq_data/atac_OCR_annotation.txt
# Differential binding analysis: "atac_dba_transitions.txt"
  # available at OSF files: 01_atac_seq_data/atac_dba_transitions.txt
# B cell subpopulation colors: "bcell_colors.rds"


rna_norm <- read.table("rna_norm_data.txt",header=TRUE, dec=".", sep="\t")
rna_metadata <- read.table("rna_metadata.txt",header=TRUE, dec=".", sep="\t")
rna_annotation <- read.table("rna_gene_annotation.txt",header=TRUE, dec=".", sep="\t")
rna_dea <- read.table("rna_dea_transitions.txt",header=TRUE, dec=".", sep="\t")
rna_clust <- read.table("rna_deg_clustering.txt",header=TRUE, dec=".", sep="\t")

transitions <- grep("logFC",colnames(rna_dea),value=TRUE)
transitions_sig <- grep("sig_",colnames(rna_dea),value=TRUE)
rna_dea[,c(transitions,transitions_sig)] <- rna_dea[,c(transitions,transitions_sig)]*-1


atac_norm <- read.table("atac_norm_data.txt",header=TRUE, dec=".", sep="\t") 
atac_metadata <- read.table("atac_metadata.txt",header=TRUE, dec=".", sep="\t") 
atac_annotation <- read.table("atac_OCR_annotation.txt",header=TRUE, dec=".", sep="\t", quote="")
atac_dba <- read.table("atac_dba_transitions.txt",header=TRUE, dec=".", sep="\t")


bcell_colors <- readRDS("bcell_colors.rds")


bcell_tf_ocr_gene <- read.table("tf_ocr_gene_interactions_curated.txt",header=TRUE, dec=".", sep="\t")

bcell_tf_ocr_gene2 <- unique(bcell_tf_ocr_gene[,c("ocr","gene_ensembl")]) #28,091 relations
CREs_enhancer <- bcell_tf_ocr_gene2[!bcell_tf_ocr_gene2$ocr%in%rownames(atac_annotation)[atac_annotation$annotation_simplified=="Promoter"],] #21,781

# Identify increasing and decreasing genes over differentiation
summary(rna_clust$pattern_id)

clust_up <- vector()
clust_down <- vector()
for(i in 1:max(rna_clust$pattern_id)){
  sel_genes <- rownames(rna_clust)[rna_clust$pattern_id==i]
  aa <- rna_dea[sel_genes,grep("sig",colnames(rna_dea))]
  direction <- median(rowSums(aa))
  if(direction>0){clust_up <- c(clust_up,i)}
  if(direction<0){clust_down <- c(clust_down,i)}
}

rna_up <- rna_clust[rna_clust$pattern_id%in%clust_up,]
rna_down <- rna_clust[rna_clust$pattern_id%in%clust_down,]


# TRANSITION WITH TOP SIGNIFICANT CHANGE

# Enhancer Increasing
sel_rna_transition <- vector()
sel_atac_transition <- vector()
transition_gene <- vector()
transition_OCR <- vector()

for(i in 1:nrow(rna_up)){
  if(rownames(rna_up)[i]%in%CREs_enhancer$gene_ensembl){
    sel_peaks <- CREs_enhancer$ocr[which(CREs_enhancer$gene_ensembl==rownames(rna_up)[i])]
    for(j in 1:length(sel_peaks)){
      ss <- which(atac_dba[sel_peaks[j],transitions_sig]==1)
      ss_rna <- which(rna_dea[rownames(rna_up)[i],transitions_sig]==1)
      if(length(ss)==0){next}
      if(length(ss_rna)==0){next}
      if(length(ss)>1){ ss <- ss[which(atac_dba[sel_peaks[j],transitions[ss]]==max(atac_dba[sel_peaks[j],transitions[ss]]))]}
      if(length(ss_rna)>1){ss_rna <- ss_rna[which(rna_dea[rownames(rna_up)[i],transitions[ss_rna]]==max(rna_dea[rownames(rna_up)[i],transitions[ss_rna]]))]}
      sel_rna_transition <- c(sel_rna_transition, ss_rna) #the transition where a the maximum sig increase of RNA occurs
      sel_atac_transition <- c(sel_atac_transition, ss) #the transition where a the maximum sig increaes of promoter accessibility occur
      transition_gene <- c(transition_gene,rownames(rna_up)[i])
      transition_OCR <- c(transition_OCR,sel_peaks[j])
    }
  }
}

regulation_up <- data.frame(gene=transition_gene,OCR=transition_OCR,rna=sel_rna_transition,atac=sel_atac_transition)

#Plot
regulation_up_plot <- data.frame(rna=factor(rep(gsub("logFC_","",transitions)[1:max(sel_rna_transition)],max(sel_atac_transition)), levels=gsub("logFC_","",transitions)[1:max(sel_rna_transition)]),
                                 atac=factor(rep(gsub("logFC_","",transitions)[1:max(sel_atac_transition)], each=max(sel_rna_transition)), levels=gsub("logFC_","",transitions)[1:max(sel_atac_transition)]),
                                 values=as.vector(table(regulation_up$rna, regulation_up$atac)))


p_up <- ggplot(regulation_up_plot, aes(x=atac, y=rna, size = values, color=values)) +
  geom_point(alpha=0.7) +
  scale_size(range = c(.1, 18), name="Genes") +
  scale_y_discrete(limits = rev) +
  scale_color_gradient(low="blue", high="red") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


# Enhancer Decreasing
sel_rna_transition <- vector()
sel_atac_transition <- vector()
transition_gene <- vector()
transition_OCR <- vector()


for(i in 1:nrow(rna_down)){
  if(rownames(rna_down)[i]%in%CREs_enhancer$gene_ensembl){
    sel_peaks <- CREs_enhancer$ocr[which(CREs_enhancer$gene_ensembl==rownames(rna_down)[i])]
    for(j in 1:length(sel_peaks)){
      ss <- which(atac_dba[sel_peaks[j],transitions_sig]==-1)
      ss_rna <- which(rna_dea[rownames(rna_down)[i],transitions_sig]==-1)
      if(length(ss)==0){next}
      if(length(ss_rna)==0){next}
      if(length(ss)>1){ss <- ss[which(atac_dba[sel_peaks[j],transitions[ss]]==min(atac_dba[sel_peaks[j],transitions[ss]]))]}
      if(length(ss_rna)>1){ss_rna <- ss_rna[which(rna_dea[rownames(rna_down)[i],transitions[ss_rna]]==min(rna_dea[rownames(rna_down)[i],transitions[ss_rna]]))]}
      sel_rna_transition <- c(sel_rna_transition, ss_rna) #the first transition where a sig idecreases of RNA occurs
      sel_atac_transition <- c(sel_atac_transition, ss) #the first transition where a sig decreaes of promoter accessibility occur
      transition_gene <- c(transition_gene,rownames(rna_down)[i])
      transition_OCR <- c(transition_OCR,sel_peaks[j])
    }
  }
}


regulation_down <- data.frame(gene=transition_gene,OCR=transition_OCR,rna=sel_rna_transition,atac=sel_atac_transition)

#Plot
regulation_down_plot <- data.frame(rna=factor(rep(gsub("logFC_","",transitions)[1:max(sel_rna_transition)],max(sel_atac_transition)), levels=gsub("logFC_","",transitions)[1:max(sel_rna_transition)]),
                                   atac=factor(rep(gsub("logFC_","",transitions)[1:max(sel_atac_transition)], each=max(sel_rna_transition)), levels=gsub("logFC_","",transitions)[1:max(sel_atac_transition)]),
                                   values=as.vector(table(regulation_down$rna, regulation_down$atac)))


p_down <- ggplot(regulation_down_plot, aes(x=atac, y=rna, size = values, color=values)) +
  geom_point(alpha=0.7) +
  scale_size(range = c(.1, 18), name="Genes") +
  scale_y_discrete(limits = rev) +
  scale_color_gradient(low="blue", high="red") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#***************************************************




#***************************************************
## Figure 5 C
#***************************************************
## UMAP single-cell annotation and pseudotime
#***************************************************








#***************************************************
## Figure 5 D
#***************************************************
## UMAP single-cell annotation and pseudotime
#***************************************************




#***************************************************
## Figure 5 E
#***************************************************
## Representation of ATAC-Seq signal for CFH promoter region for the 8 cell subpopulations

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
# B cell subpopulation colors: "bcell_colors.rds"

# Auxiliar functions
source("auxiliar_functions.R")


load("bcell_bwplots.RData", verbose = TRUE)
consensus_OCRs_granges <- import("consensus_OCRs_all_cell_subpopulations.bed", format = "BED")
h3k4me1 <- import("H3K4me1_regions.bed", format = "BED")
h3k4me3 <- import("H3K4me3_regions.bed", format = "BED")
h3k27ac <- import("H3K27ac_regions.bed", format = "BED")

# Bigwig plot
#CFH promoter: chr1:196,644,004-196,665,227
gene <- "CFH"
bw_plot2(accDT_by_group, peaks=peak_granges, gene=gene, ylims=c(0,80), start=196644004, end=196665227, chr="chr1",
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

target_gene <- "CFH"
data_to_plot <- data.frame(value=as.numeric(by(as.numeric(rna_norm[rna_annotation$ensembl_gene_id[rna_annotation$hgnc_symbol==target_gene],]), rna_metadata$CellType, median)),
                           class=factor(levels(rna_metadata$CellType),levels = levels(rna_metadata$CellType)))
p <- ggplot(data_to_plot, aes(x=class, y=value, fill=class)) + 
  geom_bar(stat="identity")+
  scale_fill_manual(values=bcell_colors) +
  theme_bw()

#***************************************************




#***************************************************
## Figure 5 F
#***************************************************
## UMAP single-cell annotation and pseudotime
#***************************************************




#***************************************************
## Figure 5 G
#***************************************************
## Representation of ATAC-Seq signal for the cCRE of the CDH7 gene.

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
#CHD7 cCRE: chr8:60,916,193-60,917,014
gene <- "CHD7"
bw_plot2(accDT_by_group, peaks=peak_granges, gene=gene, ylims=c(0,80), start=60916193, end=60917014, chr="chr8",
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

target_gene <- "CHD7"
data_to_plot <- data.frame(value=as.numeric(by(as.numeric(rna_norm[rna_annotation$ensembl_gene_id[rna_annotation$hgnc_symbol==target_gene],]), rna_metadata$CellType, median)),
                           class=factor(levels(rna_metadata$CellType),levels = levels(rna_metadata$CellType)))
p <- ggplot(data_to_plot, aes(x=class, y=value, fill=class)) + 
  geom_bar(stat="identity")+
  scale_fill_manual(values=bcell_colors) +
  theme_bw()

#***************************************************




#***************************************************
## Figure 5 I
#***************************************************
## UMAP single-cell RNA-Seq marker signatures
#***************************************************




#***************************************************
## Figure 5 J
#***************************************************
## UMAP single-cell ATAC-Seq showing signature of promoters from marker genes
#***************************************************





#***************************************************
## Figure 5 K
#***************************************************
## UMAP single-cell ATAC-Seq showing signature of enhancers from marker genes
#***************************************************
