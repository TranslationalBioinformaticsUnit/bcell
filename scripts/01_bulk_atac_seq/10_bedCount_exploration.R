#*************************************#
# ATAC-Seq pipeline - PAIR-END        #
# script: 10_bedCount_exploration.sh  #
# script: 10_bedCount_exploration.R   #
#*************************************#

# 10. Count the number of reads mapping to each feature - QUALITY CONTROL

##First read in the arguments listed at the command line
args=(commandArgs(TRUE))

##args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.
if(length(args)==0){
    stop("No arguments supplied.")
}else{
    for(i in 1:length(args)){
         eval(parse(text=args[[i]]))
    }
}


# Setting working directory
setwd(workdir)

# Library required
library(Rtsne)
library(edgeR)
library(ggplot2)
library(umap)


set.seed(12345)

## Count table
d <- read.table(filename)
rownames(d) <- paste(d[,1],d[,2],d[,3],sep="_")
counts <- d[,-c(1:6)]

#Name samples
colnames(counts)<- gsub(".clean.bam","",grep("clean.bam$",dir(workdir_bams), value=TRUE))

# Metadata
metadata <- read.table(filename2, sep=";", header=TRUE)
rownames(metadata) <- metadata$SampleID
metadata <- metadata[colnames(counts),]


#remove samples from Memory cell type, excluded from the analysis
sel_memory <- which(metadata$Tissue=="Memory_Class_Switch") #8
counts <- counts[,-sel_memory]
 
metadata <- metadata[colnames(counts),]
metadata$Tissue <- droplevels(metadata$Tissue)

write.table(counts, "count_table_with_colnames.txt", sep="\t", dec=".", quote=FALSE)



# Relevant covariables
Cell <- factor(metadata$Tissue)
Donor <- factor(metadata$Factor)

### Differential binding using edgeR

# We are interested in determine differential binding in B-cell differentiation
y <- DGEList(counts=counts)

## Filtering
# low counts across all libraries provide little evidence for differential binding
# criteria: count of 5-10 in a library to be considered as a peak in that library
# Roughly speaking, the strategy keeps genes that have at least ‘min.count’ reads in a
# worthwhile number samples. More precisely, the filtering keeps genes that have count-per-million (CPM) above _k_ in _n_ samples, where _k_ is determined by ‘min.count’ and by the sample library sizes and _n_ is determined by the design matrix.
# _n_ is essentially the smallest group sample size or, more precisely, the minimum inverse leverage for any fitted value. If all the group sizes are large, then this is relaxed slightly, but with _n_ always greater than 70% of the smallest group size.

keep <- filterByExpr(y, group=Cell, min.count=10, min.total.count=15) 
table(keep)
# only 151 peaks that will be filtered out 

# keep de peaks that pass the filter and recalculate the library size.
y <- y[keep, , keep.lib.sizes=FALSE]


## Normalization
# The most obvious technical factor that affects the read counts, other than gene expression
# levels, is the sequencing depth of each RNA sample. edgeR adjusts any differential expression
# analysis for varying sequencing depths as represented by differing library sizes.

# Following Buenrostor paper, data is normalized using the TMM function in EdgeR.
y <- calcNormFactors(y, method="TMM")
tmm_table <- cpm(y)
write.table(tmm_table, "tmm_filtered_data.txt", sep="\t", dec=".", quote=FALSE)


### Exploration of count/TMM data matrix

# List of colors
# https://www.stat.ubc.ca/~jenny/STAT545A/block14_colors.html
library(RColorBrewer)

# Labels colors
library(pals)
colors <- tol()

#order the cell type labels
metadata$Tissue2 <- factor(as.character(metadata$Tissue),levels=c("HSC","CLP","proB","preB","ImmatureB","Transitional_B","Naive_CD5pos","Naive_CD5neg"))

color_types <- colors[c(1,2,4,5,6,8,9,12)]
names(color_types) <- levels(metadata$Tissue2)


# Transformation
lcounts <- log2(counts+1)
ltmm <- log2(tmm_table+1)

# PCA
raw_pca <- prcomp(t(lcounts)) 
tmm_pca <- prcomp(t(ltmm))

pr <- summary(raw_pca)$importance[,1:5]
pt <- summary(tmm_pca)$importance[,1:5]

# TSNE
raw_tsne <- Rtsne(t(lcounts), perplexity=20)
tmm_tsne <- Rtsne(t(ltmm), perplexity=20)

# UMAP
raw_umap <- umap(t(lcounts))
tmm_umap <- umap(t(ltmm))



mplot <- data.frame(PC1=raw_pca$x[,1], PC2=raw_pca$x[,2], TSNE1=raw_tsne$Y[,1], TSNE2=raw_tsne$Y[,2], umap1=raw_umap$layout[,1], umap2=raw_umap$layout[,2], Cell_type=metadata$Tissue2)

p1 <- ggplot(mplot, aes(x=PC1, y=PC2, color = Cell_type, label=rownames(mplot))) + 
  geom_point(size=3) +
  scale_color_manual(values=color_types) + 
  labs(title="PCA - RAW Count table",x = paste("PC1 (",round(pr[2,1]*100,0),"%)", sep=""), y = paste("PC2 (",round(pr[2,2]*100,0),"%)", sep=""))

p2 <- ggplot(mplot, aes(x=TSNE1, y=TSNE2, color = Cell_type, label=rownames(mplot))) + 
  geom_point(size=3) +
  scale_color_manual(values=color_types) + 
  labs(title="TSNE - RAW Count table")

p3 <- ggplot(mplot, aes(x=umap1, y=umap2, color = Cell_type)) + 
  geom_point(size=3) +
  scale_color_manual(values=color_types) + 
  labs(title="UMAP - RAW Count table")
  
ggsave(plot = p1, width = 9, height = 7, dpi = 100, filename = "raw_count_table_pca.pdf")
ggsave(plot = p2, width = 9, height = 7, dpi = 100, filename = "raw_count_table_tsne.pdf")
ggsave(plot = p3, width = 9, height = 7, dpi = 100, filename = "raw_count_table_umap.pdf")


ggsave(plot = p1+geom_text(size=2, aes(label=rownames(mplot))), width = 9, height = 7, dpi = 100, filename = "raw_count_table_pca_annotated.pdf")
ggsave(plot = p2+geom_text(size=2, aes(label=rownames(mplot))), width = 9, height = 7, dpi = 100, filename = "raw_count_table_tsne_annotated.pdf")
ggsave(plot = p3+geom_text(size=2, aes(label=rownames(mplot))), width = 9, height = 7, dpi = 100, filename = "raw_count_table_umap_annotated.pdf")


library(directlabels)

mplot <- data.frame(PC1=tmm_pca$x[,1], PC2=tmm_pca$x[,2], TSNE1=tmm_tsne$Y[,1], TSNE2=tmm_tsne$Y[,2], umap1=tmm_umap$layout[,1], umap2=tmm_umap$layout[,2], Cell_type=metadata$Tissue)

p1 <- ggplot(mplot, aes(x=PC1, y=PC2, color = Cell_type, label=rownames(mplot))) + 
  geom_point(size=3) +
  scale_color_manual(values=color_types) + 
  labs(title="PCA: ATAC - Normalized Count table (TMM)",x = paste("PC1 (",round(pt[2,1]*100,0),"%)", sep=""), y = paste("PC2 (",round(pt[2,2]*100,0),"%)", sep=""))

p2 <- ggplot(mplot, aes(x=TSNE1, y=TSNE2, color = Cell_type)) + 
  geom_point(size=3) +
  scale_color_manual(values=color_types) + 
  labs(title="TSNE - Normalized Count table (TMM)")

p3 <- ggplot(mplot, aes(x=umap1, y=umap2, color = Cell_type)) + 
  geom_point(size=3) +
  scale_color_manual(values=color_types) + 
  labs(title="UMAP - Normalized Count table")

  
ggsave(plot = p1, width = 9, height = 7, dpi = 100, filename = "tmm_count_table_pca.pdf")
ggsave(plot = p2, width = 9, height = 7, dpi = 100, filename = "tmm_count_table_tsne.pdf")
ggsave(plot = p3, width = 9, height = 7, dpi = 100, filename = "tmm_count_table_umap.pdf")


ggsave(plot = p1+geom_text(size=2, aes(label=rownames(mplot))), width = 9, height = 7, dpi = 100, filename = "tmm_count_table_pca_annotated.pdf")
ggsave(plot = p2+geom_text(size=2, aes(label=rownames(mplot))), width = 9, height = 7, dpi = 100, filename = "tmm_count_table_tsne_annotated.pdf")
ggsave(plot = p3+geom_text(size=2, aes(label=rownames(mplot))), width = 9, height = 7, dpi = 100, filename = "tmm_count_table_umap_annotated.pdf")



######### There are any outlier or strange sample within our data???

#If yes, explude identified samples and re-normalized and explore the data again.


## SAMPLE TO REMOVE:
sample_to_remove <- "ImmatureB_ATAC_D215_S5"

#remove outlier samples from the analysis
counts <- counts[,!colnames(counts)%in%sample_to_remove]
metadata <- metadata[colnames(counts),]


# Relevant covariables
Cell <- factor(metadata$Tissue)
Donor <- factor(metadata$Factor)

### Differential binding using edgeR

# We are interested in determine differential binding in B-cell differentiation
y <- DGEList(counts=counts)

## Filtering
# low counts across all libraries provide little evidence for differential binding
# criteria: count of 5-10 in a library to be considered as a peak in that library
# Roughly speaking, the strategy keeps genes that have at least ‘min.count’ reads in a
# worthwhile number samples. More precisely, the filtering keeps genes that have count-per-million (CPM) above _k_ in _n_ samples, where _k_ is determined by ‘min.count’ and by the sample library sizes and _n_ is determined by the design matrix.
# _n_ is essentially the smallest group sample size or, more precisely, the minimum inverse leverage for any fitted value. If all the group sizes are large, then this is relaxed slightly, but with _n_ always greater than 70% of the smallest group size.

keep <- filterByExpr(y, group=Cell, min.count=10, min.total.count=15) 
table(keep)
# only 151 peaks that will be filtered out 

# keep de peaks that pass the filter and recalculate the library size.
y <- y[keep, , keep.lib.sizes=FALSE]


## Normalization
# The most obvious technical factor that affects the read counts, other than gene expression
# levels, is the sequencing depth of each RNA sample. edgeR adjusts any differential expression
# analysis for varying sequencing depths as represented by differing library sizes.

# Following Buenrostor paper, data is normalized using the TMM function in EdgeR.
y <- calcNormFactors(y, method="TMM")
tmm_table <- cpm(y)
write.table(tmm_table, "tmm_filtered_data.txt", sep="\t", dec=".", quote=FALSE)


### Exploration of count/TMM data matrix

color_types <- colors[1:length(levels(metadata$Tissue))]
names(color_types) <- levels(metadata$Tissue)


# Transformation
lcounts <- log2(counts+1)
ltmm <- log2(tmm_table+1)

# PCA
tmm_pca <- prcomp(t(ltmm))
pt <- summary(tmm_pca)$importance[,1:5]

# TSNE
tmm_tsne <- Rtsne(t(ltmm), perplexity=20)

# UMAP
tmm_umap <- umap(t(ltmm))



mplot <- data.frame(PC1=tmm_pca$x[,1], PC2=tmm_pca$x[,2], TSNE1=tmm_tsne$Y[,1], TSNE2=tmm_tsne$Y[,2], umap1=tmm_umap$layout[,1], umap2=tmm_umap$layout[,2], Cell_type=metadata$Tissue)

p1 <- ggplot(mplot, aes(x=PC1, y=PC2, color = Cell_type, label=rownames(mplot))) + 
  geom_point(size=3) +
  scale_color_manual(values=color_types) + 
  labs(title="PCA - Normalized and Filtered Count table (TMM)",x = paste("PC1 (",round(pt[2,1]*100,0),"%)", sep=""), y = paste("PC2 (",round(pt[2,2]*100,0),"%)", sep=""))

p2 <- ggplot(mplot, aes(x=TSNE1, y=TSNE2, color = Cell_type)) + 
  geom_point(size=3) +
  scale_color_manual(values=color_types) + 
  labs(title="TSNE - Normalized and Filtered Count table(TMM)")

p3 <- ggplot(mplot, aes(x=umap1, y=umap2, color = Cell_type)) + 
  geom_point(size=3) +
  scale_color_manual(values=color_types) + 
  labs(title="UMAP - Normalized and Filtered Count table")

  
ggsave(plot = p1, width = 9, height = 7, dpi = 100, filename = "tmm_filtered_count_table_pca.pdf")
ggsave(plot = p2, width = 9, height = 7, dpi = 100, filename = "tmm_filtered_count_table_tsne.pdf")
ggsave(plot = p3, width = 9, height = 7, dpi = 100, filename = "tmm_filtered_count_table_umap.pdf")


ggsave(plot = p1+geom_text(size=2, aes(label=rownames(mplot))), width = 9, height = 7, dpi = 100, filename = "tmm_filtered_count_table_pca_annotated.pdf")
ggsave(plot = p2+geom_text(size=2, aes(label=rownames(mplot))), width = 9, height = 7, dpi = 100, filename = "tmm_filtered_count_table_tsne_annotated.pdf")
ggsave(plot = p3+geom_text(size=2, aes(label=rownames(mplot))), width = 9, height = 7, dpi = 100, filename = "tmm_filtered_count_table_umap_annotated.pdf")


##plot for publication
require(BuenColors)

#>>>>>>>>>>>> Cell types:
bcell_colors <- jdb_color_map(c("MPP","LMPP","CD8","CD4","GMP-B","Ery","mono","MEP"))

names(bcell_colors) <- c("HSC","CLP","proB","preB","ImmatureB","Transitional_B","Naive_CD5pos","Naive_CD5neg")


library(directlabels)
library(pals)

mds <- data.frame(PC1=tmm_pca$x[,1],PC2=tmm_pca$x[,2], CellType=atac_metadata$Tissue2)

p3 <- ggplot(mds, aes(x = PC1, y = PC2, color = CellType)) +
  geom_point(size = 3, show.legend=FALSE) + coord_fixed(ratio=1.35) + scale_color_manual(values = bcell_colors) + directlabels::geom_dl(aes(label = CellType), method = list("smart.grid", cex=1.5)) + 
  labs(title="PCA: ATAC - Normalized and Filtered Count table (TMM)",x = paste("PC1 (",round(pt[2,1]*100,0),"%)", sep=""), y = paste("PC2 (",round(pt[2,2]*100,0),"%)", sep=""))
  
ggsave(plot = p3, dpi = 300, filename = "atac_pca_tmm_filtered_tunned.pdf")



##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Correlation plot
##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
atac <- read.table("tmm_filtered_data.txt",sep="\t",dec=".",header=TRUE, check.names=FALSE)

atac_metadata <- read.csv("sampleSheet_diffbind_atacseq_bcells_only_merged.csv", sep=";")
atac_metadata <- atac_metadata[atac_metadata$SampleID%in%colnames(atac),]
atac <- atac[,as.character(atac_metadata$SampleID)]

M <- cor(log2(atac+1))

atac_metadata$Tissue2 <- factor(as.character(atac_metadata$Tissue),levels=c("HSC","CLP","proB","preB","ImmatureB","Transitional_B","Naive_CD5pos","Naive_CD5neg"))

names(bcell_colors) <- levels(atac_metadata$Tissue2)

ha = HeatmapAnnotation(stage = atac_metadata$Tissue2,
                       col = list(stage =bcell_colors))

row_ha = rowAnnotation(stage = atac_metadata$Tissue2,
                       col = list(stage =bcell_colors))


pdf("peak_correlation.pdf")
Heatmap(M, top_annotation = ha, show_row_names = FALSE, show_column_names =FALSE, column_names_gp = gpar(fontsize = 6), col=f1, 
        width = unit(8, "cm"), height = unit(8, "cm"), name="rho",
        column_title = "ATAC: Pearson Correlation")

Heatmap(M, top_annotation = ha, left_annotation=row_ha, show_row_names = FALSE, show_column_names =FALSE, column_names_gp = gpar(fontsize = 6), col=f1, 
        width = unit(8, "cm"), height = unit(8, "cm"), name="rho",
        column_title = "ATAC: Pearson Correlation")
dev.off()

pdf("peak_correlation_ordered.pdf")
Heatmap(M, top_annotation = ha, show_row_names = FALSE, show_column_names =FALSE, column_names_gp = gpar(fontsize = 6), col=f1, 
        width = unit(8, "cm"), height = unit(8, "cm"), name="rho",
        column_title = "ATAC: Pearson Correlation",
        column_order = order(atac_metadata$Tissue2),
        row_order = order(atac_metadata$Tissue2))

Heatmap(M, top_annotation = ha, left_annotation=row_ha, show_row_names = FALSE, show_column_names =FALSE, column_names_gp = gpar(fontsize = 6), col=f1, 
        width = unit(8, "cm"), height = unit(8, "cm"), name="rho",
        column_title = "ATAC: Pearson Correlation",
        column_order = order(atac_metadata$Tissue2),
        row_order = order(atac_metadata$Tissue2))
dev.off()


##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Correlation differentiating by peak type
##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# !! This part is run after peak annotation !! (script: 11_peak_annotation.sh)
#Annotation based on chipseeker 1Kb promoter

atac_anno <- read.table("consensus_peaks_annotation_chipseeker_1Mb_promoter.txt", header=TRUE, sep=",")
rownames(atac_anno) <- paste(atac_anno[,2],atac_anno[,3],atac_anno[,4],sep="_")
atac_anno2 <- atac_anno[rownames(atac),]

##>>>> PROMOTER
atac_promoter <- atac[atac_anno2$annotation2=="Promoter",]
#14589 x 78
atac_promoter <- atac_promoter[,as.character(atac_metadata$SampleID)]

M <- cor(log2(atac_promoter+1))
M_pro <- M

pdf("promoter_peak_correlation.pdf")
Heatmap(M, top_annotation = ha, show_row_names = FALSE, show_column_names =FALSE, column_names_gp = gpar(fontsize = 6), col=f1, 
        width = unit(8, "cm"), height = unit(8, "cm"), name="rho",
        column_title = "ATAC-Promoters: Pearson Correlation")

Heatmap(M, top_annotation = ha, left_annotation=row_ha, show_row_names = FALSE, show_column_names =FALSE, column_names_gp = gpar(fontsize = 6), col=f1, 
        width = unit(8, "cm"), height = unit(8, "cm"), name="rho",
        column_title = "ATAC-Promoters: Pearson Correlation")
dev.off()

pdf("promoter_peak_correlation_ordered.pdf")
Heatmap(M, top_annotation = ha, show_row_names = FALSE, show_column_names =FALSE, column_names_gp = gpar(fontsize = 6), col=f1, 
        width = unit(8, "cm"), height = unit(8, "cm"), name="rho",
        column_title = "ATAC-Promoters: Pearson Correlation",
        column_order = order(atac_metadata$Tissue2),
        row_order = order(atac_metadata$Tissue2))

Heatmap(M, top_annotation = ha, left_annotation=row_ha, show_row_names = FALSE, show_column_names =FALSE, column_names_gp = gpar(fontsize = 6), col=f1, 
        width = unit(8, "cm"), height = unit(8, "cm"), name="rho",
        column_title = "ATAC-Promoters: Pearson Correlation",
        column_order = order(atac_metadata$Tissue2),
        row_order = order(atac_metadata$Tissue2))
dev.off()



##>>>> INTRONIC
atac_intron <- atac[atac_anno2$annotation2=="Intron",]
#44,298 x 78
atac_intron <- atac_intron[,as.character(atac_metadata$SampleID)]

M <- cor(log2(atac_intron+1))
M_intro <- M

pdf("intron_peak_correlation.pdf")
Heatmap(M, top_annotation = ha, show_row_names = FALSE, show_column_names =FALSE, column_names_gp = gpar(fontsize = 6), col=f1, 
        width = unit(8, "cm"), height = unit(8, "cm"), name="rho",
        column_title = "ATAC-Intronic: Pearson Correlation")

Heatmap(M, top_annotation = ha, left_annotation=row_ha, show_row_names = FALSE, show_column_names =FALSE, column_names_gp = gpar(fontsize = 6), col=f1, 
        width = unit(8, "cm"), height = unit(8, "cm"), name="rho",
        column_title = "ATAC-Intronic: Pearson Correlation")
dev.off()

pdf("intron_peak_correlation_ordered.pdf")
Heatmap(M, top_annotation = ha, show_row_names = FALSE, show_column_names =FALSE, column_names_gp = gpar(fontsize = 6), col=f1, 
        width = unit(8, "cm"), height = unit(8, "cm"), name="rho",
        column_title = "ATAC-Intronic: Pearson Correlation",
        column_order = order(atac_metadata$Tissue2),
        row_order = order(atac_metadata$Tissue2))

Heatmap(M, top_annotation = ha, left_annotation=row_ha, show_row_names = FALSE, show_column_names =FALSE, column_names_gp = gpar(fontsize = 6), col=f1, 
        width = unit(8, "cm"), height = unit(8, "cm"), name="rho",
        column_title = "ATAC-Intronic: Pearson Correlation",
        column_order = order(atac_metadata$Tissue2),
        row_order = order(atac_metadata$Tissue2))
dev.off()



##>>>> INTERGENIC
atac_de <- atac[atac_anno2$annotation2=="Distal Intergenic",]
#31446 x 78
atac_de <- atac_de[,as.character(atac_metadata$SampleID)]

M <- cor(log2(atac_de+1))
M_inter <- M

pdf("intergenic_peak_correlation.pdf")
Heatmap(M, top_annotation = ha, show_row_names = FALSE, show_column_names =FALSE, column_names_gp = gpar(fontsize = 6), col=f1, 
        width = unit(8, "cm"), height = unit(8, "cm"), name="rho",
        column_title = "ATAC-Intergenic: Pearson Correlation")

Heatmap(M, top_annotation = ha, left_annotation=row_ha, show_row_names = FALSE, show_column_names =FALSE, column_names_gp = gpar(fontsize = 6), col=f1, 
        width = unit(8, "cm"), height = unit(8, "cm"), name="rho",
        column_title = "ATAC-Intergenic: Pearson Correlation")
dev.off()

pdf("intergenic_peak_correlation_ordered.pdf")
Heatmap(M, top_annotation = ha, show_row_names = FALSE, show_column_names =FALSE, column_names_gp = gpar(fontsize = 6), col=f1, 
        width = unit(8, "cm"), height = unit(8, "cm"), name="rho",
        column_title = "ATAC-Intergenic: Pearson Correlation",
        column_order = order(atac_metadata$Tissue2),
        row_order = order(atac_metadata$Tissue2))

Heatmap(M, top_annotation = ha, left_annotation=row_ha, show_row_names = FALSE, show_column_names =FALSE, column_names_gp = gpar(fontsize = 6), col=f1, 
        width = unit(8, "cm"), height = unit(8, "cm"), name="rho",
        column_title = "ATAC-Intergenic: Pearson Correlation",
        column_order = order(atac_metadata$Tissue2),
        row_order = order(atac_metadata$Tissue2))
dev.off()




##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
## Exploration the peak type distribution by T-SNE and relation with Gini index
##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#all peaks
atac <- read.table("tmm_filtered_data.txt",sep="\t",dec=".",header=TRUE, check.names=FALSE)

atac_metadata <- read.csv("sampleSheet_diffbind_atacseq_bcells_only_merged.csv", sep=";")
atac_metadata <- atac_metadata[atac_metadata$SampleID%in%colnames(atac),]
atac <- atac[,as.character(atac_metadata$SampleID)]

anno <- read.table("consensus_peaks_annotation_chipseeker_1Mb_promoter.txt", header=TRUE, sep=",")
rownames(atac_anno) <- paste(atac_anno[,2],atac_anno[,3],atac_anno[,4],sep="_")


#UMAP
library(umap)
load("umap_atac.RData")
##set.seed(64756)
#set.seed(123)
#atac.umap <- umap(log2(atac+1))

require(ggplot2)

##Density UMAP

#Promoters
promoter <- rep(0, nrow(atac_anno2))
promoter[atac_anno2$annotation2=="Promoter"] <- 1

m <- data.frame(umap1=atac.umap$layout[,1],umap2=atac.umap$layout[,2],Promoter=as.factor(promoter))

require(BuenColors)

gg <- ggplot(m[m$Promoter==1,], aes(x=umap1, y=umap2)) + 
  stat_density_2d(aes(fill=..density..),geom="raster",contour=FALSE, n = 500, h=1) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  scale_fill_gradientn(colors = jdb_palette("horizon")) + 
  labs(subtitle="ATAC-Seq peaks - Promoters", 
       y="UMAP2", 
       x="UMAP1", 
       title="UMAP") +
  guides(colour = guide_legend(override.aes = list(size=2)))

gg

ggsave("umap_promoters.pdf")


#Intronic
promoter <- rep(0, nrow(atac_anno2))
promoter[atac_anno2$annotation2=="Intron"] <- 1

m <- data.frame(umap1=atac.umap$layout[,1],umap2=atac.umap$layout[,2],Intron=as.factor(promoter))

gg <- ggplot(m[m$Intron==1,], aes(x=umap1, y=umap2)) + 
  stat_density_2d(aes(fill=..density..),geom="raster",contour=FALSE, n = 500, h=1) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  scale_fill_gradientn(colors = jdb_palette("horizon")) + 
  labs(subtitle="ATAC-Seq peaks - Introns", 
       y="UMAP2", 
       x="UMAP1", 
       title="UMAP") +
  guides(colour = guide_legend(override.aes = list(size=2)))

gg

ggsave("umap_intron.pdf")


#Distal intergenic
promoter <- rep(0, nrow(atac_anno2))
promoter[atac_anno2$annotation2=="Distal Intergenic"] <- 1

m <- data.frame(umap1=atac.umap$layout[,1],umap2=atac.umap$layout[,2],Intergenic=as.factor(promoter))
gg <- ggplot(m[m$Intergenic==1,], aes(x=umap1, y=umap2)) + 
  stat_density_2d(aes(fill=..density..),geom="raster",contour=FALSE, n = 500, h=1) + 
  scale_x_continuous(expand = c(0, 0)) + 
  scale_y_continuous(expand = c(0, 0)) + 
  scale_fill_gradientn(colors = jdb_palette("horizon")) + 
  labs(subtitle="ATAC-Seq peaks - Distal intergenic", 
       y="UMAP2", 
       x="UMAP1", 
       title="UMAP") +
  guides(colour = guide_legend(override.aes = list(size=2)))

gg

ggsave("umap_intergenic.pdf")



###Plot by Gini Index
require(REAT)

gini_output <- vector()
for(i in 1:nrow(atac)){
  r <- gini(atac[i,]) 
  gini_output <- c(gini_output,r)
}

m <- data.frame(umap1=atac.umap$layout[,1],umap2=atac.umap$layout[,2],gini_index=gini_output)

gg <- ggplot(m, aes(x=umap1, y=umap2, color=gini_index)) + 
  geom_point(size=0.1) + 
  #  scale_color_gradientn(colours = rainbow(5)) +
  scale_color_gradientn(colours = rainbow(3)) +
  labs(subtitle="ATAC-Seq peaks", 
       y="UMAP2", 
       x="UMAP1", 
       title="UMAP") + 
  theme_void()

ggsave("umap_gini_index.pdf")
