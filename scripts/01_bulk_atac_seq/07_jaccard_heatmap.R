#***********************************#
# ATAC-Seq pipeline - PAIR-END      #
# script: 07_peaksQC.sh             #
# script: 07_jaccard_heatmap.R      #
#***********************************#

# 7. Quality Control of Peak Calling - Measuring dataset similarity

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


print(workdir)
print(filename)
print(filename2)


#Setting working directory
setwd(workdir)

library(circlize)
library(ComplexHeatmap)
library(RColorBrewer)


###########################################
# PLOT JACCARD SCORE - SIMILARITY MEASURE #
###########################################

# Jaccard
jaccard <- read.table(filename, header=TRUE, sep=" ", row.names=1, check.names=FALSE)

# Metadata
metadata <- read.table(filename2, sep=";", header=TRUE)
rownames(metadata) <- metadata$SampleID
metadata <- metadata[rownames(jaccard),]


# List of colors
# https://www.stat.ubc.ca/~jenny/STAT545A/block14_colors.html

# Labels colors
colors <- c(brewer.pal(n=8, name = 'Dark2'),brewer.pal(n=12, name = 'Set3'),brewer.pal(n=12, name = 'Paired'))

color_types <- colors[1:length(levels(metadata$Tissue))]
names(color_types) <- levels(metadata$Tissue)

color_case <- colors[1:length(levels(metadata$Factor))]
names(color_case) <- levels(metadata$Factor)

color_batch <- colors[1:length(levels(metadata$Batch))]
names(color_batch) <- levels(metadata$Batch)

ha = HeatmapAnnotation(
    Cell_type = metadata$Tissue,
    Case = metadata$Factor,  
    Batch = metadata$Batch,  
    col = list(Cell_type = color_types,
    Case = color_case,
    Batch = color_batch),
    gp = gpar(col = "black"),  #black border per cell
    annotation_legend_param = list(Cell_type = list(border="black"), Case = list(border="black"), Batch = list(border="black"))
)


# Heatmap colors
col_heatmap = colorRamp2(seq(0, 1, length = 3), c("white", "gray", "red"), space = "RGB")


nr <- nrow(jaccard)
png("jaccard.png", width = 0.15*nr + 1.3150, height = 0.15*nr + 1.3150, units = "in", res = 100)
  draw(Heatmap(as.matrix(jaccard), top_annotation = ha, show_row_names = FALSE, col=col_heatmap, column_names_gp = gpar(fontsize = 4), border=TRUE, heatmap_legend_param = list(title = "jaccard", border="black"), height = unit(3, "mm")*nr))
dev.off()


# MDS
library(ggplot2)

#Remove samples with missings
sel <- which(colSums(is.na(jaccard))!=0)
if(length(sel)!=0){
  jaccard <- jaccard[-sel,-sel]
  metadata <- metadata[-sel,]
}


mds <- cmdscale(1-jaccard,eig=TRUE, k=2)

mplot <- data.frame(MDS1=mds$points[,1], MDS2=mds$points[,2], Cell_type=metadata$Tissue, Batch=metadata$Batch)


p <- ggplot(mplot, aes(x=MDS1, y=MDS2, color = Cell_type, label=rownames(jaccard))) + 
  geom_point(size=2) +
  scale_color_manual(values=color_types) + 
  labs(title="MDS - Measure of similarity (Jaccard)", x ="Coordinate 1", y="Coordinate 2")

ggsave(plot = p, width = 0.06*nr, height = 0.04*nr, dpi = 100, filename = "jaccard_mds_tissue.pdf")


p <- ggplot(mplot, aes(x=MDS1, y=MDS2, color = Batch, label=rownames(jaccard))) + 
  geom_point(size=2) +
  scale_color_manual(values=color_batch) + 
  labs(title="MDS - Measure of similarity (Jaccard)", x ="Coordinate 1", y="Coordinate 2")

ggsave(plot = p, width = 0.06*nr, height = 0.04*nr, dpi = 100, filename = "jaccard_mds_batch.pdf")


p <- ggplot(mplot, aes(x=MDS1, y=MDS2, color = Cell_type)) + 
  geom_point(size=0.5) +
  geom_text(size=2, aes(label=rownames(jaccard))) +
  scale_color_manual(values=color_types) + 
  labs(title="MDS - Measure of similarity (Jaccard)", x ="Coordinate 1", y="Coordinate 2")

ggsave(plot = p, width = 0.06*nr, height = 0.04*nr, dpi = 100, filename = "jaccard_mds_annotated.pdf")
