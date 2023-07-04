#***************************************#
# RNA-Seq pipeline - SINGLE-END         #
# script: 05_count_table_exploration.sh #
# script: 05_count_table_exploration.R  #
#***************************************#

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


####################### DATA PREPARATION #######################

### STEP 1: READ DATA

#count matrix
load(count_table)

counts <- df
rownames(counts) <- counts[, 1]
counts[, 1] <- NULL

#remove the first 5 rows with the alignment information
counts <- counts[-c(1:5),] #58,884 x 213

#metadata
metadata <- read.table(metadata,sep="\t", dec=".", header=TRUE, check.names=FALSE) #213 x 23
summary(metadata$CellType)

#select samples of interest
sel <- which(metadata$CellType%in%c("HSC","CLP","ProB","PreB","Immature.B","Transitional.B","Naive.CD5pos","Naive.CD5neg")) #90 samples

metadata2 <- metadata[sel, ] #90 x 23
counts2 <- counts[,as.character(metadata2$Alias)] #58,884 x 90

# We will delete those HSC samples that came from donors without more cell types
sel <- which(metadata2$Donor%in%c("D196","D197","D203","D207","D208","D209","D216","D217","D218","D219")) #9

metadata2 <- metadata2[-sel, ] #80 x 23
counts2 <- counts2[,as.character(metadata2$Alias)] #58,884 x 80

metadata2$CellType <- droplevels(metadata2$CellType)
metadata2$Donor <- droplevels(metadata2$Donor)
# -------------------------------



### STEP 2: LIMMA VOOM NORMALIZATION 
require(edgeR)

#create DGEList object
d0 <- DGEList(counts=counts2)

#filter: remove rows that consistently have zero or very low counts
group <- metadata2$CellType
donor <- metadata2$Donor
design <- model.matrix(~0 + donor + group)

keep <- filterByExpr(d0, design, group=group, min.count = 1) #23,460

d1 <- d0[keep,,keep.lib.sizes=FALSE]

#calculate normalization factors (TMM)
d1_norm <- calcNormFactors(d1)

#voom transformation
v <- voom(d1_norm, design, plot=TRUE)

norm_data <- v$E

####################### EXPLORATORY ANALYSIS #######################

png("mds_voom.png")
plotMDS(v, col = as.numeric(group))
dev.off()

## Our MDS
dist_go <- 1-cor(norm_data, method="spearman")
fit <- cmdscale(dist_go, eig = TRUE, k = 2)
mds <- as.data.frame(fit$points)
mds$CellType <- metadata2$CellType

library(ggplot2)
library(RColorBrewer)

# Labels colors
colors <- c(brewer.pal(n=8, name = 'Dark2'),brewer.pal(n=12, name = 'Set3'),brewer.pal(n=12, name = 'Paired'))

color_types <- colors[1:length(levels(metadata2$CellType))]
names(color_types) <- levels(metadata2$CellType)

p <- ggplot(mds, aes(x = mds[, 1], y = mds[, 2], color = CellType)) +
  geom_point(size = 2) + coord_fixed() + scale_color_manual(values = color_types)

ggsave(plot = p, width = 9, height = 7, dpi = 100, filename = "mds_voom_tunned.pdf")
##There is one strange sample: D213_4
# -------------------------------



### STEP 3: REMOVE OUTLIERS

##Remove samples D213_4 and re-normalized data
sel <- which(metadata2$Alias=="D213_4") #9
metadata3 <- metadata2[-sel,] #79 x 23
counts3 <- counts2[,-sel] #58,884 x 79
# -------------------------------



### STEP 4: GENE ANNOTATION USING biomaRt

library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
biotypes <- getBM(filters="ensembl_gene_id",attributes=c('ensembl_gene_id', 'hgnc_symbol','chromosome_name','start_position','end_position','strand','gene_biotype','percentage_gene_gc_content'), values =rownames(counts3), mart = mart)  #58,778 x 8

#there are two entries for two ensembl gene id
biotypes[biotypes$ensembl_gene_id=="ENSG00000187510",] #keep PLEKHG7
biotypes[biotypes$ensembl_gene_id=="ENSG00000230417",] #keep LINC00595
biotypes[biotypes$ensembl_gene_id=="ENSG00000276085",] #keep CCL3L3

biotypes2 <- biotypes[-which(biotypes$ensembl_gene_id=="ENSG00000187510" & biotypes$hgnc_symbol!="PLEKHG7"),]
biotypes2 <- biotypes2[-which(biotypes2$ensembl_gene_id=="ENSG00000230417" & biotypes2$hgnc_symbol!="LINC00595"),]
biotypes2 <- biotypes2[-which(biotypes2$ensembl_gene_id=="ENSG00000276085" & biotypes2$hgnc_symbol!="CCL3L3"),]
#58,775 x 8

sum(!(rownames(counts3) %in% biotypes2$ensembl_gene_id))
#there are 109 genes from counts table without annotation
#Check what type of ensembl ID are...
rownames(counts3)[!(rownames(counts3) %in% biotypes2$ensembl_gene_id)]
#seems that are out of the actual database

#remove those genes without annotation
counts3_filtered <- counts3[rownames(counts3) %in% biotypes2$ensembl_gene_id,] #58,884 x 79

write.table(biotypes2,"gene_annotation.txt", sep="\t", dec=".", quote=FALSE, row.names=FALSE, col.names=TRUE)
# -------------------------------



### STEP 5: REDO THE ANALYSIS WITHOUT OUTLIERS -> LIMMA VOOM NORMALIZATION 

#create DGEList object
d0 <- DGEList(counts=counts3_filtered)

#filter: remove rows that consistently have zero or very low counts
group <- metadata3$CellType
donor <- metadata3$Donor
design <- model.matrix(~0 + donor + group)

keep <- filterByExpr(d0, design, group=group, min.count = 1) #23,454

d1 <- d0[keep,,keep.lib.sizes=FALSE]

#calculate normalization factors (TMM)
d1_norm <- calcNormFactors(d1)

#voom transformation
v <- voom(d1_norm, design, plot=TRUE)

norm_data <- v$E #23,454 x 79

write.table(norm_data,"VOOM_norm_filtered_data.txt", sep="\t", dec=".", quote=FALSE, row.names=TRUE, col.names=TRUE)


####################### EXPLORATORY ANALYSIS #######################

png("mds_voom_filtered.png")
plotMDS(v, col = as.numeric(group))
dev.off()


## Our MDS
dist_go <- 1-cor(norm_data, method="spearman")
fit <- cmdscale(dist_go, eig = TRUE, k = 2)
mds <- as.data.frame(fit$points)
mds$CellType <- metadata3$CellType

library(ggplot2)
library(RColorBrewer)

# Labels colors
colors <- c(brewer.pal(n=8, name = 'Dark2'),brewer.pal(n=12, name = 'Set3'),brewer.pal(n=12, name = 'Paired'))

color_types <- colors[1:length(levels(metadata3$CellType))]
names(color_types) <- levels(metadata3$CellType)

p <- ggplot(mds, aes(x = mds[, 1], y = mds[, 2], color = CellType)) +
  geom_point(size = 2) + coord_fixed() + scale_color_manual(values = color_types)


##>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Plot for publication
require(BuenColors)
library(directlabels)
library(pals)

#>>>>>>>>>>>> Cell types:
bcell_colors <- jdb_color_map(c("LMPP","CLP","CD8","CD4","B","NK","mono","MEP"))
names(bcell_colors) <- c("HSC","CLP","ProB","PreB","Immature.B","Transitional.B","Naive.CD5pos","Naive.CD5neg")

#order the cell type labels
metadata3$CellType2 <- factor(as.character(metadata3$CellType),levels=c("HSC","CLP","ProB","PreB","Immature.B","Transitional.B","Naive.CD5pos","Naive.CD5neg"))

p3 <- ggplot(mds, aes(x = mds[, 1], y = mds[, 2], color = CellType)) +
  geom_point(size = 2, show.legend=FALSE) + coord_fixed() + scale_color_manual(values = color_types) + directlabels::geom_dl(aes(label = CellType), method = "smart.grid")

ggsave(plot = p, width = 9, height = 7, dpi = 100, filename = "mds_voom_filtered_tunned.pdf")


##Plot the PCA instead of MDS
pca_rna <- prcomp(t(norm_data))
var_prcomp <- pca_rna$sdev^2

# Plot PCs variance
pcvar <- data.frame(var=var_prcomp/sum(var_prcomp), 
                    pc=c(1:length(var_prcomp)))
                    
mds <- data.frame(PC1=pca_rna$x[,1],PC2=pca_rna$x[,2], CellType=metadata3$CellType)

p3 <- ggplot(mds, aes(x = PC1, y = PC2, color = CellType)) +
  geom_point(size = 3, show.legend=FALSE) + coord_fixed(ratio=1.35) + scale_color_manual(values = bcell_colors) + directlabels::geom_dl(aes(label = CellType), method = list("smart.grid", cex=1.5)) + 
  labs(title="PCA: RNA - Log2(Normalized and Filtered Count table (VOOM))",x = paste("PC1 (",round(pcvar[1,1]*100,0),"%)", sep=""), y = paste("PC2 (",round(pcvar[2,1]*100,0),"%)", sep="")) 

ggsave(plot = p3, dpi = 300, filename = "rna_pca_voom_filtered_tunned.pdf")

save(list=c("v","biotypes2","norm_data"), file="voom_to_dea.RData")
# -------------------------------