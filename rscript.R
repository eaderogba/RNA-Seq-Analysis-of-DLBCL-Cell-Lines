# load libraries
library(DESeq2)
library(vsn)
library(dplyr)
library(tidyverse)
library(ggrepel)
library(DEGreport)
library(pheatmap)
library(knitr)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(genefilter)
library(grDevices)

# Create Sample Info
sample_info <- data.frame(
  batch = c("SRR28420795", "SRR28420796", "SRR28420797", "SRR28420798"),
  condition = c("treated", "control", "treated", "control"))

# Set row names to the Accession numbers
rownames(sample_info) <- sample_info$batch
sample_info$batch <- NULL
# Create a data table
write.table(sample_info, "sample_info.txt", sep = "\t", row.names = TRUE)

# load data
meta <- read.table("~/Bionformatics_Projects/TNBC_Lespedeza_RNASeq_Analysis/sample_info.txt", header = TRUE)
feature_count <- read.table("~/Bionformatics_Projects/TNBC_Lespedeza_RNASeq_Analysis/counts/featureCounts_results.mod.txt", header = TRUE, row.names = 1)

# remove the first 5 columns
data <- feature_count[,6:9]

# check to make sure that all rows labels in the counts_data are colums in data
all(colnames(data) == rownames(meta))

# Create dataset and Run DESeq2 analysis
dds <- DESeqDataSetFromMatrix(countData = data, colData = meta, design = ~ condition)
dds <- DESeq(dds)
vsd <- vst(dds)

# Principal component analysis
rld <- rlog(dds, blind = TRUE)
plotPCA(rld, intgroup = "condition") + geom_text(aes(label=name))

# creating contrasts and running a Wald test
contrast <- c("condition", "treated", "control")
res_unshrunken <- results(dds, contrast = contrast)
summary(res_unshrunken)

# Shrinkage of the log2 fold changes
res <- lfcShrink(dds, coef = 2, type = "apeglm", res = res_unshrunken)
summary(res)
head(res)

# Filtering to find significant genes
padj.cutoff <- 0.05 # False Discovery Rate cutoff
significant_results <- res[which(res$padj < padj.cutoff),]

significant_results
file_name = 'significant_padj_0.05.txt'
write.table(significant_results, file_name, quote = FALSE)


# Visualisation
# simple plot for a single gene 
par(mar = c(5, 5, 2, 2))
plotCounts(dds, gene="ENSG00000284616", intgroup = "condition")
 
# heatmap: plot top 5o genes in a heatmap
# Sample distance heatmap
mycols <- brewer.pal(8, "Dark2")[1:length(unique(sample_info$condition))]
sampleDists <- as.matrix(dist(t(assay(rld))))
png("heatmap-samples.png", w=1000, h=1000, pointsize = 20)
heatmap.2(as.matrix(sampleDists), key = F, trace="none",
          col=colorpanel(100, "green", "yellow"),
          ColSideColors=mycols[sample_info$condition], RowSideColors=mycols[sample_info$condition],
          margin=c(15, 15), main="Sample Distance Matrix")
dev.off()

# Principal Component Analysis (PCA)
png("PCA.png", w=1000, h=1000)
DESeq2::plotPCA(rld, intgroup="condition")
dev.off()

# heatmap
# extract the counts from the rlog transformed object
rv <- rowVars(assay(rld))
select <- order(rv, decreasing = T)[seq_len(min(500,length(rv)))]
pc <- prcomp(t(assay(vsd)[select,]))
scores <- data.frame(pc$x, sample_info$condition)
.
.
.
.
.

# Heatmap of data
select <- order(rowMeans(counts(dds, normalized = T)), decreasing = T)[1:1000]

# Check the select object if necessary
print(select)

# Define the palette with a reduced number of colors
my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 100)

# Start the PNG device
png("Heatmap.png", w = 1000, h = 1000, pointsize = 20)

# Generate the heatmap
heatmap.2(assay(vsd)[select,], col = my_palette,
          scale = "row", key = TRUE, keysize = 1, symkey = TRUE,
          density.info = "none", trace = "none",
          cexCol = 0.6, labRow = FALSE,
          main = "1000 Top Expressed Genes Heatmap")

# Close the PNG device
dev.off()








