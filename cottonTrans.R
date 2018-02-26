
## cotton.R
## Author: Parampreet Sodhi, Mohammad Najaf-Panah, Jesse Sombrano
## Created: May 3, 2017
## RNA-seq analysis with DESeq2 for three cotton genotypes in response 
## to salt over time.

# Import & pre-process ----------------------------------------------------
require("DESeq2")
require("RColorBrewer")
require("gplots")
require("ggplot2")


# Import data from transcriptcount.csv
# Previously ran at command line something like this:
# ./prepDE.py -i /home/tsunami1/scratch/509/mpanah/ -g /home/tsunami1/scratch/509/mpanah/genecount.csv 
# -t /home/tsunami1/scratch/509/mpanah/transcriptcount.csv
countdata <- read.csv("/home/tsunami1/scratch/509/mpanah/transcriptcount.csv", header=TRUE)



# Assign condition (first two are controls, second four contain the treated by salt and so on)
condition <- factor(c(rep("Control", 2), rep("NaCl", 4), rep("Control", 2), rep("NaCl", 4), rep("Control", 2), rep("NaCl", 4)))
response <- factor(c(rep("Salt Tolerant",12), rep("Salt Sensitive",6)))
# Analysis with DESeq2 ----------------------------------------------------

dds <- DESeqDataSetFromMatrix(countData=countdata, colData= DataFrame(condition, response), design = ~condition + response + condition:response, tidy = TRUE)
dds

# Run the DESeq pipeline
dds <- DESeq(dds)
resultsNames(dds)
res <- results(dds, contrast = list("conditionNaCl.responseSalt.Tolerant"))

resOrdered <- res[order(res$padj),]
head(resOrdered)

summary(res)

plotMA(res, main="DESeq2",
       ylim=c(-2,2))

sel <- which(abs(res$log2FoldChange) > 2 & res$padj < 1e-7)

#o <- order(res$padj)

for(i in sel) {
  plotCounts(
    dds, gene=i, #o[i],
    intgroup=c("response", "condition"))
  legend("left", c(i, paste(res$log2FoldChange[i]), paste(res$padj[i])))
}

plotCounts(
  dds, gene=which.min(res$padj),
  intgroup=c("response", "condition"))

write.csv(
  as.data.frame(resOrdered), 
  file="/home/tsunami1/scratch/509/mpanah/condition_treated_gene_results.csv"
)

d <- plotCounts(
  dds, gene=which.min(res$padj), 
  intgroup=c("condition", "response"), 
  returnData=TRUE)

ggplot(d, aes(x=condition, y=count)) +
  geom_point(position=position_jitter(w=0.1, h=0)) +
  scale_y_log10(breaks=c(25,100,400))

ggplot(d1, aes(x=response ,y=count)) +
  geom_point(position=position_jitter(w=0.1, h=0)) +
  scale_y_log10(breaks=c(25,100,400))

select <- order(
  rowMeans(counts(dds,normalized=TRUE)), 
  decreasing=TRUE
)[1:30]

hmcol <- colorRampPalette(
  brewer.pal(9, "GnBu"))(100)

heatmap.2(
  counts(dds,normalized=TRUE)[select,], 
  col = hmcol, Rowv = FALSE, 
  Colv = FALSE, scale="none",
  dendrogram="none", trace="none", 
  margin=c(10,6)
)

rld <- rlog(dds) # regularized log transformation
rld1 <- rlog(dds) # regularized log transformation

heatmap.2(
  assay(rld)[select,], col = hmcol, 
  Rowv = FALSE, Colv = FALSE, scale="none", 
  dendrogram="none", trace="none", 
  margin=c(10, 6)
)

vsd <- varianceStabilizingTransformation(dds)

heatmap.2(
  assay(vsd)[select,], col = hmcol,
  Rowv = FALSE, Colv = FALSE, scale="none",
  dendrogram="none", trace="none", margin=c(10, 6)
)

rlogMat <- assay(rld)
vstMat <- assay(vsd)

distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- 
  colnames(mat) <- 
  with(colData(dds), 
       paste(condition, sep=" : "))

distsRL <- dist(t(assay(rld)))
mat <- as.matrix(distsRL)
rownames(mat) <- 
  colnames(mat) <- 
  with(colData(dds), 
       paste(response, sep=" : "))

heatmap.2(mat, trace="none", 
          col = rev(hmcol), 
          margin=c(13, 13))


plotPCA(
  rld, 
  intgroup=c("condition")
)

plotPCA(
  rld1, 
  intgroup=c("response")
)

data <- plotPCA(
  rld, intgroup=c("condition"), 
  returnData=TRUE)

data1 <- plotPCA(
  rld1, intgroup=c("response"), 
  returnData=TRUE)

percentVar <- round(
  100 * attr(data, "percentVar"))

percentVar1 <- round(
  100 * attr(data1, "percentVar"))

ggplot(
  data, 
  aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))



 ## volcano plot  ######################################################

 library("ggrepel") #Avoid overlapping labels
 library("plyr")

 df <- as.data.frame(res)

 #Will have different colors depending on significance:
 mutateddf <- mutate(
  df, sig=ifelse(
    padj<=1e-7 & abs(log2FoldChange) >=2, 
    "padj≤0.05, |lfc|≥1", "Insignificant")) 

 #convert the rownames to a column
 input <- cbind(gene=rownames(mutateddf ), mutateddf ) 
 volc <- ggplot(input, aes(log2FoldChange, -log10(padj))) + 
  #volcanoplot with log2Foldchange versus pvalue
  geom_point(aes(col=sig)) + 
  #add points colored by significance
  scale_color_manual(values=c("gray", "chocolate")) + 
  ggtitle("Vocano plot: eeTreatment effect and statistical significance")

 #adding text for the top 20 genes
 top <- subset(input, padj<=1e-12 & abs(log2FoldChange) >= 2)
 volc <- volc + geom_text_repel(data=top, aes(label=gene)) 

 #ggsave("Volcanoplot.jpeg", device="jpeg") #In case you want to easily save to disk

 print(volc)

