#install.packages("Vennerable", repos="http://R-Forge.R-project.org") 
#source("https://bioconductor.org/biocLite.R")
#biocLite("affy")
#biocLite("DESeq2")
#biocLite("pheatmap")

setwd("~/Desktop/IBD/DESeq2")
require("DESeq2")

cts <- read.csv(file="miRNAs_expressed_all_samples_26_04_2018_t_13_43_32.csv", sep="\t")
#remove duplicate mature miRNAs 
cts<-cts[order(cts$X.miRNA,-cts$read_count),]
cts<-cts[!duplicated(cts$X.miRNA),]
#cts <- as.matrix(cts)
cts<-cts[,c(1,5:41)]
# Assign the first column as row names
cts2 <- cts[,-1]
cts2 <- as.matrix(cts2)
rownames(cts2) <- cts[,1]

coldata <- read.csv(file="config2.txt", sep="\t", header=FALSE)
coldata <- coldata[,c(2,3,4)]
coldata2 <- coldata[,-1]
rownames(coldata2) <- coldata[,1]
colnames(coldata2) <- c("condition","type")

#colnames(coldata2) <- c("condition")

#arrange columns of cts in the same order as rows of colnames 
all(rownames(coldata2) %in% colnames(cts2))
cts2 <- cts2[, rownames(coldata2)]
all(rownames(coldata2) == colnames(cts2))

dds <- DESeqDataSetFromMatrix(countData = cts2,
                              colData = coldata2,
                              design = ~ condition)
dds <- DESeq(dds)

res_UC_DD <- results(dds, contrast=c("condition","UC","DD"))
res_CD_DD <- results(dds, contrast=c("condition","CD","DD"))
res_UC_CD <- results(dds, contrast=c("condition","UC","CD"))

resOrdered_UC_DD <- as.data.frame(res_UC_DD[order(res_UC_DD$pvalue),])
resOrdered_CD_DD <- as.data.frame(res_CD_DD[order(res_CD_DD$pvalue),])
resOrdered_UC_CD <- as.data.frame(res_UC_CD[order(res_UC_CD$pvalue),])

write.table(resOrdered_UC_DD, file="condition_UC_vs_DD.csv", sep="\t")
write.table(resOrdered_CD_DD, file="condition_CD_vs_DD.csv", sep="\t")
write.table(resOrdered_UC_CD, file="condition_UC_vs_CD.csv", sep="\t")

resOrdered_UC_DD$mirna <-rownames(resOrdered_UC_DD)
resOrdered_CD_DD$mirna <-rownames(resOrdered_CD_DD)
resOrdered_UC_CD$mirna <-rownames(resOrdered_UC_CD)

thr<-0.05

resOrdered_UC_DD <- subset(resOrdered_UC_DD, padj < thr)
resOrdered_CD_DD <- subset(resOrdered_CD_DD, padj < thr)
resOrdered_UC_CD <- subset(resOrdered_UC_CD, padj < thr)

resOrdered_UC_DD$common_id <-gsub("_.*","",resOrdered_UC_DD$mirna)
resOrdered_CD_DD$common_id <-gsub("_.*","",resOrdered_CD_DD$mirna)
resOrdered_UC_CD$common_id <-gsub("_.*","",resOrdered_UC_CD$mirna)

common_1_2 <- as.data.frame(merge(resOrdered_UC_DD[c(1,2,6,7)], resOrdered_CD_DD[c(1,2,6,7)], by=c('mirna','mirna'))[c(1,2,3,4,6,7)])
common_1_3 <- as.data.frame(merge(resOrdered_UC_DD[c(1,2,6,7)], resOrdered_UC_CD[c(1,2,6,7)], by=c('mirna','mirna'))[c(1,2,3,4,6,7)])
common_2_3 <- as.data.frame(merge(resOrdered_CD_DD[c(1,2,6,7)], resOrdered_UC_CD[c(1,2,6,7)],  by=c('mirna','mirna'))[c(1,2,3,4,6,7)])
common_1_2_3 <- as.data.frame(merge(common_1_2, resOrdered_UC_CD[c(1,2,6,7)], by=c('mirna','mirna')))



vennD=Venn(SetNames = c("Samp1", "Samp2","Samp3"), Weight=c(`100`=x,`010`=x,`110`=x,`001`=x,`101`=x,`011`=x,`111`=x))


summary(resOrdered)
sum(res$padj < 0.1, na.rm=TRUE)
plotMA(res)
resNorm <- lfcShrink(dds, coef=2, type="normal")
plotMA(resNorm)
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")

library("pheatmap")
ntd <- normTransform(dds)
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("condition","type")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=FALSE, annotation_col=df)


vsd <- vst(dds, blind=FALSE)

library("RColorBrewer")
vsd<-ntd
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

vsd<-ntd
library("affy")
plotPCA(vsd, intgroup=c("condition", "type"))
