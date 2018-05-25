
setwd("~/Desktop/IBD/geo_microarray_analysis")

thr=1
source("R/deseq_GSE89667.R")

######
#Draw volcano plot with labels for most significantly DE genes for deseq obj

pdf("figures/volcano_plot_GSE89667.pdf",width=6,height=4) 
res<-res_UC_DD
alpha <- 0.05 # Threshold on the adjusted p-value
cols <- densCols(res$log2FoldChange, -log10(res$pvalue))
plot(res$log2FoldChange, -log10(res$padj), col=cols, panel.first=grid(),
     main="", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
     pch=20, cex=0.6)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")
gn.selected <- abs(res$log2FoldChange) > 2.5 & res$padj < alpha 
text(res$log2FoldChange[gn.selected],
     -log10(res$padj)[gn.selected],
     lab=rownames(res)[gn.selected ], cex=0.4)
dev.off()


thr = 1
source("R/geo2R_GSE48957.R")

pdf("figures/volcano_plot_GSE48957.pdf",width=6,height=4) 
res<-T_UCa_vs_C
rownames(res) <- res$common_id
alpha <- 0.05 # Threshold on the adjusted p-value
cols <- densCols(res$logFC, -log10(res$adj.P.Val))
plot(res$logFC, -log10(res$adj.P.Val), col=cols, panel.first=grid(),
     main="", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
     pch=20, cex=0.6)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")
gn.selected <- abs(res$logFC) > 2 & res$adj.P.Val < alpha 
text(res$logFC[gn.selected],
     -log10(res$adj.P.Val)[gn.selected],
     lab=rownames(res)[gn.selected ], cex=0.4)
dev.off()


#####
#Another type of volcano plot
par(mar=c(5,5,5,5), cex=1.0, cex.main=1.4, cex.axis=1.4, cex.lab=1.4)
topT <- as.data.frame(res_UC_DD)
#Adjusted P values (FDR Q values)
with(topT, plot(log2FoldChange, -log10(padj), pch=20, main="Volcano plot", cex=1.0, xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~Q~value)))
with(subset(topT, padj<0.05 & abs(log2FoldChange)>2), points(log2FoldChange, -log10(padj), pch=20, col="red", cex=0.5))
#with(subset(topT, padj<0.05 & abs(log2FoldChange)>2), text(log2FoldChange, -log10(padj), labels=subset(rownames(topT), topT$padj<0.05 & abs(topT$log2FoldChange)>2), cex=0.8, pos=3))
#Add lines for absolute FC>2 and P-value cut-off at FDR Q<0.05
abline(v=0, col="black", lty=3, lwd=1.0)
abline(v=-2, col="black", lty=4, lwd=2.0)
abline(v=2, col="black", lty=4, lwd=2.0)
abline(h=-log10(max(topT$pvalue[topT$padj<0.05], na.rm=TRUE)), col="black", lty=4, lwd=2.0)




summary(res_UC_DD)
sum(res_UC_DD$padj < 0.1, na.rm=TRUE)
plotMA(res_UC_DD)
resNorm <- lfcShrink(dds, coef=2, type="normal")
plotMA(resNorm)
plotCounts(dds, gene=which.min(res_UC_DD$padj), intgroup="condition")

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
