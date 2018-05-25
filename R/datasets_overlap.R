##################################################
# Project: Analysis of public miRNA data
# Author(s): Oana Palasca
# Date: 25.05.2018
# Script purpose: Explore the overlap between comparisons coming from same/different studies 
##################################################

#install.packages('VennDiagram')
library(VennDiagram)
require(reshape2)
setwd("~/Desktop/IBD/geo_microarray_analysis")

# Run the analysis files for each dataset, for creating the results tables of interest 
thr = 1 # threshold on the adj. p value; works as a parameter for the scripts below
source("R/deseq_GSE89667.R")
source("R/geo2R_GSE48957.R")
source("R/geo2R_GSE32273.R")

add_comp <- function(acc, cont_names, meta_comp){
  #create a table containing study accession numbers and comparison names
  tmp <- melt(data.frame(rep(acc,3), cont_names))
  meta_comp <- rbind(meta_comp,tmp)
  return(meta_comp)
}

#create a meta-data table containing study accession numbers and comparison names
mc <- data.frame()
acc1 = "GSE89667"
cont_names1 <- c("UC_vs_DD","CD_vs_DD", "UC_vs_CD")
acc2 = "GSE48957"
cont_names2 <- c("UCa_vs_C", "UCa_vs_UCi", "UCi_vs_C")
acc3 = "GSE32273"
cont_names3 <- c("UCpl_vs_Cpl", "UCm_vs_Cm", "UCP_vs_CP")
mc <- add_comp(acc1, cont_names1, mc)
mc <- add_comp(acc2, cont_names2, mc)
mc <- add_comp(acc3, cont_names3, mc)
colnames(mc) <- c("study","comparison") 

for(i in 1:nrow(mc)) {
    acc <- mc[i,"study"]
    cp <- mc[i,"comparison"]
    rname=paste(acc, "_", cp, sep="")
    res <- read.table(file=paste("results/",acc,"_", cp, "_p_", thr, ".tsv", sep=""), sep="\t", row.names=1, header=TRUE)
    assign(rname, res)
}

for(i in 1:nrow(mc)) {
  for (j in (i+1):nrow(mc))
    acc1 <- mc[i,"study"]
    cp1 <- mc[i,"comparison"]
    acc2 <- mc[j,"study"]
    cp2 <- mc[j,"comparison"] 
    res1 <- read.table(file=paste("results/",acc1,"_", cp1, "_p_", thr, ".tsv", sep=""), sep="\t", row.names=1, header=TRUE)
    res2 <- read.table(file=paste("results/",acc2,"_", cp2, "_p_", thr, ".tsv", sep=""), sep="\t", row.names=1, header=TRUE)
    
}

#TO DO: create a function for intersecting any two result tables 
#and returning the correlation plot; run it in a loop for all result tables; 
#the "good" comparisons should give more meaningful results than the others  

GSE89667_UC_DD <- read.table(file=paste("results/",acc1,"_","UC_vs_DD","_p_",thr, ".tsv", sep=""), sep="\t", row.names=1, header=TRUE)
GSE48957_UCa_vs_C <- read.table(file=paste("results/",acc2,"_","UCa_vs_C","_p_",thr, ".tsv",sep=""),sep="\t", row.names=1, header=TRUE)  
GSE32273_UCpl_vs_Cpl <- read.table(file=paste("results/",acc3,"_","T_UCpl_vs_Cpl","_p_",thr, ".tsv",sep=""),sep="\t", row.names=1, header=TRUE)  

common_GSE48957_GSE89667 <- as.data.frame(merge(T_UCa_vs_C_merged, resOrdered_UC_DD, by.x='common_id', by.y='common_id', all.x=TRUE, all.y=TRUE)[c(1,2,7,8,12,13,9,17)])
#common_GSE48957_GSE89667 <- as.data.frame(merge(T_UCa_vs_C_merged, resOrdered_UC_DD, by.x='common_id', by.y='common_id')[c(1,2,7,8,12,13,9,17)])
#rownames(common_GSE48957_GSE89667)=common_GSE48957_GSE89667$common_id

common_GSE48957_GSE32273 <- as.data.frame(merge(GSE48957_UCa_vs_C, GSE32273_UCpl_vs_Cpl, by.x='ID', by.y='ID')[c(1,12,13,16,29,30,33,18)])
rownames(common_GSE48957_GSE32273) <- common_GSE48957_GSE32273$common_id

pdf("figures/logFC_correlation_GSE48957_GSE89667.pdf",width=6,height=5) 
plot(common_GSE48957_GSE89667$log2FoldChange, common_GSE48957_GSE89667$logFC, col="darkgrey", panel.first=grid(),
     main="", xlab="log2FC GSE89667", ylab="log2FC GSE48957 ",
     pch=20, cex=0.6)
abline(v=0)
abline(h=0)
with(subset(common_GSE48957_GSE89667, padj<0.1), points(log2FoldChange, logFC, pch=20, col="red", cex=0.7))
with(subset(common_GSE48957_GSE89667, adj.P.Val<0.1), points(log2FoldChange, logFC, cex=1.1))
#gn.selected <- (res$log2FoldChange > 1 & res$logFC < -1.2 ) 
#|| (res$log2FoldChange < -2 & res$logFC > 1.5 )
#text(common_GSE48957_GSE89667$log2FoldChange[gn.selected],
 #    common_GSE48957_GSE89667$logFC[gn.selected],
  #   lab=rownames(res)[gn.selected ], cex=0.4)
dev.off()

pdf("figures/logFC_correlation_GSE48957_GSE32273.pdf",width=6,height=6) 
res<-common_GSE48957_GSE32273
plot(res$logFC.y, res$logFC.x, col="darkgrey", panel.first=grid(),
     main="", xlab="log2FC GSE32273", ylab="log2FC GSE48957 ",
     pch=20, cex=0.6)
abline(v=0)
abline(h=0)
with(subset(res, adj.P.Val.y<0.1), points(logFC.y, logFC.x, pch=20, col="red", cex=0.7))
with(subset(res, adj.P.Val.x<0.1), points(logFC.y, logFC.x, cex=1.1))
#gn.selected <- (res$log2FoldChange > 1 & res$logFC < -1.2 ) 
#|| (res$log2FoldChange < -2 & res$logFC > 1.5 )
#text(common_GSE48957_GSE89667$log2FoldChange[gn.selected],
#    common_GSE48957_GSE89667$logFC[gn.selected],
#   lab=rownames(res)[gn.selected ], cex=0.4)
dev.off()



plot(common_GSE48957_GSE89667$logFC,common_GSE48957_GSE89667$log2FoldChange)
plot(common_GSE48957_GSE89667$adj.P.Val ,common_GSE48957_GSE89667$padj)

common_GSE48957_GSE89667_filtered <- as.data.frame(merge(T_UCa_vs_C_merged, resOrdered_UC_DD, by.x='common_id', by.y='common_id')[c(1,2,7,8,12,13,9,17)])
write.table(common_GSE48957_GSE89667_filtered, file="common_GSE48957_GSE89667_filtered.tsv", sep="\t")

common_GSE48957_GSE89667_all <- as.data.frame(merge(T_UCa_vs_C, resOrdered_UC_DD, by.x='common_id', by.y='common_id')[c(1,2,7,8,12,13,9,17)])

common_GSE48957_GSE89667 <- as.data.frame(merge(T_UCa_vs_C_merged, resOrdered_UC_DD, by.x='common_id', by.y='common_id', all.x=TRUE, all.y=TRUE)[c(1,2,7,8,12,13,9,17)])
common_GSE48957_GSE89667_filtered <- as.data.frame(merge(T_UCa_vs_C_merged, resOrdered_UC_DD, by.x='common_id', by.y='common_id')[c(1,2,7,8,12,13,9,17)])
write.table(common_GSE48957_GSE89667_filtered, file="common_GSE48957_GSE89667_filtered.tsv", sep="\t")

common_GSE48957_GSE89667_all <- as.data.frame(merge(T_UCa_vs_C, resOrdered_UC_DD, by.x='common_id', by.y='common_id')[c(1,2,7,8,12,13,9,17)])

  



a <- dim(resOrdered_UC_DD)[1]
b <- dim(T_UCa_vs_C_merged)[1]
ab <- dim(common_GSE48957_GSE89667_filtered)[1]

a<- 106
b<- 100
ab<- 21 
pdf(overlap_GSE48957_GSE89667.pdf,width=6,height=4) 

grid.newpage()

plot1 <- draw.pairwise.venn(a, b, ab, category = c("GSE89667 (RNA-Seq, colon)", "GSE48957 (array, colon)"), 
                   lty = rep("blank", 2), 
                   cex=2,
                   fill = c("light blue", "skyblue"), alpha = rep(0.5, 2), 
                   cat.pos = c(0, 0), cat.dist = rep(0.025, 2))


pdf("overlap_GSE48957_GSE89667.pdf",width=6,height=4) 
grid.draw(plot1);
dev.off()


a <- dim(T_UCpl_vs_Cpl)[1]
b <- dim(T_UCa_vs_C)[1]
ab <- dim(common_GSE48957_GSE32273)[1]
plot2 <- draw.pairwise.venn(a, b, ab, category = c("GSE32273 (array, blood)", "GSE48957 (array, colon)"), 
                            lty = rep("blank", 2), cex=2,
                            fill = c("sea green", "skyblue"), alpha = rep(0.5, 2), 
                            cat.pos = c(0, 0), cat.dist = rep(0.025, 2))

pdf("overlap_GSE32273_GSE48957.pdf",width=6,height=4) 
grid.draw(plot2);
dev.off()

#7 of the common ones are downregulated in colon, but upregulated in blood
common_GSE48957_GSE32273 <- as.data.frame(merge(T_UCa_vs_C, T_UCpl_vs_Cpl, by.x='ID', by.y='ID')[c(1,12,13,16,29,30,33)])
