##################################################
# Project: Analysis of public miRNA data
# Author(s): Oana Palasca
# Date: 25.05.2018
# Script purpose: Explore the overlap between comparisons coming from same/different studies 
##################################################

library(VennDiagram)
require(reshape2)
setwd("~/Desktop/IBD/integrated_miRNA_analysis")

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

draw_logFC_corr <- function(acc1, acc2, cp1, cp2, common, p1, p2){
  pdf(file=paste("figures/logFC_correlation_",acc1,"_",cp1,"_",acc2,"_",cp2,".pdf",sep=""),width=6,height=5) 
  plot(common$logFC.x, common$logFC.y, col="darkgrey", panel.first=grid(),
       main=paste(acc1,"_",cp1,"_",acc2,"_",cp2,sep=''), xlab=paste("logFC ",acc1,"_",cp1,sep=''), ylab=paste("logFC ",acc2,"_",cp2,sep=''),
       pch=20, cex=0.6)
  abline(v=0)
  abline(h=0)
  with(subset(common, padj.x<p1), points(logFC.x, logFC.y, pch=20, col="red", cex=0.7))
  with(subset(common, padj.y<p2), points(logFC.x, logFC.y, cex=1.1))
  selected <- abs(common$logFC.x) > 1.2 & abs(common$logFC.y) > 1.2 #& common$logFC.x/common$logFC.y > 0
  selected
  if (selected){
  text(common$logFC.x[selected], common$logFC.y[selected],
     lab=common$common_id[selected], cex=0.4)
  }
  dev.off()
}

setdiffDF <- function(A, B){ 
  f <- function(A, B) 
    A[!duplicated(rbind(B, A))[nrow(B) + 1:nrow(A)], ] 
  df1 <- f(A, B) 
  df2 <- f(B, A) 
  rbind(df1, df2) 
} 

draw_volcano_highlight_uncommon <- function(acc1, acc2, cp1, cp2, common_full, res1, res2, p1, p2){
  pdf(file=paste("figures/volcano_highlighted_",acc1,"_",cp1,"_",acc2,"_",cp2,".pdf",sep=""),width=6,height=5) 
  pthr <- 0.05 # Threshold on the adjusted p-value
  cols <- densCols(res1$logFC, -log10(res1$padj))
  plot(res1$logFC, -log10(res1$padj), col=cols, panel.first=grid(),
       main=paste(acc1,"_",cp1,sep=''), xlab="log2(fold-change)", ylab="-log10(adjusted p-value)",
       pch=20, cex=0.8)
#  with(subset(common, padj.x<pthr), points(logFC.x, -log10(padj.y), cex=1.1))
#  with(subset(common, padj.x<0.1 & padj.y<0.1), points(logFC.x, -log10(padj.x), cex=1.1))
  with(subset(common_full, padj.x<0.1 & is.na(padj.y)), points(logFC.x, -log10(padj.x), cex=1.1))
  
  abline(v=0)
  abline(v=c(-1,1), col="brown")
  abline(h=-log10(pthr), col="brown")
  #selected <- abs(res$logFC) > lfcthr & res$padj < pthr 
  #if (selecteds){
   # text(res$logFC[selected],-log10(res$padj)[selected],
    #     lab=res$common_id[selected], cex=0.4)
  #}
  dev.off()
}


draw_volcano <- function(acc, cp, res, pthr, lfcthr){
  pthr <- 0.05 # Threshold on the adjusted p-value
  cols <- densCols(res$logFC, -log10(res$padj))
  plot(res$logFC, -log10(res$padj), col=cols, panel.first=grid(),
     main=paste(acc,"_",cp,sep=''), xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
     pch=20, cex=0.6)
  abline(v=0)
  abline(v=c(-1,1), col="brown")
  abline(h=-log10(pthr), col="brown")
  selected <- abs(res$logFC) > lfcthr & res$padj < pthr 
  if (selected){
  text(res$logFC[selected],-log10(res$padj)[selected],
       lab=res$common_id[selected], cex=0.4)
  }
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

pdf(file=paste("figures/volcano_plots.pdf")) 
for(i in 1:nrow(mc)) {
    acc <- mc[i,"study"]
    cp <- mc[i,"comparison"]
    rname=paste(acc, "_", cp, sep="")
    res <- read.table(file=paste("results/",acc,"_", cp, "_p_", thr, ".tsv", sep=""), sep="\t", row.names=1, header=TRUE)
    assign(rname, res)
    draw_volcano(acc, cp, res, 0.05, 2)
    #TO DO draw volcano plots; colour by avg.expr parameter
}
dev.off()

pdf(file=paste("figures/logFC_correlation_all.pdf")) 

for(i in 1:(nrow(mc)-1)){
  for (j in (i):nrow(mc)){
    common <- data.frame()
    head(common)
    acc1 <- mc[i,"study"]
    cp1 <- mc[i,"comparison"]
    acc2 <- mc[j,"study"]
    cp2 <- mc[j,"comparison"] 
    print(paste(i,"_",j,"_",acc1,"_",cp1,"_",acc2,"_",cp2,sep="")) 
    res1 <- read.table(file=paste("results/",acc1,"_", cp1, "_p_", thr, ".tsv", sep=""), sep="\t", row.names=1, header=TRUE)
    res2 <- read.table(file=paste("results/",acc2,"_", cp2, "_p_", thr, ".tsv", sep=""), sep="\t", row.names=1, header=TRUE)
    #remove rows with NAs
    res1<-res1[complete.cases(res1),]
    res2<-res2[complete.cases(res2),]
    #intersection of results
    common <- as.data.frame(merge(res1, res2, by='common_id'))
    common$common_id <- sub("hsa-", "", common$common_id)
    common$common_id_trunc <- sub("-5p", "", common$common_id)
    #use full intersection file for plotting the disjoint sets as well
    common_full <- as.data.frame(merge(res1, res2, by='common_id', all.x=TRUE,all.y=TRUE))
    #rownames(common)=common$common_id
    #save files with results
    #TO DO: save publication-like tables with top common miRNAs
    write.table(common, file=paste("results/common_",acc1,"_",cp1,"_",acc2,"_",cp2,".tsv",sep=""), sep="\t")
    write.table( subset(common,(padj.x<0.05 & padj.y<0.2) | (padj.x<0.2 & padj.y<0.5 )), file=paste("results/common_",acc1,"_",cp1,"_",acc2,"_",cp2,"_significant.tsv",sep=""), sep="\t")
    #write.table( subset(common,(padj.x<0.05 & padj.y<0.15) ), file=paste("results/common_",acc1,"_",cp1,"_",acc2,"_",cp2,"_significant.tsv",sep=""), sep="\t")
    
    #draw the correlation plot with names near the points;
    #TO DO: draw volcano-like plots side by side showing
    #the miRNAs missing from the intersection in each of the two datasets; 
    #use different colour for the missing points
    draw_logFC_corr(acc1, acc2, cp1, cp2, common, 0.1, 0.1)
    draw_volcano_highlight_uncommon(acc1, acc2, cp1, cp2, common_full, res1, res2, 0.1, 0.1)
      
  }
}  
dev.off()

#TO DO: create a function for intersecting any two result tables 
#and returning the correlation plot; run it in a loop for all result tables; 
#the "good" comparisons should give more meaningful results than the others  
common_GSE89667_GSE48957_sig <- read.table(file="results/common_GSE89667_UC_vs_DD_GSE48957_UCa_vs_C_significant.tsv", row.names=1, header=TRUE)
common_GSE89667_GSE48957 <- read.table(file="results/common_GSE89667_UC_vs_DD_GSE48957_UCa_vs_C.tsv", row.names=1, header=TRUE)
curated_UC_up <-read.table(file="data/curated_UC_up.txt", sep=" ", header=FALSE)
curated_UC_up<- curated_UC_up[c(5,4)]
curated_UC_up$common_id <- sub("-5p", "", curated_UC_up$V5 )

curated_UC_down <-read.table(file="data/curated_UC_down.txt", sep=" ", header=FALSE)
curated_UC_down<- curated_UC_down[c(5,4)]
common_all <- as.data.frame(merge(common_GSE89667_GSE48957, curated_UC_up, by.x='common_id_trunc',by.y='common_id'))
common_all_sig <- as.data.frame(merge(common_GSE89667_GSE48957_sig, curated_UC_up, by.x='common_id_trunc',by.y='common_id'))

#rownames(common_GSE48957_GSE89667)=common_GSE48957_GSE89667$common_id
acc1="GSE89667"
acc2="GSE48957"
cp1="UC_vs_DD"
cp2="UCa_vs_C"

res1 <- read.table(file=paste("results/",acc1,"_", cp1, "_p_", thr, ".tsv", sep=""), sep="\t", row.names=1, header=TRUE)
res2 <- read.table(file=paste("results/",acc2,"_", cp2, "_p_", thr, ".tsv", sep=""), sep="\t", row.names=1, header=TRUE)
#remove rows with NAs
res1$common_id <- sub("hsa-", "", res1$common_id)
res1$common_id_trunc <- sub("-5p", "", res1$common_id)
res2$common_id <- sub("hsa-", "", res2$common_id)
res2$common_id_trunc <- sub("-5p", "", res2$common_id)
res1sig <- subset(res1, padj<0.1)
res2sig <- subset(res2, padj<0.1)

common_1 <- as.data.frame(merge(res1sig, curated_UC_up, by.x='common_id_trunc',by.y='common_id'))
common_2 <- as.data.frame(merge(res2sig, curated_UC_up, by.x='common_id_trunc',by.y='common_id'))








