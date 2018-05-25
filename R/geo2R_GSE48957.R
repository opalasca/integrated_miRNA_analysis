# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Sun Apr 15 08:55:48 EDT 2018


# http://bioinformatics-core-shared-training.github.io/cruk-bioinf-sschool/Day1/bioc-intro.pdf

################################################################
#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)
library("dplyr")


setwd("~/Desktop/IBD/geo_microarray_analysis")

# load series and platform data from GEO
acc="GSE48957"
gset <- getGEO("GSE48957", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL14613", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
#gsms <- "222222222200000000001111111"
gsms <- "000000000011111111112222222"
#           G0          G1      G2
sml <- c()
labels <- c("Control","UC_active","UC_inactive")
conditions <- c("C","UCa","UCi")
for (i in 1:nchar(gsms)) { sml[i] <- conditions[as.integer(substr(gsms,i,i))+1] }

annot <- fData(gset)
annot <- annot[annot$Sequence!="",]
annot <- subset(annot, select=c("ID","Species.Scientific.Name","Sequence"))
write.table(annot, file=paste("data/",acc,"_probe_annot.txt",sep=""), row.names=F, col.names=F,sep="\t")

#probe sequences from all organisms have been mapped with bowtie to the human mature miRNA set 
#from miRBase v.22, allowing 1 mismatch
annot_mapping <- read.csv(file="data/probe_mapping/GSE48957_probe_annot_aligned_modif.arf", sep="\t", header=FALSE)
#select only mappings without mismatches
annot_mapping <- annot_mapping[annot_mapping$V13=="0",]
#when same sequence corresponding to different miRNAs map to the same miRBase sequence, keep only one human miRNA, where available
annot_mapping <- annot_mapping[order(annot_mapping$V6,annot_mapping$V1),]
annot_mapping <- annot_mapping[!duplicated(annot_mapping$V6),][,-1]
#remove trailing species name (added by bowtie?)
annot_mapping$V2 <-gsub("(.*)_.*","\\1",annot_mapping$V2)

#remove control probes and other probes not of interest (only keep human related probesets)
fData(gset)<-fData(gset)[fData(gset)$Sequence!="",]
fData(gset)<-fData(gset)[fData(gset)$Species.Scientific.Name=="Homo sapiens",]
fData(gset)<-fData(gset)[fData(gset)$Sequence.Type=="miRNA",]
rows_to_keep<-rownames(fData(gset))
gset<-gset[rows_to_keep,]

fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
#cont.matrix <- makeContrasts(G2-G0, G1-G0, G1-G2, levels=design)
cont.matrix <- makeContrasts(UCa-C, UCa-UCi, UCi-C, levels=design)
cont_names<-c("UCa_vs_C", "UCa_vs_UCi", "UCi_vs_C")
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)

#thr = 0.1
for (i in 1:length(cont_names)) {
  tname=cont_names[i]
  tname_merged=paste(cont_names[i],"merged",sep="_")
  t=data.frame()
  t=topTable(fit2, coef=i, adjust="fdr", sort.by="B", number=5000)
  t$Species.Scientific.Name[t$Species.Scientific.Name == "Homo sapiens"] <- "AA Homo sapiens"
  t<-t[order(t$Sequence,t$Species.Scientific.Name),]
  t<-t[!duplicated(t$Sequence),]
  t<-t[t$adj.P.Val<thr,]
  t<-t[order(-t$B),]
  #t<- t[t$Sequence.Type=="miRNA",]
  t_merged <- as.data.frame(merge( t[c(1,4,6,9,10,13,12,16)], annot_mapping[c(1,5,6)], by.x='ID',by.y='V2', all.x=TRUE))
  t_merged <- t_merged[order(t_merged$adj.P.Val),]
  t_merged <- t_merged[t_merged$Sequence.Type=="miRNA",]
  t_merged$common_id <- tolower(gsub("_.*","",t_merged$V7))
  t$common_id <- tolower(gsub("(.*)_.*","\\1",t$ID))
  names(t)[names(t) == "adj.P.Val"] = "padj"
  names(t_merged)[names(t_merged) == "adj.P.Val"] = "padj"
  names(t)[names(t) == "AveExpr"] = "avg.expr"
  names(t_merged)[names(t_merged) == "AveExpr"] = "avg.expr"
  assign(tname,t)
  assign(tname_merged,t_merged)
  write.table(t[c(1,13,12,16,1)], file=paste("results/",acc,"_",tname,"_p_",thr, ".tsv",sep=""), row.names=F, sep="\t")
  write.table(t_merged[c(1,6,7,8,11)], file=paste("results/",acc,"_",tname,"_p_",thr, ".tsv",sep=""), row.names=F, sep="\t")
}
  
#summ <- summary(decideTests(fit2, adjust.method="BH", p.value=0.05))
#summ

probe_boxplot <- function(myrow){
  mygene <- exprs(gset)[myrow,]
  groups <- pData(gset)$characteristics_ch1.2
  df <- data.frame(values = mygene,
                   vars = groups)
  print(groups)
  par(cex.axis=0.8) 
  boxplot(values ~ vars, data = df)
}
#probe_boxplot("characteristics_ch1.1","SNORD123_st")
#probe_boxplot("SNORD123_st")
#probe_boxplot("hsa-miR-378c_st")
#exprs(gset)["hsa-miR-4284_st",]



