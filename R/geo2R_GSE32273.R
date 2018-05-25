# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Thu Apr 12 05:31:16 EDT 2018

################################################################
#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)

setwd("/Users/jmx139/Desktop/IBD/geo_microarray_analysis")
#load external functions
#source(functions.R)

# load series and platform data from GEO
acc="GSE32273"
gset <- getGEO("GSE32273", GSEMatrix =TRUE, AnnotGPL=FALSE)
if (length(gset) > 1) idx <- grep("GPL8786", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- paste0("00000000000000000000XX11111111111111111111XX222222",
               "22222222222222XX33333333333333333333XX444444444444",
               "44444444XX55555555555555555555XX")

sml <- c()
labels <- c("UC_platelets","control_platelets", "UC_microvesicles","control_microvesicles", "UC_PBMC", "control_PBMC")
conditions <- c("UCpl","Cpl","UCm","Cm","UCP","CP")
for (i in 1:nchar(gsms)) { sml[i] <- conditions[as.integer(substr(gsms,i,i))+1] }
#for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

annot <- fData(gset)
annot <- annot[annot$SEQUENCE!="",]
annot <- subset(annot, select=c("ID","Species.Scientific.Name","SEQUENCE"))
write.table(annot, file=paste("results/",acc,"_probe_annot.txt",sep=""), row.names=F, col.names=F,sep="\t")


#probe sequences from all organisms have been mapped with bowtie to the human mature miRNA set 
#from miRBase v.22, allowing 1 mismatch
annot_mapping <- read.csv(file="data/probe_mapping/GSE32273_probe_annot_aligned_modif.arf", sep="\t", header=FALSE)
#select only mappings without mismatches
annot_mapping <- annot_mapping[annot_mapping$V13=="0",]
#when same sequence corresponding to different miRNAs map to the same miRBase sequence, keep only one human miRNA, where available
annot_mapping <- annot_mapping[order(annot_mapping$V6,annot_mapping$V1),]
annot_mapping <- annot_mapping[!duplicated(annot_mapping$V6),][,-1]
#remove trailing species name (added by bowtie?)
annot_mapping$V2 <-gsub("(.*)_.*","\\1",annot_mapping$V2)

#remove control probes and other probes not of interest (only keep human related probesets)
fData(gset)<-fData(gset)[fData(gset)$SEQUENCE!="",]
fData(gset)<-fData(gset)[fData(gset)$Species.Scientific.Name=="Homo sapiens",]
fData(gset)<-fData(gset)[fData(gset)$Sequence.Type=="miRNA",]
rows_to_keep<-rownames(fData(gset))
gset<-gset[rows_to_keep,]

# eliminate samples marked as "X"
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# set up the data and proceed with analysis
#sml <- paste("G", sml, sep="")    # set group names
#gname <-("")
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
#cont.matrix <- makeContrasts(G1-G0, G3-G2, G5-G4, levels=design)
cont.matrix <- makeContrasts(UCpl-Cpl, UCm-Cm, UCP-CP, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
cont_names <-c("UCpl_vs_Cpl", "UCm_vs_Cm", "UCP_vs_CP")

#thr = 0.1
for (i in 1:length(cont_names)) {
  tname=cont_names[i]
  tname_merged=paste(cont_names[i],"merged",sep="_")
  t=data.frame()
  t=topTable(fit2, coef=i, adjust="fdr", sort.by="B", number=5000)
  t$Species.Scientific.Name[t$Species.Scientific.Name == "Homo sapiens"] <- "AA Homo sapiens"
  t<-t[order(t$SEQUENCE,t$Species.Scientific.Name),]
  t<-t[!duplicated(t$SEQUENCE),]
  t<-t[t$adj.P.Val<thr,]
  t<-t[order(-t$B),]
  #t<- t[t$Sequence.Type=="miRNA",]
  t_merged <- as.data.frame(merge( t[c(1,4,6,9,10,13,12,16)], annot_mapping[c(1,5,6)], by.x='ID',by.y='V2', all.x=TRUE))
  t_merged <- t_merged[order(t_merged$adj.P.Val),]
  t_merged <- t_merged[t_merged$Sequence.Type=="miRNA",]
  t_merged$common_id <- tolower(gsub("_.*","",t_merged$V7))
  #change column names such as to match with deseq results 
  names(t)[names(t) == "adj.P.Val"] = "padj"
  names(t_merged)[names(t_merged) == "adj.P.Val"] = "padj"
  names(t)[names(t) == "AveExpr"] = "avg.expr"
  names(t_merged)[names(t_merged) == "AveExpr"] = "avg.expr"
  assign(tname, t)
  assign(tname_merged,t_merged)
  write.table(t[c(1,13,12,16,1)], file=paste("results/",acc,"_",tname,"_p_",thr, ".tsv",sep=""), row.names=F, sep="\t")
  write.table(t_merged[c(1,6,7,8,11)], file=paste("results/",acc,"_",tname,"_p_",thr, ".tsv",sep=""), row.names=F, sep="\t")
}

#common_GSE48957_GSE32273 <- as.data.frame(merge(T_UCa_vs_C, T_UCpl_vs_Cpl, by.x='ID', by.y='ID')[c(1,12,13,16,29,30,33)])

probe_boxplot <- function(myrow){
  mygene <- exprs(gset)[myrow,]
  #groups <- pData(gset)[,40]
  groups <- paste (pData(gset)[ ,40], pData(gset)[ ,41], sep = "" )
  df <- data.frame(values = mygene, vars = groups)
  print(groups)
  par(cex.axis=0.8) 
  boxplot(values ~ vars, data = df)
}
#probe_boxplot("characteristics_ch1.1","SNORD123_st")
#probe_boxplot("hsa-miR-501-5p_st")
#probe_boxplot("hsa-miR-27a-star_st")
#probe_boxplot("hsa-miR-4284_st")
#exprs(gset)["hsa-miR-4284_st",]

#summ <- summary(decideTests(fit2, adjust.method="BH", p.value=0.05))
#summ

