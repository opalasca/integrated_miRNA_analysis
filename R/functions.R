probe_boxplot <- function(gset, myrow){
  mygene <- exprs(gset)[myrow,]
  groups <- pData(gset)$characteristics_ch1.2
  df <- data.frame(values = mygene,
                   vars = groups)
  print(groups)
  par(cex.axis=0.8) 
  boxplot(values ~ vars, data = df)
}

write_tables <- function(cont_names, fit2, acc) {
  for (i in 1:length(cont_names)) {
    tname=paste("T",cont_names[i],sep="_")
    print(i)
    t=data.frame()
    t=topTable(fit2, coef=i, adjust="fdr", sort.by="B", number=5000)
    t$Species.Scientific.Name[t$Species.Scientific.Name == "Homo sapiens"] <- "AA Homo sapiens"
    t<-t[order(t$Sequence,t$Species.Scientific.Name),]
    t<-t[!duplicated(t$Sequence),]
    t<-t[t$adj.P.Val<0.1,]
    assign(tname, t)
    write.table(t, file=paste(acc,"_",tname,".txt",sep=""), row.names=F, sep="\t")
  }
}

# log2 transform - already transformed
#ex <- exprs(gset)
#qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
#LogC <- (qx[5] > 100) ||
#  (qx[6]-qx[1] > 50 && qx[2] > 0) ||
#  (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
#if (LogC) { ex[which(ex <= 0)] <- NaN
#exprs(gset) <- log2(ex) }
