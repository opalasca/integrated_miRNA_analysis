#install required packages

source("https://bioconductor.org/biocLite.R")
biocLite("GEOquery")
source("https://bioconductor.org/biocLite.R")
biocLite("limma")

#Datasets
# GSE22307 mouse DSS, 23 samples - normal, day2, day4, day6
# GSE67106 human lncRNA expression microarray, 96 samples, normal, CD, UC, inflamed/non-inflamed regions 

library(GEOquery)
require(Biobase)

gse_acc="GSE22307"

#  GSE class 
gse <- getGEO(gse_acc,GSEMatrix=FALSE)
head(Meta(gse))
names(GSMList(gse))
GSMList(gse)[[1]]
names(GPLList(gse))

# Getting GSE Series Matrix files as an ExpressionSet
gse2 <- getGEO(gse_acc,GSEMatrix=TRUE)
show(gse2)
show(pData(phenoData(gse2[[1]]))[1:5,c(1,6,8)])

#Convert GSE to ExpressionSet
gsmplatforms <- lapply(GSMList(gse),function(x) {Meta(x)$platform_id})
head(gsmplatforms)

platf='GPL1261'
#filter samples corresponding to given platform
gsmlist = Filter(function(gsm) {Meta(gsm)$platform_id==platf},GSMList(gse))
length(gsmlist)

#check single GSM
Table(gsmlist[[1]])[1:5,]
Columns(gsmlist[[1]])[1:5,]

# Make a matrix from the values column
# get the probeset ordering
probesets <- Table(GPLList(gse)[[1]])$ID
# make the data matrix from the VALUE columns from each GSM
# being careful to match the order of the probesets in the platform
# with those in the GSMs
data.matrix <- do.call('cbind',lapply(gsmlist,function(x) 
                                     {tab <- Table(x)
                                      mymatch <- match(probesets,tab$ID_REF)
                                      return(tab$VALUE[mymatch])
                                    }))
data.matrix <- apply(data.matrix,2,function(x) {as.numeric(as.character(x))})
data.matrix <- log2(data.matrix)
data.matrix[1:5,]

# go through the necessary steps to make a compliant ExpressionSet
rownames(data.matrix) <- probesets
colnames(data.matrix) <- names(gsmlist)
pdata <- data.frame(samples=names(gsmlist))
rownames(pdata) <- names(gsmlist)
pheno <- as(pdata,"AnnotatedDataFrame")
eset <- new('ExpressionSet',exprs=data.matrix,phenoData=pheno)
eset

getGEOSuppFiles
library(limma)
targets <- readTargets("mouse_GSE22307/targets.txt")

lev <- c("D0","D2","D4","D6")
f <- factor(targets$Target, levels=lev)
design <- model.matrix(~0+f)
colnames(design) <- lev
fit <- lmFit(eset, design)

cont <- makeContrasts(
        "D2-D0",
        "D4-D2",
        "D6-D4",
   levels=design)
fit2 <- contrasts.fit(fit, cont)
fit2 <- eBayes(fit2)
topTableF(fit2, adjust="BH")

biocLite()

