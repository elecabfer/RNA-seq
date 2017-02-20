source("https://bioconductor.org/biocLite.R")
#biocLite("limma")
#biocLite("edgeR")
library(limma)
library(edgeR)
biocLite("DESeq2")

#read data
fc=read.delim("~/Desktop/Karin/fC/chromosome_counts_fC.txt", row.names= 1)
v =NA
dge=NA
fit=NA

#Create a edgeR list. 6 is counts.DD.1, 14 is counts.wt.3
dge<- DGEList(counts=fc[, 6:14], group=groups) 
#Filtering
dim(dge)
isexpr <- rowSums(cpm(dge) > 10) >= 3 #3 replicates/condition
dge <- dge[isexpr,]
dim(dge)
write.table(dge$counts,file="/Users/admin/Desktop/Karin/fC/filtered_chromosome.xls",sep="\t",quote=F)
#Normalization by trimmed mean of M values (TMM)
dge <- calcNormFactors(dge)

################################################

design <- model.matrix(~ 0+factor(c("DD", "DD", "DD", "fur", "fur", "fur", "wt", "wt", "wt")))
colnames(design) <- c("groupDD", "groupfur", "groupwt")
#Normalization
logCPM <- cpm(dge, log=TRUE, prior.count=3)
fit <- lmFit(logCPM, design)
fit <- eBayes(fit, trend=TRUE)

#####intentando sacar cada par
contr.matrix<- makeContrasts(
  furvswt = groupfur-groupwt,
  DDvsfur = groupDD-groupfur,
  DDvswt = groupDD-groupwt,
  levels = design
)
fit2 <-contrasts.fit(fit, contr.matrix )
fit2 <- eBayes(fit2, trend=TRUE)

s1<-topTable(fit2, coef=1, number=4000)
write.table(s1, file="/Users/admin/Desktop/Karin/fC/DE_furvswt.xls",sep="\t",quote=F)

s2<-topTable(fit2, coef=2, number=4000)
write.table(s2, file="/Users/admin/Desktop/Karin/fC/DE_DDvsfur.xls",sep="\t",quote=F)

s3<-topTable(fit2, coef=3, number=4000)
write.table(s3, file="/Users/admin/Desktop/Karin/fC/DE_DDvswt.xls",sep="\t",quote=F)


######################### FINI!!! 
#RPKM
#dge_rpkm <- rpkm(dge,  fc$Length)
#head(dge_rpkm)

############################################################################## MEGAPLASMID
#read data
fcm=read.delim("~/Desktop/Karin/fC/megaplasmid_counts_fC.txt", row.names= 1)
v =NA
dge=NA
fit=NA
m=NA
#Create a edgeR list. 6 is counts.DD.1, 14 is counts.wt.3
design <- model.matrix(~ 0+factor(c("pKM033", "pKM033", "pKM033", "pRyhB", "pRyhB", "pRyhB")))
colnames(design) <- c("grouppKM033", "grouppRyhB")
dgem<- DGEList(counts=fcm[, 6:11]) 
#Filtering
dim(dgem)
isexpr <- rowSums(cpm(dgem) > 10) >= 3 #3 replicates/condition
dgem <- dgem[isexpr,]
dim(dgem)
#Normalization by trimmed mean of M values (TMM)
dgem <- calcNormFactors(dgem)
contr.matrix<- makeContrasts(
  pKM033vspRyhB = grouppKM033-grouppRyhB,
  levels = design
)
logCPM <- cpm(dgem, log=TRUE, prior.count=3)
fitm <- lmFit(logCPM, design)
fit2m <-contrasts.fit(fitm, contr.matrix )
fit2m <- eBayes(fit2m, trend=TRUE)

m<-topTable(fit2m, coef=1, number= 4000)
write.table(m, file="/Users/admin/Desktop/Karin/fC/DE_pKM033_vs_pRyhB.xls",sep="\t",quote=F)
