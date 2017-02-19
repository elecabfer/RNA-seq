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
fit <- lmFit(logCPM, design,)
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
fc=read.delim("~/Desktop/Karin/fC/megaplasmid_counts_fC.txt", row.names= 1)
v =NA
dge=NA
fit=NA

#Create a edgeR list. 6 is counts.DD.1, 14 is counts.wt.3
design <- model.matrix(~ 0+factor(c("pKM033", "pKM033", "pKM033", "pRyhB", "pRyhB", "pRyhB")))
dge<- DGEList(counts=fc[, 6:11]) 
#Filtering
dim(dge)
isexpr <- rowSums(cpm(dge) > 10) >= 3 #3 replicates/condition
dge <- dge[isexpr,]
dim(dge)
#Normalization by trimmed mean of M values (TMM)
dge <- calcNormFactors(dge)
colnames(design) <- c("grouppKM033", "grouppRyhB")
contr.matrix<- makeContrasts(
  pKM033vspRyhB = grouppKM033-grouppRyhB,
  levels = design
)
logCPM <- cpm(dge, log=TRUE, prior.count=3)
fitm <- lmFit(logCPM, design,)
fitm <- eBayes(fitm, trend=TRUE)
m<-topTable(fitm, coef=ncol(design), number= 4000)
write.table(m, file="/Users/admin/Desktop/Karin/fC/DE_pKM033_vs_pRyhB.xls",sep="\t",quote=F)

#####################################################################################################################
#####################################################################################################################

#Sample clustering Plot MD, repeat with all the columns
plotMD(v, column=9)
abline(h=0, col="red", lty=2, lwd=2)

#GRAPH cluster groups
pch <- c(0,1,2,15,16,17, 20, 22, 24)
colors <- rep(c("darkgreen", "red", "blue"), 2)
plotMDS(dge, col=colors[group], pch=pch[group])
legend("topleft", legend=levels(group), pch=pch, col=colors, ncol=2)

#linear model fitting and DEA
fit <- eBayes(lmFit(v,groups))
topTable(fit,coef=2)

#Alternative: differential expression voom
v <- voom(dge, group=~DD+fur, plot=TRUE)
voomDat <- voom(dge, group=DD+fur)
fit <- eBayes(lmFit(voomDat, group=~DD+fur))
topTable(fit, coef=2, n=10)
