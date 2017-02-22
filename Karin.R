#source("https://bioconductor.org/biocLite.R")
#biocLite("limma")
#biocLite("edgeR")
library(limma)
library(edgeR)
#biocLite("DESeq2")


fc=read.delim("~/Desktop/Karin/chr_megaplasmid_counts_fC.txt", row.names= 1)
dge=NA
fit=NA
design<-NA
#Create a edgeR list. 6 is counts.DD.1, 14 is counts.wt.3
dge<- DGEList(counts=fc[, 6:20], group=groups) 
#Filtering
dim(dge)
isexpr <- rowSums(cpm(dge) > 10) >= 3 #3 replicates/condition
dge <- dge[isexpr,]
dim(dge)
head(dge$counts)
#Normalization by trimmed mean of M values (TMM)
dge <- calcNormFactors(dge)
write.table(dge$counts,file="/Users/admin/Desktop/Karin/chr_megaplasmid_normalized.xls",sep="\t",quote=F)
head(dge$counts)
################################################
#Normalization
logCPM <- cpm(dge, log=TRUE, prior.count=3)

groups=c("DD", "DD", "DD", "fur", "fur", "fur", "wt", "wt", "wt","pKM033", "pKM033", "pKM033", "pRyhB", "pRyhB", "pRyhB")
#design <- model.matrix(~ 0+factor(c("DD", "DD", "DD", "fur", "fur", "fur", "wt", "wt", "wt", "pKM033", "pKM033", "pKM033", "pRyhB", "pRyhB", "pRyhB")))

design <- model.matrix(~0+factor(c(rep("DD", 3), rep("fur", 3), rep("wt", 3), rep("pKM033", 3), rep("pRyhB", 3))))
colnames(design) <- c("groupDD", "groupfur", "grouppKM033", "grouppRyhB", "groupwt")
design

fit <- lmFit(logCPM, design)
fit_b <- eBayes(fit, trend=TRUE)
logCPM["ryhB",]
write.table(logCPM,file="/Users/admin/Desktop/Karin/chr_megaplasmid_normalized.xls",sep="\t",quote=F)

#####intentando sacar cada par
contr.matrix<- makeContrasts(
  furvswt = groupfur-groupwt, #1
  DDvsfur = groupDD-groupfur, #2
  DDvswt = groupDD-groupwt, #3
  DDvspRyhB = groupDD-grouppRyhB, #4
  DDvspKM033 = groupDD-grouppKM033, #5
  furvspRyhB = groupfur-grouppRyhB, #6
  furvspKM033 = groupfur-grouppKM033, #7
  pRyhBvswt = grouppRyhB-groupwt, #8
  pKM033vswt = grouppKM033-groupwt, #9
  pRyhBvspKM033 = grouppRyhB-grouppKM033, #10
  levels = design
)

fit2 <-contrasts.fit(fit_b, contr.matrix )
fit3 <- eBayes(fit2, trend=TRUE)
# p.value= 0.05, lfc=1,
s1<-topTable(fit3, p.value= 0.05, lfc=1, coef=1, number=4000)
write.table(s1, file="/Users/admin/Desktop/Karin/DE_limma/DE_furvswt.xls",sep="\t",quote=F)

s2<-topTable(fit3,  p.value= 0.05, lfc=1, coef=2, number=4000)
write.table(s2, file="/Users/admin/Desktop/Karin/DE_limma/DE_DDvsfur.xls",sep="\t",quote=F)

s3<-topTable(fit3,  p.value= 0.05, lfc=1, coef=3, number=4000)
write.table(s3, file="/Users/admin/Desktop/Karin/DE_limma/DE_DDvswt.xls",sep="\t",quote=F)

s4<-topTable(fit3, p.value= 0.05, lfc=1, coef=4, number=5000)
write.table(s4, file="/Users/admin/Desktop/Karin/DE_limma/DE_DDvspRyhBT.xls",sep="\t",quote=F)

s5<-topTable(fit3,  p.value= 0.05, lfc=1, coef=5, number=5000)
write.table(s5, file="/Users/admin/Desktop/Karin/DE_limma/DE_DDvspKM033.xls",sep="\t",quote=F)

s6<-topTable(fit3, p.value= 0.05, lfc=1, coef=6, number=5000)
write.table(s6, file="/Users/admin/Desktop/Karin/DE_limma/DE_furvspRyhB.xls",sep="\t",quote=F)

s7<-topTable(fit3,  p.value= 0.05, lfc=1, coef=7, number=5000)
write.table(s7, file="/Users/admin/Desktop/Karin/DE_limma/DE_furvspKM033.xls",sep="\t",quote=F)

s8<-topTable(fit3,  p.value= 0.05, lfc=1, coef=8, number=5000)
write.table(s8, file="/Users/admin/Desktop/Karin/DE_limma/DE_pRyhBvsWT.xls",sep="\t",quote=F)

s9<-topTable(fit3, p.value= 0.05, lfc=1, coef=9, number=5000)
write.table(s9, file="/Users/admin/Desktop/Karin/DE_limma/DE_pKM033vsWT.xls",sep="\t",quote=F)

s10<-topTable(fit3, p.value= 0.05, lfc=1, coef=10, number=5000)
write.table(s10, file="/Users/admin/Desktop/Karin/DE_limma/DE_pRyhBvspKM033.xls",sep="\t",quote=F)

s10["ryhB",]

s11<-topTable(fit3, p.value= 0.05, lfc=1,  number=5000)
write.table(s11, file="/Users/admin/Desktop/Karin/DE_all_contrasts.xls",sep="\t",quote=F)

s11["ryhB",]

plotMD(dge, column=3)
abline(h=0, col="red", lty=2, lwd=2)


#Sample clustering Plot MD, repeat with all the columns
for (i in 1:15){
  filename=paste("~/Desktop/Karin/", i,"_plotMD.jpeg", sep="" )
  jpeg(filename)
  plotMD(dge, column=i)
  abline(h=0, col="red", lty=2, lwd=2)
  dev.off()}

#GRAPH cluster groups
pch <- c(0,1,2,7,8 , 9, 10, 11, 12,15,16,17, 20, 22, 24)
colors <- rep(c("darkgreen", "red", "blue"), 2)
plotMDS(dge)#, col=colors[group], pch=pch[group])
#legend("topleft", legend=levels(design), pch=pch, col=colors, ncol=2)

############################################## ORIGINAL

# p.value= 0.05, lfc=1,
s1<-topTable(fit3, coef=1, number=4000)
write.table(s1, file="/Users/admin/Desktop/Karin/DE_limma/DE_furvswt.xls",sep="\t",quote=F)

s2<-topTable(fit3,   coef=2, number=4000)
write.table(s2, file="/Users/admin/Desktop/Karin/DE_limma/DE_DDvsfur.xls",sep="\t",quote=F)

s3<-topTable(fit3,  coef=3, number=4000)
write.table(s3, file="/Users/admin/Desktop/Karin/DE_limma/DE_DDvswt.xls",sep="\t",quote=F)

s4<-topTable(fit3, coef=4, number=5000)
write.table(s4, file="/Users/admin/Desktop/Karin/DE_limma/DE_DDvspRyhBT.xls",sep="\t",quote=F)

s5<-topTable(fit3, coef=5, number=5000)
write.table(s5, file="/Users/admin/Desktop/Karin/DE_limma/DE_DDvspKM033.xls",sep="\t",quote=F)

s6<-topTable(fit3, coef=6, number=5000)
write.table(s6, file="/Users/admin/Desktop/Karin/DE_limma/DE_furvspRyhB.xls",sep="\t",quote=F)

s7<-topTable(fit3, coef=7, number=5000)
write.table(s7, file="/Users/admin/Desktop/Karin/DE_limma/DE_furvspKM033.xls",sep="\t",quote=F)

s8<-topTable(fit3,  coef=8, number=5000)
write.table(s8, file="/Users/admin/Desktop/Karin/DE_limma/DE_pRyhBvsWT.xls",sep="\t",quote=F)

s9<-topTable(fit3, coef=9, number=5000)
write.table(s9, file="/Users/admin/Desktop/Karin/DE_limma/DE_pKM033vsWT.xls",sep="\t",quote=F)

s10<-topTable(fit3, coef=10, number=5000)
write.table(s10, file="/Users/admin/Desktop/Karin/DE_limma/DE_pRyhBvspKM033.xls",sep="\t",quote=F)

