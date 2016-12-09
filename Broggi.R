if (!require("RColorBrewer")){
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}
setwd("~/Desktop/Broggi")
data <- read.delim("lymph_blood_miRNA_comparison.txt", row.names = 1, header = TRUE)

### PCA 1: Lymph vs Plasma
values<-as.matrix(data[4:69,])
lymph<-as.vector(values[, 1:21])
plasma=as.vector(values[,22:42])
values.types<-data[1,]
matriz<-cbind(lymph, plasma)
#row.names(matriz)<-data[,1[2:]]
#colnames(matriz)<-data[1,]
####colores
colores<-values
colores[lymph]<-"yellow"
colores[plasma]<-"red"

######
matriz.pca <- prcomp(as.numeric(matriz))
summary(matriz.pca)
##pairs
pairs(values.pca$x,col=colores)
#
plot(values.pca, type="l")
summary(values.pca)

View(values)
############## Positive and Negative
pos_index<- which(data[3,]=="positive")
neg_index<- which(data[3,]=="negative")

positive <- as.vector(values[,pos_index])
negative <- as.vector(values[,neg_index])
matriz_posneg<-as.matrix(cbind(positive, negative))
pca.posneg<-prcomp(matriz_posneg)
