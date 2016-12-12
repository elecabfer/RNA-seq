if (!require("RColorBrewer")){
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}
setwd("~/Desktop/Broggi")
data <- read.delim("lymph_blood_miRNA_comparison.txt", row.names = 1, header = TRUE)
values<-as.matrix(data[4:69,])
d <- colnames(values)
values <- t(apply(values,1,as.numeric))
head(values)
values.types<-data[1,]
head(matriz)
##############  PCA

matriz.pca <- prcomp(t(values))
head(matriz.pca$x)

############## PLOT Lymph/Plasma
colores=c(rep("blue", 21), rep("red", 21))
formas=c(rep(16, 21), rep(16, 21))  ## 16 para circulos, 15 para cuadrados
plot(matriz.pca$x, col=colores, pch=formas, main="PCA Lymph-Plasma")
legend("topleft", c("Lymph", "Plasma"), col=c("blue", "red"), box.lty =0, pch=c(16,16), cex=0.9,pt.cex=0.7)

############### PLOT Neg/Pos
colores=c(rep("red", 9), rep("green", 12), rep("red", 9), rep("green", 12))
formas=c(rep(16, 9), rep(16, 12), rep(16, 9), rep(16, 12))  ## 16 para circulos, 15 para cuadrados, 17 para triangulos
plot(matriz.pca$x, col=colores, main= "PCA Positive-Negative" , pch=formas)
legend("topleft", c("Negative", "Positive"), col=c("red","green"), box.lty =0, pch=c(16,16), cex=0.9,pt.cex=0.7)


############### PLOT Todo
colores=c(rep("blue", 9), rep("blue", 12), rep("red", 9), rep("red", 12))
formas=c(rep(1, 9), rep(16, 12), rep(1, 9), rep(16, 12))  ## 16 para circulos, 15 para cuadrados
plot(matriz.pca$x, col=colores, pch=formas)
legend("topleft", c("Lymph negative", "Lymph positive", "Plasma negative", "Plasma positive" ), 
       col=c("blue", "blue", "red", "red"),  
       box.lty =0, pch=c(1, 16, 1, 16), cex=0.9,pt.cex=0.7)

