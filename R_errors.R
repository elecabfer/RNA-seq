####ERROR
> matriz_meta[i,]<-as.vector(metas)
Error in matriz_meta[i, ] <- as.vector(metas) : 
  number of items to replace is not a multiple of replacement length
###SOLUTION: matriz_meta[,i]
################################################
> heatmap.2(matriz_prim,
+           #cellnote = r_data[1:2000,],
+           main = "Spearman Primary",     # heat map title                              #### CHANGE ###
+           density.info="none",  # turns off density plot inside color legend
+           trace="none",         # turns off trace lines inside the heat map
+           margins =c(12,9),     # widens margins around plot
+           na.color="black",
+           #       col=my_palette,       # use on color palette defined earlier
+           breaks=col_breaks,    # enable color transition at specified limits
+           dendrogram="both")#,     # only draw a row dendrogram
Error in seq.default(min.raw, max.raw, by = min(diff(breaks)/100)) : 
  wrong sign in 'by' argument
#Solution
QUitar breaks porque no estaba definido
  
#################################################
Error in heatmap.2(matriz_meta, main = "Spearman Primary", density.info = "none",  : 
  `x' must be a numeric matrix
#Solution
as.numeric
