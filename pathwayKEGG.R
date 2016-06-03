source("http://bioconductor.org/biocLite.R")
setwd("~/Desktop/RNA/Selene/results/limma_means/")
XXX=read.delim("limma_sting-cgas.txt")
options(stringsAsFactors=F)
library(pathview)
expr1 = XXX$logFC #vector of log fold-changes
names(expr1) = XXX[,1]#Gene_ID
library(biomaRt)
ensembl = useMart("ensembl",dataset="mmusculus_gene_ensembl")
#ensembl = useMart("ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl", host="www.ensembl.org")
table1 = getBM(
  attr=c("ensembl_gene_id","entrezgene"),
  filter=c("ensembl_gene_id"),
  values=list(names(expr1)),
  mart=ensembl)
# Here you need the ids of the KEGG pathways you're looking at, find here: http://www.genome.jp/kegg-bin/get_htext
#04110: cell cycle
#04115: p53
#04210: apoptosis
library(org.Mm.eg.db)## Mus musculus
for (pid in c("04110", "04115", "04140", "04350")) {     #"00540", "04140", "04064",
  pv.out = pathview(gene.data=expr1, gene.idtype="ENSEMBL", pathway.id=pid,
                    species="mmu", kegg.native=T, same.layer=F, out.suffix="sting-cgas",
                    limit=list(gene=4,cpd=4))
  write.table(merge(XXX,pv.out$plot.data.gene[,c("kegg.names","labels")],by.x=2,by.y=1),
              file=paste(pid,"_genes.txt",sep=""),
              sep="\t",col.names=F,row.names=F,quote=F)
}
