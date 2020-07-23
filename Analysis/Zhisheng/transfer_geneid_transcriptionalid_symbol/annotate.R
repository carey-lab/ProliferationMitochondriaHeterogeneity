
#使用org.Mm.eg.db，没有用biomart


data <- read.delim("RNASeqESC_Replicate1.tab.tsv",header=T,row.names=1,check.names=FALSE)


library(org.Mm.eg.db)

symb <- mapIds(org.Mm.eg.db, keys=rownames(sce), keytype="ENSEMBL", column="SYMBOL")
rowData(sce)$ENSEMBL <- rownames(sce)
rowData(sce)$SYMBOL <- symb
head(rowData(sce))