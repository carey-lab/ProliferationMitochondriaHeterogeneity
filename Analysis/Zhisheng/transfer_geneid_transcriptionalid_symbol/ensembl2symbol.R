require("biomaRt")
listMarts()
datasets <- listDatasets()
mart <- useMart("ENSEMBL_MART_ENSEMBL")
listmart<-listDatasets(mart)
View(listmart)
mart <- useDataset("mmusculus_gene_ensembl", mart)

esc1ens <- rownames(esc1)

esc1Lookup <- gsub("\\.[0-9]*$", "", esc1ens)


annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_transcript_id","ensembl_gene_id", "gene_biotype", "external_gene_name"),
  filter="ensembl_transcript_id",
  values=esc1Lookup,
  uniqueRows=T)

annotLookup <- data.frame(
  esc1ens[match(annotLookup$ensembl_transcript_id, esc1Lookup)],
  annotLookup)

colnames(annotLookup) <- c(
  "original_id",
  c("ensembl_transcript_id","ensembl_gene_id", "gene_biotype", "external_gene_name"))

write.table(annotLookup,file = "./ESC1NOTUNIQ.table")

head(annotLookup)
length(rownames(annotLookup))



#----本来是按照表达量排序然后做GO分析，但是这个只是排序---------
#---------------------------ESC1---------------------------------
esc1 <- read.delim("RNASeqESC_Replicate1.tab.tsv",header=T,row.names=1,check.names=FALSE,sep=" ")

esc1log2 <- log(esc1+0.1,2)
write.table(esc1log2,file = "esc1log2.table")


esc1log2slow <- as.data.frame(esc1log2$ESC_slow_TPM)
esc1log2slow$tnames <- row.names(esc1log2)
colnames(esc1log2slow) <- c("esc1log2slowTPM","tname")

head(esc1log2slow)
esc1log2slow <- esc1log2slow[order(esc1log2slow$esc1log2slowTPM,decreasing = T),]

name <- annotLookup2$external_gene_name[match(annotLookup2$original_id, esc1log2slow$tname)]

head(name)

write.table(name,file = "esc1log2slowTPM.list",row.names = F,quote = F)



colnames()


#--------------------------FIB1-----------------------------------

fib1 <- read.delim("RNASeqFIB_Replicate1.tab.tsv",header=T,row.names=1,check.names=FALSE,sep=" ")

#将表达量加上0.1，然后log2
fib1log2 <- log(fib1+0.1,2)
#保存到table
write.table(fib1log2,file = "fib1log2.table")

#取出slow column-----------如何取出带行列名的subtable?
fib1log2slow <- as.data.frame(fib1log2$FIB_slow_TPM)
#将行名赋值给取出列
fib1log2slow$tnames <- row.names(fib1log2)
#列名
colnames(fib1log2slow) <- c("fib1log2slowTPM","tname")



#将表达量按照从高到低排序
fib1log2slow <- fib1log2slow[order(fib1log2slow$fib1log2slowTPM,decreasing = T),]
head(fib1log2slow)

#从annotate table 中取出 gene symbol按照 表达量排序
name <- annotLookup2$external_gene_name[match(annotLookup2$original_id, fib1log2slow$tname)]

head(name)

write.table(name,file = "fib1log2slowTPM.list",row.names = F,quote = F)
















#--------------------------------------------------------
annotLookup2 <- read.csv("ESC1.table",sep=" ")

symbol <- annotLookup2$external_gene_name[match(annotLookup2$ensembl_transcript_id, rownames(esc1))]


rownames(esc1) <- symbol

head(esc1)


length(rownames(esc1))
length(symbol)
