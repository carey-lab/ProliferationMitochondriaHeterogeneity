# Dec 31 2018
#merge 2016 and 2018 tpm data to do GSEA
#USE GENE ID

#-----read 2016 data---------------
fib1 <- read.table("./raw counts data from abundance h5/FIB_1_h5.s",header = T)
#remove all characters behind decimal
rownames(fib1) <- gsub("\\..*", "", rownames(fib1))
head(fib1)
nrow(fib1)
#get correspondent human gene symbol
require("biomaRt")

mart <- useMart("ENSEMBL_MART_ENSEMBL")
# listmart<-listDatasets(mart)
# View(listmart)
mart <- useDataset("mmusculus_gene_ensembl", mart)

geneid <- getBM(
  mart=mart,
  attributes=c("ensembl_transcript_id","ensembl_gene_id"),
  filter="ensembl_transcript_id",
  values=rownames(fib1),
  uniqueRows=F)

head(geneid)

head(fib1)
length(unique(rownames(fib1)))
nrow(geneid)
nrow(fib1)
#transfer transcripts id to gene id and merge transcripts id

tx2gene <- data.frame(geneid)
nrow(tx2gene)
# merge transcripts id
#只提取有老鼠基因名基因名的对应的老鼠转录ID
out <- lapply(unique(tx2gene$ensembl_gene_id),
              function(x) {
                tmp       <- tx2gene[tx2gene$ensembl_gene_id == x,]
                tmp_count <- fib1[match(tmp$ensembl_transcript_id,
                                        rownames(fib1)),]
                tmp_out   <- colSums(tmp_count)
                return(tmp_out)
              })

gene_counts <- matrix(unlist(out), 
                      ncol  = ncol(fib1), 
                      byrow = T)
rownames(gene_counts) <- unique(tx2gene$ensembl_gene_id)
colnames(gene_counts) <- colnames(fib1)
head(gene_counts)
fib1 <- gene_counts
#加入附加信息
nrow(fib1)
ncol(fib1)-2

#-----read 2018 data---------------
fib2 <- read.table("./FIB_2_h5.s",header = T)
#remove all characters behind decimal
rownames(fib2) <- gsub("\\..*", "", rownames(fib2))
head(fib2)
nrow(fib2)
#get correspondent human gene symbol
require("biomaRt")

mart <- useMart("ENSEMBL_MART_ENSEMBL")
# listmart<-listDatasets(mart)
# View(listmart)
mart <- useDataset("mmusculus_gene_ensembl", mart)

geneid <- getBM(
  mart=mart,
  attributes=c("ensembl_transcript_id","ensembl_gene_id"),
  filter="ensembl_transcript_id",
  values=rownames(fib2),
  uniqueRows=F)

head(geneid)

head(fib2)
length(unique(rownames(fib2)))
nrow(geneid)
nrow(fib2)
#transfer transcripts id to gene id and merge transcripts id

tx2gene <- data.frame(geneid)
nrow(tx2gene)
# merge transcripts id
#只提取有老鼠基因名基因名的对应的老鼠转录ID
out <- lapply(unique(tx2gene$ensembl_gene_id),
              function(x) {
                tmp       <- tx2gene[tx2gene$ensembl_gene_id == x,]
                tmp_count <- fib2[match(tmp$ensembl_transcript_id,
                                        rownames(fib2)),]
                tmp_out   <- colSums(tmp_count)
                return(tmp_out)
              })

gene_counts <- matrix(unlist(out), 
                      ncol  = ncol(fib2), 
                      byrow = T)
rownames(gene_counts) <- unique(tx2gene$ensembl_gene_id)
colnames(gene_counts) <- colnames(fib2)
head(gene_counts)
fib2 <- gene_counts
#加入附加信息
nrow(fib2)
ncol(fib2)-2

setdiff(rownames(fib2),rownames(fib1))
rownames(fib2)[! rownames(fib2) %in% rownames(fib1)]

# [1] "OR5BS1P"    "FAM237A"    "TMEM249"    "THSD8"      "HIST1H2BD"  "AL034430.1" "DPEP2NB"   
# [8] "HIST1H2BN"  "LMLN2"      "AC092442.1" "GPBP1"      "BX255925.3" "MARCOL"     "RNF212B"   
# [15] "AL032819.3" "AL353572.3" "IQCN"       "AL022312.1" "AL034430.2" "AC008481.3" "EPPK1"     
#merge fib fib2
fib2 <- fib2[match(rownames(fib1),rownames(fib2)),]
FIB <- cbind(fib1,fib2)
head(FIB)

FIB <- as.data.frame(FIB)
FIB$NAME <- rownames(FIB)
FIB$Description <- rep("na",nrow(FIB))
head(FIB)

FIB <- FIB[,c(9,10,1,2,5,6,3,4,7,8)]

write.table(FIB,file="./FIB_GSEA_1231.gct",sep = "\t",row.names = T,col.names = T,quote = F)
nrow(FIB)
ncol(FIB)-2