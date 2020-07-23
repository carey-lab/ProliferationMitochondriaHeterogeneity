# Jan 03 2019
#extract gene id to do PCA
# to see the number of gene id and the number of matched human gene symbol and 
# use chip to do GSEA as you have gene ensembl id

#fib
#-----read 2016 data---------------
fib1 <- read.table("./raw counts data from abundance h5/FIB_1_h5.s",header = T)
#remove all characters behind decimal
rownames(fib1) <- gsub("\\..*", "", rownames(fib1))
head(fib1)
nrow(fib1) #117667
# read 2018 data-------------------
fib2 <- read.table("./raw counts data from abundance h5/FIB_2_h5.s",header = T)
rownames(fib2) <- gsub("\\..*", "", rownames(fib2))
head(fib2)
nrow(fib2) #136535
temp <- match(rownames(fib1),rownames(fib2))
fib2 <- fib2[match(rownames(fib1),rownames(fib2)),]
nrow(fib2) #117667

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

#
nrow(geneid)  #116932
nrow(fib1)   #117667
# > nrow(geneid)
# [1] 116932
# > nrow(fib1)
# [1] 117667


# merge fib1 and fib2
FIB <- cbind(fib1,fib2)
FIB <- FIB[,c(1,2,5,6,3,4,7,8)]
head(FIB)



#transfer transcripts id to gene id and merge transcripts id
tx2gene <- data.frame(geneid)
nrow(tx2gene)

# merge transcripts id
#只提取有老鼠基因名基因名的对应的老鼠转录ID
out <- lapply(unique(tx2gene$ensembl_gene_id),
              function(x) {
                tmp       <- tx2gene[tx2gene$ensembl_gene_id == x,]
                tmp_count <- FIB[match(tmp$ensembl_transcript_id,
                                        rownames(FIB)),]
                tmp_out   <- colSums(tmp_count)
                return(tmp_out)
              })

gene_counts <- matrix(unlist(out), 
                      ncol  = ncol(FIB), 
                      byrow = T)
rownames(gene_counts) <- unique(tx2gene$ensembl_gene_id)
colnames(gene_counts) <- colnames(FIB)
head(gene_counts)
FIB <- gene_counts
#
nrow(FIB) #48178
head(FIB)
# save fib gene id counts(not normalized counts)
write.table(FIB,file="./PCA/fib gene id counts.csv",sep = ",")


#
library(DESeq2)
library(gplots)
library(RColorBrewer)
library(ggplot2)

coldata <- read.table(file="./PCA/fib gene id counts coldata.csv", header=T,sep = ',')
fib <- read.csv("./PCA/fib gene id counts.csv",header = T)
fib <- na.omit(fib)
fib <- round(fib)
nrow(fib) #48177
head(fib)
dds <- DESeqDataSetFromMatrix(fib,colData=coldata,design=~condition+year+condition:year)

dds<-DESeq(dds)
head(counts(dds))
head(counts(dds_rLogTransformed))

dds_rLogTransformed <- rlog( dds ) # regularized-log transform the data , for PCA & clustering

# result of differential expression analyse
#resCFSE = results(dds, contrast=c("condition", "slow", "fast"))

# apply PCA
plotPCA( dds_rLogTransformed , intgroup = c( "condition","year") ) 

# PCA 2016VS2018
# coldata2 <- read.table(file="./PCA/fib gene id counts coldata_2016VS2018.csv", header=T,sep = ',')
# 
# dds2 <- DESeqDataSetFromMatrix(fib,colData=coldata2,design=~condition)
# 
# dds2 <-DESeq(dds2)
# 
# dds_rLogTransformed2 <- rlog( dds2 )
# 
# #resCFSE2 = results(dds, contrast=c("condition", "slow", "fast"))
# 
# plotPCA( dds_rLogTransformed2 , intgroup = c( "condition") ) 




#ESC
#-----read 2016 data---------------
ESC1 <- read.table("./raw counts data from abundance h5/ESC_1_h5.s",header = T)
#remove all characters behind decimal
rownames(ESC1) <- gsub("\\..*", "", rownames(ESC1))
head(ESC1)
nrow(ESC1) #117667
# read 2018 data 1-------------------
ESC2_1 <- read.table("./raw counts data from abundance h5/ESC_2_1_h5.s",header = T)
rownames(ESC2_1) <- gsub("\\..*", "", rownames(ESC2_1))
nrow(ESC2_1) #136535
# read 2018 data 2-------------------
ESC2_2 <- read.table("./raw counts data from abundance h5/ESC_2_2_h5.s",header = T)
rownames(ESC2_2) <- gsub("\\..*", "", rownames(ESC2_2))
#merge ESC2_1 ESC2_2 
ESC2 <- cbind(ESC2_1,ESC2_2)

head(ESC2)
nrow(ESC2) #136535
temp <- match(rownames(ESC1),rownames(ESC2))
ESC2 <- ESC2[match(rownames(ESC1),rownames(ESC2)),]
nrow(ESC2) #117667


# merge ESC1 and ESC2
ESC <- cbind(ESC1,ESC2)
ESC <- ESC[,c(1,2,5,6,3,4,7,8)]
head(ESC)

# check whether rownames(ESC)=rownames(fib1), if so we do not need to use biomart again, if not we shell do
setdiff(rownames(ESC),rownames(fib1))



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
  values=rownames(ESC1),
  uniqueRows=F)

head(geneid)
head(ESC1)
length(unique(rownames(ESC1)))

#
nrow(geneid)  #116932
nrow(ESC1)   #117667
# > nrow(geneid)
# [1] 116932
# > nrow(ESC1)
# [1] 117667


#transfer transcripts id to gene id and merge transcripts id
tx2gene <- data.frame(geneid)
nrow(tx2gene)

# merge transcripts id
#只提取有老鼠基因名基因名的对应的老鼠转录ID
out <- lapply(unique(tx2gene$ensembl_gene_id),
              function(x) {
                tmp       <- tx2gene[tx2gene$ensembl_gene_id == x,]
                tmp_count <- ESC[match(tmp$ensembl_transcript_id,
                                       rownames(ESC)),]
                tmp_out   <- colSums(tmp_count)
                return(tmp_out)
              })

gene_counts <- matrix(unlist(out), 
                      ncol  = ncol(ESC), 
                      byrow = T)
rownames(gene_counts) <- unique(tx2gene$ensembl_gene_id)
colnames(gene_counts) <- colnames(ESC)
head(gene_counts)
ESC <- gene_counts
#
nrow(ESC) #48178
head(ESC)
# save ESC gene id counts(not normalized counts)
write.table(ESC,file="./PCA/ESC gene id counts.csv",sep = ",")


#
library(DESeq2)
library(gplots)
library(RColorBrewer)
library(ggplot2)

coldata <- read.table(file="./PCA/ESC gene id counts coldata.csv", header=T,sep = ',')
ESC <- read.csv("./PCA/ESC gene id counts.csv",header = T)
ESC <- na.omit(ESC)
ESC <- round(ESC)
nrow(ESC) #48177
head(ESC)
dds <- DESeqDataSetFromMatrix(ESC,colData=coldata,design=~condition+year+condition:year)

dds<-DESeq(dds)
head(counts(dds))


dds_rLogTransformed <- rlog( dds ) # regularized-log transform the data , for PCA & clustering

# result of differential expression analyse
#resCFSE = results(dds, contrast=c("condition", "slow", "fast"))

# apply PCA
plotPCA( dds_rLogTransformed , intgroup = c( "condition","year") ) 

# # PCA 2016VS2018
# coldata2 <- read.table(file="./PCA/ESC gene id counts coldata_2016VS2018.csv", header=T,sep = ',')
# 
# dds2 <- DESeqDataSetFromMatrix(ESC,colData=coldata2,design=~condition)
# 
# dds2 <-DESeq(dds2)
# 
# dds_rLogTransformed2 <- rlog( dds2 )
# 
# #resCFSE2 = results(dds, contrast=c("condition", "slow", "fast"))
# 
# plotPCA( dds_rLogTransformed2 , intgroup = c( "condition") ) 
# #PCA 2016VS2018 discard 2018 ESC fast
# coldata3 <- read.table(file="./PCA/ESC gene id counts coldata_2016VS2018_discard fast2018.csv", header=T,sep = ',')
# 
# ESC3 <- ESC[,c(1,2,3,4,5,6)]
# 
# dds3 <- DESeqDataSetFromMatrix(ESC3,colData=coldata3,design=~condition)
# 
# dds3 <-DESeq(dds3)
# 
# dds_rLogTransformed3 <- rlog( dds3 )
# 
# #resCFSE2 = results(dds, contrast=c("condition", "slow", "fast"))
# 
# plotPCA( dds_rLogTransformed3 , intgroup = c( "condition") ) 
# 
# 
# # what if I exchange 2018 solw with 2016 fast
# 
# coldata4 <- read.table(file="./PCA/ESC gene id counts coldata_2016VS2018.csv", header=T,sep = ',')
# 
# head(ESC4)
# 
# ESC4 <- ESC[,c(1,2,5,6,3,4,7,8)]
# colnames(ESC4) <- c("ESC_slow1", "ESC_slow2", "ESC_slow3", "ESC_slow4", "ESC_fast1", "ESC_fast2", "ESC_fast3","ESC_fast4")
# 
# 
# dds4 <- DESeqDataSetFromMatrix(ESC4,colData=coldata4,design=~condition)
# 
# dds4 <-DESeq(dds4)
# 
# dds_rLogTransformed4 <- rlog( dds4 )
# 
# #resCFSE2 = results(dds, contrast=c("condition", "slow", "fast"))
# 
# plotPCA( dds_rLogTransformed4 , intgroup = c( "condition") ) 
# 

# PCA ORIGINAL 2016 SLOW VS 2016 FAST

coldata5 <- read.table(file="./PCA/ESC gene id counts coldata_ORIGINAL 2016 SLOW VS 2016 FAST.csv", header=T,sep = ',')


ESC5 <- ESC[,c(1,2,5,6)]
head(ESC5)


dds5 <- DESeqDataSetFromMatrix(ESC5,colData=coldata5,design=~condition)

dds5 <-DESeq(dds5)

dds_rLogTransformed5 <- rlog( dds5 )

#resCFSE2 = results(dds, contrast=c("condition", "slow", "fast"))

plotPCA( dds_rLogTransformed5 , intgroup = c( "condition") ) 

#PCA ORIGINAL 2018 SLOW VS 2018 fast

coldata6 <- read.table(file="./PCA/ESC gene id counts coldata_ORIGINAL 2018 SLOW VS 2018 FAST.csv", header=T,sep = ',')


ESC6 <- ESC[,c(3,4,7,8)]
head(ESC6)


dds6 <- DESeqDataSetFromMatrix(ESC6,colData=coldata6,design=~condition)

dds6 <-DESeq(dds6)

dds_rLogTransformed6 <- rlog( dds6 )

#resCFSE2 = results(dds, contrast=c("condition", "slow", "fast"))

plotPCA( dds_rLogTransformed6 , intgroup = c( "condition") ) 


#因为2016的行数和2018不同， 所以2016和2018不可能搞混。



# save all the image
for(d in dev.list()) {
  dev.set(d)
  Name = paste("Image", d, ".jpg", sep="")
  dev.copy(jpeg, Name)
  dev.off()
}

