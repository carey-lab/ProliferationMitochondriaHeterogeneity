#Dec 31 2018
# extract tpm data and counts data from h5 file
# 
# library(rhdf5)
# fib_hh <- h5read("G:/proliferation/DataCFSE/Fibro_HH_28879_GTGAAA_abundance.h5","fib_hh")
# h5ls("G:/proliferation/DataCFSE/Fibro_HH_28879_GTGAAA_abundance.h5")
# 
# 
# h5dump ("abundance.h5",
#         recursive = TRUE,
#         load = TRUE,
#         all = FALSE,
#         index_type = h5default("H5_INDEX"),
#         order = h5default("H5_ITER"), native = FALSE)
# 
# 
# 
# h5ls("G:/proliferation/DataCFSE/Fibro_HH_28879_GTGAAA_abundance.h5",
#         recursive = TRUE,
#         all = FALSE, 
#         datasetinfo = TRUE,
#         index_type = h5default("H5_INDEX"),
#         order = h5default("H5_ITER"), native = FALSE)


###################import data########################
#------------------
library(tximport)
dir <- print(getwd())
##########################TxDb.Mmusculus.UCSC.mm10.knownGene
# library(TxDb.Mmusculus.UCSC.mm10.knownGene)
# txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
# k <- keys(txdb, keytype = "TXNAME")
# tx2gene2 <- select(txdb, k, "GENEID", "TXNAME")
# head(tx2gene2)
# library(EnsDb.Hsapiens.v86)
# head( keys(EnsDb.Hsapiens.v86, keytype="TXID") )
###################ts2gene-----------------------
library(EnsDb.Mmusculus.v79)
EnsDb <- EnsDb.Mmusculus.v79
#Allowed choices are: 'ENTREZID', 'EXONID', 'GENEBIOTYPE', 'GENEID', 'GENENAME', 'PROTDOMID', 'PROTEINDOMAINID', 'PROTEINDOMAINSOURCE', 'PROTEINID', 'SEQNAME', 'SEQSTRAND', 'SYMBOL', 'TXBIOTYPE', 'TXID', 'TXNAME', 'UNIPROTID'.
k <- keys(EnsDb, keytype = "TXNAME")
tx2gene <- select(EnsDb, k, "GENEID", "TXNAME")
tx2gene <- tx2gene[,1:2]

columns(EnsDb)
keytypes(EnsDb)


###############FIB 1 DATA---------------------

samples1 <- read.table(file.path(dir, "FIB_slow2fast_1.tsv"), header = TRUE)

files1 <- file.path(dir, samples1$run, "abundance.h5")

file.exists(files1)
names(files1) <- c("FIB_slow1","FIB_slow2","FIB_fast1","FIB_fast2"
                  )
FIB_1 <- tximport(files1, type = "kallisto", txOut = TRUE,tx2gene=tx2gene,ignoreAfterBar = TRUE)
head(FIB_1$counts)
FIB_1_h5 <- FIB_1$counts
write.table(FIB_1_h5,file="./FIB_1_h5.s",sep="\t",col.names = T,row.names = T)

#deseq------use deseq2 to remove decimal in the counts

samples1$condition <- factor(c("A","A","B","B"))
rownames(samples1) <- c("FIB_slow1","FIB_slow2","FIB_fast1","FIB_fast2"
)


library("DESeq2")
FIB_1_deseq <- DESeqDataSetFromTximport(FIB_1,
                                        colData = samples1,
                                        design = ~ condition)

# save the unnormalized data from deseq
FIB_1_deseq_unnormalize <- counts(FIB_1_deseq, normalized=F)
write.table(FIB_1_deseq_unnormalize,file="./FIB_1_deseq_unnormalize.tsv",sep = "\t",row.names = T,col.names = T)


###############FIB 2 DATA---------------------

samples2 <- read.table(file.path(dir, "FIB_slow2fast_2.tsv"), header = TRUE)

files2 <- file.path(dir, samples2$run, "abundance.h5")

file.exists(files2)
names(files2) <- c("FIB_slow3","FIB_slow4","FIB_fast3","FIB_fast4"
)
FIB_2 <- tximport(files2, type = "kallisto", txOut = TRUE,tx2gene=tx2gene,ignoreAfterBar = TRUE)
head(FIB_2$counts)
FIB_2_h5 <- FIB_2$counts
write.table(FIB_2_h5,file="./FIB_2_h5.s",sep="\t",col.names = T,row.names = T)

#deseq------

# the number of gene in 2016 is different from 2018
# how to merge them togather?
#try the tx2gene
samples2$condition <- factor(c("A","A","B","B"))
rownames(samples2) <- c("FIB_slow3","FIB_slow4","FIB_fast3","FIB_fast4"
)


library("DESeq2")
FIB_2_deseq <- DESeqDataSetFromTximport(FIB_2,
                                              colData = samples2,
                                              design = ~ condition)

# save the unnormalized data from deseq
FIB_2_deseq_unnormalize <- counts(FIB_2_deseq, normalized=F)
write.table(FIB_2_deseq_unnormalize,file="./FIB_2_deseq_unnormalize.tsv",sep = "\t",row.names = T,col.names = T)

# #deseq normalize
# FIB_2_deseq2 <- DESeq(FIB_2_deseq)
# #get result
# res <- results(FIB_ESC_2_1_deseq)
# res


#################################################################
###############ESC 1 DATA---------------------


########################################

samples1 <- read.table(file.path(dir, "ESC_slow2fast_1.tsv"), header = TRUE)

files1 <- file.path(dir, samples1$run, "abundance.h5")

file.exists(files1)
names(files1) <- c("ESC_slow1","ESC_slow2","ESC_fast1","ESC_fast2"
)
ESC_1 <- tximport(files1, type = "kallisto", txOut = TRUE,tx2gene=tx2gene,ignoreAfterBar = TRUE)
head(ESC_1$counts)
ESC_1_h5 <- ESC_1$counts
write.table(ESC_1_h5,file="./ESC_1_h5.s",sep="\t",col.names = T,row.names = T)

#deseq------

samples1$condition <- factor(c("A","A","B","B"))
rownames(samples1) <- c("ESC_slow1","ESC_slow2","ESC_fast1","ESC_fast2"
)


library("DESeq2")
ESC_1_deseq <- DESeqDataSetFromTximport(ESC_1,
                                        colData = samples1,
                                        design = ~ condition)

# save the unnormalized data from deseq
ESC_1_deseq_unnormalize <- counts(ESC_1_deseq, normalized=F)
write.table(ESC_1_deseq_unnormalize,file="./ESC_1_deseq_unnormalize.tsv",sep = "\t",row.names = T,col.names = T)


###############ESC 2-1 DATA---------------------

samples2 <- read.table(file.path(dir, "ESC_slow2fast_2_1.tsv"), header = TRUE)

files2 <- file.path(dir, samples2$run, "abundance.h5")

file.exists(files2)
names(files2) <- c("ESC_slow3","ESC_slow4")

ESC_2_1 <- tximport(files2, type = "kallisto", txOut = TRUE,tx2gene=tx2gene,ignoreAfterBar = TRUE)
head(ESC_2_1$counts)
ESC_2_1_h5 <- ESC_2_1$counts
write.table(ESC_2_1_h5,file="./ESC_2_1_h5.s",sep="\t",col.names = T,row.names = T)

#deseq------

# the number of gene in 2016 is different from 2018
# how to merge them togather?
#try the tx2gene
samples2$condition <- factor(c("A","B"))
rownames(samples2) <- c("ESC_slow3","ESC_slow4"
)


library("DESeq2")
ESC_2_1_deseq <- DESeqDataSetFromTximport(ESC_2_1,
                                        colData = samples2,
                                        design = ~ condition)

# save the unnormalized data from deseq
ESC_2_1_deseq_unnormalize <- counts(ESC_2_1_deseq, normalized=F)
write.table(ESC_2_1_deseq_unnormalize,file="./ESC_2_1_deseq_unnormalize.tsv",sep = "\t",row.names = T,col.names = T)


###############ESC 2-2 DATA---------------------

samples2 <- read.table(file.path(dir, "ESC_slow2fast_2_2.tsv"), header = TRUE)

files2 <- file.path(dir, samples2$run, "abundance.h5")

file.exists(files2)
names(files2) <- c("ESC_fast3","ESC_fast4"
)
ESC_2_2 <- tximport(files2, type = "kallisto", txOut = TRUE,tx2gene=tx2gene,ignoreAfterBar = TRUE)
head(ESC_2_2$counts)
ESC_2_2_h5 <- ESC_2_2$counts
write.table(ESC_2_2_h5,file="./ESC_2_2_h5.s",sep="\t",col.names = T,row.names = T)

#deseq------

# the number of gene in 2016 is different from 2018
# how to merge them togather?
#try the tx2gene
samples2$condition <- factor(c("A","B"))
rownames(samples2) <- c("ESC_fast3","ESC_fast4"
)


library("DESeq2")
ESC_2_2_deseq <- DESeqDataSetFromTximport(ESC_2_2,
                                        colData = samples2,
                                        design = ~ condition)

# save the unnormalized data from deseq
ESC_2_2_deseq_unnormalize <- counts(ESC_2_2_deseq, normalized=F)
write.table(ESC_2_2_deseq_unnormalize,file="./ESC_2_2_deseq_unnormalize.tsv",sep = "\t",row.names = T,col.names = T)




















###############--------------END---------------------



###########################2018data----------------
samples2_1 <- read.table(file.path(dir, "samples_fast2slow2_1.tsv"), header = TRUE)

files2_1 <- file.path(dir, samples2_1$run, "abundance.h5")

file.exists(files2_1)
names(files2_1) <- c("ES_slow3","ES_fast3","ES_slow4","ES_fast4"
                  )
FIB_ESC_2_1 <- tximport(files2_1, type = "kallisto", txOut = TRUE)
FIB_ESC_2_1_h5 <- FIB_ESC_2_1$counts
write.table(FIB_ESC_2_1_h5,file="./FIB_ESC_2_1_h5.s",sep="\t",col.names = T,row.names = T)
#######-----------

samples2_2 <- read.table(file.path(dir, "samples_fast2slow2_2.tsv"), header = TRUE)

files2_2 <- file.path(dir, samples2_2$run, "abundance.h5")

file.exists(files2_2)
names(files2_2) <- c("FIB_slow3","FIB_fast3","FIB_slow4","FIB_fast4"
                     )
FIB_ESC_2_2 <- tximport(files2_2, type = "kallisto", txOut = TRUE)
head(FIB_ESC_2_2$counts)

#######------------deseq2  FIB_ESC_2_1_deseq-------------------
# the number of gene in 2016 is different from 2018
# how to merge them togather?
#try the tx2gene
samples2_1$condition <- factor(c("A","B","A","B"))
rownames(samples2_1) <- c("ES_slow3","ES_fast3","ES_slow4","ES_fast4"
)


library("DESeq2")
FIB_ESC_2_1_deseq <- DESeqDataSetFromTximport(FIB_ESC_2_1,
                                   colData = samples2_1,
                                   design = ~ condition)

FIB_ESC_2_1_deseq <- DESeq(FIB_ESC_2_1_deseq)
res <- results(FIB_ESC_2_1_deseq)
res











