# do GSEA with 2016 data and 2018 data separately
#use gene level/trans level

#load 206 data
fib1 <- read.table("./raw counts data from abundance h5/FIB_1_h5.s",header = T)
head(fib1)

# extract transcripts id and gene id
fib1_tid <- gsub("\\..*", "", rownames(fib1))

fib1_gid <- gsub("ENSMUST.{1,20}\\|", "", rownames(fib1))
fib1_gid <- gsub("\\..*", "", fib1_gid)

head(fib1_tid);head(fib1_gid)

# build tx2gene
tx2gene <- as.data.frame(cbind(fib1_tid,fib1_gid))
colnames(tx2gene) <- c("ensembl_transcript_id","ensembl_gene_id")
head(tx2gene) 

rownames(fib1) <- gsub("\\..*", "", rownames(fib1))

# merge transcripts id

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

write.table(fib1,file="./GSEA/fib2016.tsv",row.names = T,col.names = T,sep = "\t")


fib1 <- as.data.frame(fib1)
fib1$NAME <- rownames(fib1)
fib1$Description <- rep("na",nrow(fib1))
head(fib1)
fib1 <- fib1[,c(5,6,1,2,3,4)]

write.table(fib1,file="./GSEA/fib2016.gct",row.names = F,col.names = T,sep = "\t",quote = F)

# test how many corresponding human gene symbols for mouse gene id
test <- read.table("./GSEA/fib2016.tsv",header = T)
head(test)
nrow(test)

require("biomaRt")

mart <- useMart("ENSEMBL_MART_ENSEMBL")
# listmart<-listDatasets(mart)
# View(listmart)
mart <- useDataset("mmusculus_gene_ensembl", mart)

geneid <- getBM(
  mart=mart,
  attributes=c("ensembl_transcript_id","ensembl_gene_id"),
  filter="ensembl_gene_id",
  values=rownames(test),
  uniqueRows=F)

nrow(geneid)  #130230




require("biomaRt")
mart1 = useMart("ensembl", dataset="hsapiens_gene_ensembl")
mart2 = useMart("ensembl", dataset="mmusculus_gene_ensembl") 

# get human symbol for mouse transcript id 
#只有拥有对应基因名的会提取出来
human_symbol <- getLDS(attributes=c("ensembl_gene_id"),
                       filters="ensembl_gene_id", values=rownames(test), mart=mart2,
                       attributesL=c("external_gene_name","ensembl_gene_id"), martL=mart1                         ,uniqueRows = F)

nrow(human_symbol)




# do GSEA with 2016 data and 2018 data separately
#use transcription level

#load 2016 data
fib1 <- read.table("./raw counts data from abundance h5/FIB_1_h5.s",header = T)
rownames(fib1) <- gsub("\\..*", "", rownames(fib1))
head(fib1)
nrow(fib1)

fib1 <- as.data.frame(fib1)
fib1$NAME <- rownames(fib1)
fib1$Description <- rep("na",nrow(fib1))
head(fib1)
fib1 <- fib1[,c(5,6,1,2,3,4)]
write.table(fib1,file="./GSEA/fib2016_transcription_level.gct",row.names = F,col.names = T,sep = "\t",quote = F)


#use human gene symol
fib <- read.table("./GSEA/FIB_GSEA_1231.gct",header = T)

fib <- fib[,c(1,2,3,4,7,8)]
head(fib)
nrow(fib)
write.table(fib,file="./GSEA/fib2016_genesymbol_level.gct",row.names = F,col.names = T,sep = "\t",quote = F)
