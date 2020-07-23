# Jan 1 2019
#merge 2016 and 2018 counts data to do DEA
#USE HUMAN GENE SYMBOL


###################################################
# 1. merge tarscrpts level counts
# 2. do deseq2
# 3. transfer


#-----read 2016 data---------------
fib1 <- read.table("./FIB_1_deseq_unnormalize.tsv",header = T)
colSums(fib1)


#remove all characters behind decimal
rownames(fib1) <- gsub("\\..*", "", rownames(fib1))
head(fib1)
nrow(fib1)
#get correspondent human gene symbol
require("biomaRt")
mart1 = useMart("ensembl", dataset="hsapiens_gene_ensembl")
mart2 = useMart("ensembl", dataset="mmusculus_gene_ensembl") 

# get human symbol for mouse transcript id 
#只有拥有对应基因名的会提取出来
human_symbol <- getLDS(attributes=c("ensembl_transcript_id","ensembl_gene_id"),
                       filters="ensembl_transcript_id", values=rownames(fib1), mart=mart2,
                       attributesL=c("external_gene_name","ensembl_gene_id"), martL=mart1                         ,uniqueRows = F)

head(human_symbol)
colnames(human_symbol) <- c("mouse_ensembl_transcript_id","mouse_ensembl_gene_id","human_external_gene_name","human_ensembl_gene_id")

head(fib1)
length(unique(rownames(fib1)))
nrow(human_symbol)
#transfer transcripts id to gene id and merge transcripts id

tx2gene <- data.frame(human_symbol[,c(1,3)])
nrow(tx2gene)
# merge transcripts id
#只提取有人类基因名的对应的老鼠转录ID
out <- lapply(unique(tx2gene$human_external_gene_name),
              function(x) {
                tmp       <- tx2gene[tx2gene$human_external_gene_name == x,]
                tmp_count <- fib1[match(tmp$mouse_ensembl_transcript_id,
                                        rownames(fib1)),]
                tmp_out   <- colSums(tmp_count)
                return(tmp_out)
              })

gene_counts <- matrix(unlist(out), 
                      ncol  = ncol(fib1), 
                      byrow = T)
rownames(gene_counts) <- unique(tx2gene$human_external_gene_name)
colnames(gene_counts) <- colnames(fib1)
head(gene_counts)
fib1 <- gene_counts
#加入附加信息
nrow(fib1)
ncol(fib1)-2

#-----read 2018 data---------------
#the colsums of 2018 data are not 10^6
#transfer it to 10^6
fib2 <- read.table("./FIB_2_deseq_unnormalize.tsv",header = T)
colSums(fib2)

#remove all characters behind decimal
rownames(fib2) <- gsub("\\..*", "", rownames(fib2))
head(fib2)

#get correspondent human gene symbol
require("biomaRt")
mart1 = useMart("ensembl", dataset="hsapiens_gene_ensembl")
mart2 = useMart("ensembl", dataset="mmusculus_gene_ensembl") 

# get human symbol for mouse transcript id 
#只有拥有对应基因名的会提取出来
human_symbol <- getLDS(attributes=c("ensembl_transcript_id","ensembl_gene_id"),
                       filters="ensembl_transcript_id", values=rownames(fib2), mart=mart2,
                       attributesL=c("external_gene_name","ensembl_gene_id"), martL=mart1                         ,uniqueRows = F)

head(human_symbol)
colnames(human_symbol) <- c("mouse_ensembl_transcript_id","mouse_ensembl_gene_id","human_external_gene_name","human_ensembl_gene_id")

head(fib2)
length(unique(rownames(fib2)))
nrow(human_symbol)
nrow(fib2)
#transfer transcripts id to gene id and merge transcripts id
nrow(tx2gene)
tx2gene <- data.frame(human_symbol[,c(1,3)])

# merge transcripts id
#只提取有人类基因名的对应的老鼠转录ID
out <- lapply(unique(tx2gene$human_external_gene_name),
              function(x) {
                tmp       <- tx2gene[tx2gene$human_external_gene_name == x,]
                tmp_count <- fib2[match(tmp$mouse_ensembl_transcript_id,
                                        rownames(fib2)),]
                tmp_out   <- colSums(tmp_count)
                return(tmp_out)
              })

gene_counts <- matrix(unlist(out), 
                      ncol  = ncol(fib2), 
                      byrow = T)
rownames(gene_counts) <- unique(tx2gene$human_external_gene_name)
colnames(gene_counts) <- colnames(fib2)
head(gene_counts)
fib2 <- gene_counts
#加入附加信息
nrow(fib2)
ncol(fib2)-2
nrow(fib1)
setdiff(rownames(fib2),rownames(fib1))
rownames(fib2)[! rownames(fib2) %in% rownames(fib1)]

# [1] "OR5BS1P"    "FAM237A"    "TMEM249"    "THSD8"      "HIST1H2BD"  "AL034430.1" "DPEP2NB"   
# [8] "HIST1H2BN"  "LMLN2"      "AC092442.1" "GPBP1"      "BX255925.3" "MARCOL"     "RNF212B"   
# [15] "AL032819.3" "AL353572.3" "IQCN"       "AL022312.1" "AL034430.2" "AC008481.3" "EPPK1"     
#merge fib fib2
fib2 <- fib2[match(rownames(fib1),rownames(fib2)),]
FIB <- cbind(fib1,fib2)

head(FIB)
colSums(fib2)
colSums(fib1)
colSums(FIB)



FIB <- FIB[,c(1,2,5,6,3,4,7,8)]
colSums(FIB[,c(1,2,3,4,5,6,7,8)])

write.table(FIB,file="./FIB_counts_0101.tsv",sep = ",",row.names = T,col.names = T,quote = F)

#############################

# apply deseq2################




library("DESeq2")

mycount<-read.table(file='FIB_counts_0101.tsv', header=T,sep = ',')
# mycount$gene <- rownames(mycount)
# head(mycount,10)
# mycount <- mycount[,c(9,1,2,3,4,5,6,7,8)]
# mycount <- as.matrix(mycount,rownames=F)
# rownames(mycount) <- mycount[,1]
# mycount <- mycount[,-1]
# mycount <- as.data.frame(mycount)
coldata<-read.table(file="./deceq2_coldata.tsv", header=T,sep = ',')

head(coldata)

head(mycount)

dds <- DESeqDataSetFromMatrix(mycount,colData=coldata,design=~condition)

dds<-DESeq(dds)


res = results(dds, contrast=c("condition", "slow", "fast"))

res = res[order(res$padj),]

head(res)

summary(res)

table(res$padj<0.05 & abs(res$log2FoldChange) > 1)


write.csv(res,file="./FIB_DEseq2_result.csv")


#p adj <0.05  significant
res <- na.omit(res)
res <- res[res$padj<0.05,]

write.csv(res,file="./FIB_DEseq2_significant_result.csv")

#significant and abs(log2fold change)>1

res <- res[res$padj<0.05 & abs(res$log2FoldChange) > 1,]
head(res)
write.table(rownames(res),file="./FIB_DEseq2_significant_result_GO.tsv",row.names = F,col.names = F,quote = F)

#significant and log2fold change>1
res1 <- res[res$padj<0.05 & res$log2FoldChange > 1,]
head(res1)
write.table(rownames(res1),file="./FIB_DEseq2_significant_result_GO_1.tsv",row.names = F,col.names = F,quote = F)

#significant and log2fold change<1
res2 <- res[res$padj<0.05 & res$log2FoldChange < -1,]
head(res2)
write.table(rownames(res2),file="./FIB_DEseq2_significant_result_GO_2.tsv",row.names = F,col.names = F,quote = F)
# volcano plot

library("ggplot2")

mydata=read.table('./FIB_DEseq2_result.csv',header=T,sep = ',')

head(mydata)

diff_gene_deseq2 <-subset(mydata, padj < 0.05 & abs(log2FoldChange) > 1)

dim(diff_gene_deseq2)

mydata$Significant<-'Not change'

mydata$Significant[mydata$padj< 0.05 & mydata$log2FoldChange>1]<-'Up'

mydata$Significant[mydata$padj < 0.05 & mydata$log2FoldChange<(-1)]<-'Down'

head(mydata)

ggplot(data=mydata,aes(x=log2FoldChange,y=-log10(padj)))+
  
  geom_point(size=0.7,aes(color=Significant))+
  
  geom_hline(yintercept=1.3,linetype=4)+geom_vline(xintercept=c(-1,1),linetype=4)+
  
  coord_cartesian(xlim=c(-7,7),ylim=c(0,30))+
  
  scale_color_manual(values=c("blue","grey","red"),guide=guide_legend(title="Significant"))+
  
  labs(x="Log2FoldChange",y="-Log10(padj)")+
  
  theme_bw(base_size=15)





TEMP <- read.table("./FIB_counts_0101.tsv",head=T,sep = ",")
write.table(rownames(TEMP),file="./GO_backgroud_gene",row.names = F,col.names = F,quote = F)
nrow(TEMP)

















