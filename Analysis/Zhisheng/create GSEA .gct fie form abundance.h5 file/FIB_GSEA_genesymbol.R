# Dec 31 2018
#merge 2016 and 2018 tpm data to do GSEA
#USE HUMAN GENE SYMBOL


###################################################
# FIB_slow1 FIB_slow2 FIB_slow3 FIB_slow4 FIB_fast1 FIB_fast2 FIB_fast3 FIB_fast4 
# 27627981  24313890  24632993  25713148  26875099  28225032  25204377  24867463 


#-----read 2016 data---------------
fib1 <- read.table("./raw counts data from abundance h5/FIB_1_h5.s",header = T)

fib1$FIB_slow1 <- fib1$FIB_slow1/sum(fib1$FIB_slow1)
fib1$FIB_slow2 <- fib1$FIB_slow2/sum(fib1$FIB_slow2)
fib1$FIB_fast1 <- fib1$FIB_fast1/sum(fib1$FIB_fast1)
fib1$FIB_fast2 <- fib1$FIB_fast2/sum(fib1$FIB_fast2)
fib1 <- fib1*10^6
#remove all characters behind decimal
rownames(fib1) <- gsub("\\..*", "", rownames(fib1))
head(fib1)
nrow(fib1) #117667
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
# Transcript.stable.ID     Gene.stable.ID Gene.name Gene.stable.ID.1
# 1   ENSMUST00000184707 ENSMUSG00000098914   RF02109  ENSG00000278287
# 2   ENSMUST00000157094 ENSMUSG00000087719     RPPH1  ENSG00000277209
# 3   ENSMUST00000158341 ENSMUSG00000088966   RF00402  ENSG00000199392
# 4   ENSMUST00000157139 ENSMUSG00000087764     RPPH1  ENSG00000277209
# 5   ENSMUST00000158484 ENSMUSG00000089109     RPPH1  ENSG00000277209
# 6   ENSMUST00000158050 ENSMUSG00000088675     RPPH1  ENSG00000277209
colnames(human_symbol) <- c("mouse_ensembl_transcript_id","mouse_ensembl_gene_id","human_external_gene_name","human_ensembl_gene_id")

head(fib1)
length(unique(rownames(fib1)))
nrow(human_symbol)   #90655 
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
fib2 <- read.table("./FIB_2_h5.s",header = T)
colSums(fib2)
fib2$FIB_slow3 <- fib2$FIB_slow3/sum(fib2$FIB_slow3)
fib2$FIB_slow4 <- fib2$FIB_slow4/sum(fib2$FIB_slow4)
fib2$FIB_fast3 <- fib2$FIB_fast3/sum(fib2$FIB_fast3)
fib2$FIB_fast4 <- fib2$FIB_fast4/sum(fib2$FIB_fast4)
fib2 <- fib2*10^6

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

FIB <- as.data.frame(FIB)
FIB$NAME <- rownames(FIB)
FIB$Description <- rep("na",nrow(FIB))
head(FIB)

FIB <- FIB[,c(9,10,1,2,5,6,3,4,7,8)]
colSums(FIB[,c(3,4,5,6,7,8,9,10)])

write.table(FIB,file="./FIB_GSEA_1231.gct",sep = "\t",row.names = F,col.names = T,quote = F)
nrow(FIB)
ncol(FIB)-2

FIB <- read.table("./FIB_GSEA_1231.gct",header = T)
colSums(FIB[,c(3,4,5,6,7,8,9,10)])

# FIB_slow1 FIB_slow2 FIB_slow3 FIB_slow4 FIB_fast1 FIB_fast2 FIB_fast3 FIB_fast4 
# 956105.9  955899.1  953602.4  951539.7  941166.2  945500.8  934269.2  940320.3 

