# Dec 31 2018
#merge 2016 and 2018 tpm data to do GSEA
#USE HUMAN GENE SYMBOL


###################################################




#-----read 2016 data---------------
ESC1 <- read.table("./raw counts data from abundance h5/ESC_1_h5.s",header = T)
head(ESC1)
ESC1$ESC_slow1 <- ESC1$ESC_slow1/sum(ESC1$ESC_slow1)
ESC1$ESC_slow2 <- ESC1$ESC_slow2/sum(ESC1$ESC_slow2)
ESC1$ESC_fast1 <- ESC1$ESC_fast1/sum(ESC1$ESC_fast1)
ESC1$ESC_fast2 <- ESC1$ESC_fast2/sum(ESC1$ESC_fast2)
ESC1 <- ESC1*10^6
colSums(ESC1)
#remove all characters behind decimal
rownames(ESC1) <- gsub("\\..*", "", rownames(ESC1))
head(ESC1)
nrow(ESC1)
#get correspondent human gene symbol
require("biomaRt")
mart1 = useMart("ensembl", dataset="hsapiens_gene_ensembl")
mart2 = useMart("ensembl", dataset="mmusculus_gene_ensembl") 

# get human symbol for mouse transcript id 
#只有拥有对应基因名的会提取出来
human_symbol <- getLDS(attributes=c("ensembl_transcript_id","ensembl_gene_id"),
                       filters="ensembl_transcript_id", values=rownames(ESC1), mart=mart2,
                       attributesL=c("external_gene_name","ensembl_gene_id"), martL=mart1                         ,uniqueRows = F)

head(human_symbol)
colnames(human_symbol) <- c("mouse_ensembl_transcript_id","mouse_ensembl_gene_id","human_external_gene_name","human_ensembl_gene_id")

head(ESC1)
length(unique(rownames(ESC1)))
nrow(human_symbol)
#transfer transcripts id to gene id and merge transcripts id

tx2gene <- data.frame(human_symbol[,c(1,3)])
nrow(tx2gene)
# merge transcripts id, sum transcrpts level to human gene level
#只提取有人类基因名的对应的老鼠转录ID
out <- lapply(unique(tx2gene$human_external_gene_name),
              function(x) {
                tmp       <- tx2gene[tx2gene$human_external_gene_name == x,]
                tmp_count <- ESC1[match(tmp$mouse_ensembl_transcript_id,
                                        rownames(ESC1)),]
                tmp_out   <- colSums(tmp_count)
                return(tmp_out)
              })

gene_counts <- matrix(unlist(out), 
                      ncol  = ncol(ESC1), 
                      byrow = T)
rownames(gene_counts) <- unique(tx2gene$human_external_gene_name)
colnames(gene_counts) <- colnames(ESC1)
head(gene_counts) 
# ESC_slow1 ESC_slow2 ESC_fast1 ESC_fast2
# RF02109   0.0000000 0.0000000 0.0000000 0.0000000
# RPPH1     0.6835327 0.5100836 0.3943935 0.4039823
# RF00402   0.0000000 0.0000000 0.0000000 0.0000000
# RNU6-581P 0.0000000 0.0000000 0.0000000 0.0000000
# RNU1-6P   0.0000000 0.0000000 0.0000000 0.0000000
# RNU6-338P 0.0000000 0.0000000 0.0000000 0.0000000
ESC1 <- gene_counts
#加入附加信息
nrow(ESC1)
ncol(ESC1)-2

#-----read 2018 data---------------ESC_2_1-----------------------
#the colsums of 2018 data are not 10^6
#transfer it to 10^6
#merge the row data first, it will save time

ESC2_1 <- read.table("./ESC_2_1_h5.s",header = T)
head(ESC2_1)
ESC2_1$ESC_slow3 <- ESC2_1$ESC_slow3/sum(ESC2_1$ESC_slow3)
ESC2_1$ESC_slow4 <- ESC2_1$ESC_slow4/sum(ESC2_1$ESC_slow4)
ESC2_1 <- ESC2_1*10^6
colSums(ESC2_1)
#remove all characters behind decimal
rownames(ESC2_1) <- gsub("\\..*", "", rownames(ESC2_1))
head(ESC2_1)
nrow(ESC2_1)
#get correspondent human gene symbol
require("biomaRt")
mart1 = useMart("ensembl", dataset="hsapiens_gene_ensembl")
mart2 = useMart("ensembl", dataset="mmusculus_gene_ensembl") 

# get human symbol for mouse transcript id 
#只有拥有对应基因名的会提取出来
human_symbol <- getLDS(attributes=c("ensembl_transcript_id","ensembl_gene_id"),
                       filters="ensembl_transcript_id", values=rownames(ESC2_1), mart=mart2,
                       attributesL=c("external_gene_name","ensembl_gene_id"), martL=mart1                         ,uniqueRows = F)

head(human_symbol)
colnames(human_symbol) <- c("mouse_ensembl_transcript_id","mouse_ensembl_gene_id","human_external_gene_name","human_ensembl_gene_id")

head(ESC2_1)
length(unique(rownames(ESC2_1)))
nrow(human_symbol)
#transfer transcripts id to gene id and merge transcripts id

tx2gene <- data.frame(human_symbol[,c(1,3)])
nrow(tx2gene)
# merge transcripts id
#只提取有人类基因名的对应的老鼠转录ID
out <- lapply(unique(tx2gene$human_external_gene_name),
              function(x) {
                tmp       <- tx2gene[tx2gene$human_external_gene_name == x,]
                tmp_count <- ESC2_1[match(tmp$mouse_ensembl_transcript_id,
                                          rownames(ESC2_1)),]
                tmp_out   <- colSums(tmp_count)
                return(tmp_out)
              })

gene_counts <- matrix(unlist(out), 
                      ncol  = ncol(ESC2_1), 
                      byrow = T)
rownames(gene_counts) <- unique(tx2gene$human_external_gene_name)
colnames(gene_counts) <- colnames(ESC2_1)
head(gene_counts)
# ESC_slow3 ESC_slow4
# RF00493   0.0000000 0.0000000
# RPPH1     0.3224638 0.6386664
# RNU6-337P 0.0000000 0.0000000
# RNU6-47P  0.0000000 0.0000000
# RF02198   0.0000000 0.0000000
# RNU6-641P 0.0000000 0.0000000
ESC2_1 <- gene_counts
#加入附加信息
nrow(ESC2_1)
ncol(ESC2_1)-2






#########################################
############ESC_2_2

ESC2_2 <- read.table("./ESC_2_2_h5.s",header = T)
head(ESC2_2)
ESC2_2$ESC_fast3 <- ESC2_2$ESC_fast3/sum(ESC2_2$ESC_fast3)
ESC2_2$ESC_fast4 <- ESC2_2$ESC_fast4/sum(ESC2_2$ESC_fast4)
ESC2_2 <- ESC2_2*10^6
colSums(ESC2_2)
#remove all characters behind decimal
rownames(ESC2_2) <- gsub("\\..*", "", rownames(ESC2_2))
head(ESC2_2)
nrow(ESC2_2)
#get correspondent human gene symbol
require("biomaRt")
mart1 = useMart("ensembl", dataset="hsapiens_gene_ensembl")
mart2 = useMart("ensembl", dataset="mmusculus_gene_ensembl") 

# get human symbol for mouse transcript id 
#只有拥有对应基因名的会提取出来
human_symbol <- getLDS(attributes=c("ensembl_transcript_id","ensembl_gene_id"),
                       filters="ensembl_transcript_id", values=rownames(ESC2_2), mart=mart2,
                       attributesL=c("external_gene_name","ensembl_gene_id"), martL=mart1                         ,uniqueRows = F)

head(human_symbol)
colnames(human_symbol) <- c("mouse_ensembl_transcript_id","mouse_ensembl_gene_id","human_external_gene_name","human_ensembl_gene_id")

head(ESC2_2)
length(unique(rownames(ESC2_2)))
nrow(human_symbol)
#transfer transcripts id to gene id and merge transcripts id

tx2gene <- data.frame(human_symbol[,c(1,3)])
nrow(tx2gene)
# merge transcripts id
#只提取有人类基因名的对应的老鼠转录ID
out <- lapply(unique(tx2gene$human_external_gene_name),
              function(x) {
                tmp       <- tx2gene[tx2gene$human_external_gene_name == x,]
                tmp_count <- ESC2_2[match(tmp$mouse_ensembl_transcript_id,
                                          rownames(ESC2_2)),]
                tmp_out   <- colSums(tmp_count)
                return(tmp_out)
              })

gene_counts <- matrix(unlist(out), 
                      ncol  = ncol(ESC2_2), 
                      byrow = T)
rownames(gene_counts) <- unique(tx2gene$human_external_gene_name)
colnames(gene_counts) <- colnames(ESC2_2)
head(gene_counts)
# ESC_fast3  ESC_fast4
# RF00493   0.0000000 0.07445192
# RPPH1     0.5891469 0.39960679
# RNU6-337P 0.0000000 0.00000000
# RNU6-47P  0.0000000 0.17694623
# RF02198   0.0000000 0.00000000
# RNU6-641P 0.0000000 0.00000000
ESC2_2 <- gene_counts
#加入附加信息
nrow(ESC2_2)
ncol(ESC2_2)-2


# merge ESC_2_1 AND ESC_2_2
ESC2 <- cbind(ESC2_1,ESC2_2)



setdiff(rownames(ESC2),rownames(ESC1))
rownames(ESC2)[! rownames(ESC2) %in% rownames(ESC1)]

# [1] "FAM237A"    "OR5BS1P"    "THSD8"      "TMEM249"    "MARCOL"     "GPBP1"      "RNF212B"   
# [8] "AL353572.3" "AC092442.1" "HIST1H2BN"  "LMLN2"      "AL034430.1" "BX255925.3" "IQCN"      
# [15] "DPEP2NB"    "AL032819.3" "HIST1H2BD"  "AL034430.2" "AC008481.3" "EPPK1"      "AL022312.1"     
#merge ESC ESC2
ESC2 <- ESC2[match(rownames(ESC1),rownames(ESC2)),]
ESC <- cbind(ESC1,ESC2)

head(ESC2)
colSums(ESC2)
colSums(ESC1)
colSums(ESC)

ESC <- as.data.frame(ESC)
ESC$NAME <- rownames(ESC)
ESC$Description <- rep("na",nrow(ESC))
head(ESC)

ESC <- ESC[,c(9,10,1,2,5,6,3,4,7,8)]
colSums(ESC[,c(3,4,5,6,7,8,9,10)])

write.table(ESC,file="./ESC_test.gct",sep = "\t",row.names = F,col.names = T,quote = F)
nrow(ESC)
ncol(ESC)-2

ESC <- read.table("./ESC_GSEA_1231.gct",header = T)
colSums(ESC[,c(3,4,5,6,7,8,9,10)])

# ESC_slow1 ESC_slow2 ESC_slow3 ESC_slow4 ESC_fast1 ESC_fast2 ESC_fast3 ESC_fast4 
# 956105.9  955899.1  953602.4  951539.7  941166.2  945500.8  934269.2  940320.3 

