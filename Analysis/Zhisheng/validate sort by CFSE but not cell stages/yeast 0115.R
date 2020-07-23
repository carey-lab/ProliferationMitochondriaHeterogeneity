# Jan 15 2019

yeast <- read.table("../../G1_M and Fast_Slow/FitFlow__Fast_Slow.tab",header = F,sep = "\t",quote = "")  # the gene IMP' is weird
#remove -
yeast <- yeast[-1,]


imp2 <- yeast[grep("'",yeast$V2),]

yeast[grep("'",yeast$V2),2] <- "wrong"

head(yeast)

nrow(yeast) #6720
# tools::showNonASCII(yeast$V2)
# max(nchar(yeast$V2))
# get corresponding human genes 
require("biomaRt")


mart1 = useMart("ensembl", dataset="hsapiens_gene_ensembl")
mart2 = useMart("ensembl", dataset="scerevisiae_gene_ensembl") 

# human / mouse
human_gene <- getLDS(attributes=("ensembl_gene_id"),
                     filters="ensembl_gene_id", 
                     values=yeast$V1, mart=mart2,
                     attributesL=c("external_gene_name"), martL=mart1)
head(human_gene)

write.table(human_gene,file="./yeast2human id.tsv",col.names = T,row.names = F,sep = "\t",quote = F)

nrow(human_gene)

# ORINGINAL METHOD 
yeast$human_symbol <- human_gene$Gene.name[match(yeast$V1,human_gene$Gene.stable.ID)]
nrow(yeast)

yeast <- na.omit(yeast) #2766

yeast <- yeast[,c(6,5,4)]

# change method at Jan 20
yeast_change <- yeast[match(human_gene$Gene.stable.ID,yeast$V1),]
head(yeast_change)

yeast_change$human_symbol <- human_gene$Gene.name

yeast_change <- yeast_change[,c(6,1,5,4)]


colnames(yeast_change) <- c("NAME","Description","slow","fast")

head(yeast_change)

nrow(yeast_change)

write.table(yeast_change,file="../../G1_M and Fast_Slow/yeast_0120.gct",col.names = T,row.names = F,sep = "\t",quote = F)

# human_gene[human_gene$Gene.stable.ID=="YNL110C",]
# yeast[yeast$human_symbol["TSL1"]]
# # try sgd
# library(org.Sc.sgd.db)
# columns(org.Sc.sgd.db)
# symb <- mapIds(org.Sc.sgd.db, keys=as.character(yeast$V2), keytype="GENENAME", column="ENSEMBL")


# transfer to gct


############################# use ethe table ############################

yeast <- read.table("./G1_M and Fast_Slow/FitFlow__Fast_Slow.tab",header = F,sep = "\t",quote = "")  # the gene IMP' is weird
#remove -
# yeast <- yeast[-1,]
# 
# 
# imp2 <- yeast[grep("'",yeast$V2),]
# 
# yeast[grep("'",yeast$V2),2] <- "wrong"

head(yeast)

human2yeast <- read.table("./G1_M and Fast_Slow/Human2Yeast.tab",header = T,sep ="\t")

head(human2yeast)


human2yeast <- human2yeast[,c(1,3)]
nrow(human2yeast)
na.omit(human2yeast)

yeast$human_gene_id <- human2yeast$GeneStableID[match(yeast$V1,human2yeast$SaccharomycesCerevisiaeGeneStableID)]

length(na.omit(yeast$human_gene_id)) #2766

################## they are exactly the same ########################

########################


# Jan 16 2019 
# transfer yeast data to gct file

yeast$Description <- rep("na",nrow(yeast))
yeast <- yeast[,c(1,4,2,3)]
head(yeast)
colnames(yeast) <- c("NAME","Description","slow","fast")


nrow(yeast)

class(yeast$NAME)

write.table(yeast,file="./G1_M and Fast_Slow/yeast_0116.gct",col.names = T,row.names = F,sep = "\t",quote = F)
