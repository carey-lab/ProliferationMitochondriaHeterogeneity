# use the G1/M data compare with fast/slow
# data sox2

sox2_EG1_G2 <- sox2[,c(13,2,3,4,8,9,10)]
  


head(sox2_EG1_G2)
head(GES_EG1_G2)


# transfer gct format, TPM  sox2_EG1_G2
sox2_EG1_G2$Description <- rep("na",nrow(sox2_EG1_G2))
sox2_EG1_G2 <- sox2_EG1_G2[,c(1,8,2,3,4,5,6,7)]

colnames(sox2_EG1_G2) #19238
 write.table(sox2_EG1_G2,file="./G1_M and Fast_Slow/sox2_EG1_G2.gct",col.names = T,row.names = F,quote = F,sep = "\t")

 

 
 
# load GSE
 GSE <- read.table("./G1_M and Fast_Slow/GSE75066_exonic_counts.tsv",header = T,sep = "\t")
 head(GSE)
 
 

 
 # transfer mouse gene symbol to human gene symbol
 require("biomaRt")
 mart1 = useMart("ensembl", dataset="hsapiens_gene_ensembl")
 mart2 = useMart("ensembl", dataset="mmusculus_gene_ensembl") 
 
 # human / mouse
 human_gene <- getLDS(attributes=c("external_gene_name"),
                      filters="external_gene_name", values=GSE$Symbol, mart=mart2,
                      attributesL=c("external_gene_name"), martL=mart1)
 head(human_gene)
 
 # add human symbol
 GSE$human_symbol <- human_gene$Gene.name.1[match(GSE$Symbol,human_gene$Gene.name)]
 nrow(GSE)

 head(GSE)
 GSE <- GSE[,c(25,7,8,9,19,20,21)]
 
 
# for (i in 2:10) {GSE[,i]= GSE[,i]}
#  # GSE[,c(2,3,4,5,6,7,8,9,10)] <- 
 
 # divide by row
 GSE[,c(2,3,4,5,6,7)] <- t(t(GSE[,c(2,3,4,5,6,7)])/colSums(GSE[,c(2,3,4,5,6,7)])*10^6)
 colSums( GSE[,c(2,3,4,5,6,7)])
# GES_EG1_G2 
 GES_EG1_G2 <- GSE[,c(1,2,3,4,5,6,7)]
 
 GES_EG1_G2$Description <- rep("na",nrow(GES_EG1_G2))
 GES_EG1_G2 <- GES_EG1_G2[,c(1,8,2,3,4,5,6,7)]
 
 colnames(GES_EG1_G2) <- c("NAME","Description","EG1_1","EG1_2","EG1_3","G2_1","G2_2","G2_3")
 GES_EG1_G2 <- na.omit(GES_EG1_G2)
 nrow(GES_EG1_G2) #17250
 head(GES_EG1_G2)
 
 write.table(GES_EG1_G2,file="./G1_M and Fast_Slow/GES_EG1_G2.gct",col.names = T,row.names = F,quote = F,sep = "\t")
 
 
 
 
 
 
 
 
 