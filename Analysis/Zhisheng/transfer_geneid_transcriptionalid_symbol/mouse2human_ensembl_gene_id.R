
#使用biomart提取与小鼠基因id对应的人类基因id和基因名

require("biomaRt")
mart1 = useMart("ensembl", dataset="hsapiens_gene_ensembl")
mart2 = useMart("ensembl", dataset="mmusculus_gene_ensembl") 

mouse_gene_ids <- annotLookup2$ensembl_gene_id

# human / mouse
human_gene <- getLDS(attributes=c("ensembl_gene_id"),
       filters="ensembl_gene_id", values=mouse_gene_ids, mart=mart2,
       attributesL=c("ensembl_gene_id","external_gene_name"), martL=mart1)


colnames(human_gene) <- c("mouse_ensembl_gene_id","human_ensembl_gene_id","human_external_gene_name")


write.table(human_gene,file = "mouse2human_ensembl_gene_id.table",row.names = F,quote = F)


head(human_gene)
head(annotLookup2)
