# use the G1/M data compare with fast/slow
# data sox2

sox2 <- read.table("./G1_M and Fast_Slow/sox2_tpm.tsv",header = T)
head(sox2)
sox2 <- sox2[,c(1,5,6,7,11,12,13,17,18,19)]
sox2$mNAME <- toupper(gsub(".*_","",   sox2$Gene))
sox2$geneid <- gsub("_.*","",   sox2$Gene)

# transfer mouse gene id to human gene symbol
require("biomaRt")
mart1 = useMart("ensembl", dataset="hsapiens_gene_ensembl")
mart2 = useMart("ensembl", dataset="mmusculus_gene_ensembl") 

# human / mouse
human_gene <- getLDS(attributes=c("ensembl_gene_id","external_gene_name"),
                     filters="ensembl_gene_id", values=sox2$geneid, mart=mart2,
                     attributesL=c("ensembl_gene_id","external_gene_name"), martL=mart1)

exist_human_symbol <- human_gene[na.omit(match(sox2$geneid,human_gene$Gene.stable.ID)),]

sox2 <- sox2[match(exist_human_symbol$Gene.stable.ID,sox2$geneid),]
head(sox2)
nrow(sox2)
head(exist_human_symbol)
sox2$NAME <- exist_human_symbol$Gene.name.1


length(na.omit(match(sox2$geneid,human_gene$Gene.stable.ID)))
length(match(sox2$geneid,human_gene$Gene.stable.ID))
nrow(human_gene)
nrow(sox2)

# plot cell cycle genes log2(G1/M)
human_pd_symbol <- read.table("./human_periodical_gene_symbol.tsv",header = T)
head(human_pd_symbol)
nrow(human_pd_symbol)

# want to get the peak time
human_pd <- read.table("./periodical genes/human_periodic.tsv",header = F)
head(human_pd)

human_pd_symbol2 <- human_pd[match(human_pd_symbol$ensembl_peptide_id,human_pd$V2),]
human_pd_symbol2 <- human_pd_symbol2[,c(2,3,4)]
head(human_pd_symbol2)

human_pd_symbol3 <- cbind(human_pd_symbol,human_pd_symbol2)
human_pd_symbol3 <- human_pd_symbol3[,c(2,3,5,6)]
head(human_pd_symbol3)
colnames(human_pd_symbol3) <- c("ensembl_peptide_id","external_gene_name","rank","peak_time")
# order by peak time
human_pd_symbol3_sort <- human_pd_symbol3[order(human_pd_symbol3$peak_time),]
head(human_pd_symbol3_sort)

write.table(human_pd_symbol3_sort,file="./G1_M and Fast_Slow/peaktime_symbol.tsv",col.names = T,row.names = F,sep = "\t",quote = F)

# sort sox2 expression data

sox2_human_pd_symbol3_sort <- sox2[na.omit(match(human_pd_symbol3_sort$external_gene_name,sox2$NAME)),]
nrow(sox2_human_pd_symbol3_sort) #350



sox2_human_pd_symbol3_sort$peak_time <- human_pd_symbol3_sort$peak_time[na.omit(match(sox2_human_pd_symbol3_sort$NAME,human_pd_symbol3_sort$external_gene_name))]

head(sox2_human_pd_symbol3_sort)

sox2_human_pd_symbol3_sort$rank <- human_pd_symbol3_sort$rank[na.omit(match(sox2_human_pd_symbol3_sort$NAME,human_pd_symbol3_sort$external_gene_name))]

nrow(sox2_human_pd_symbol3_sort) #350

head(sox2_human_pd_symbol3_sort)
rownames(sox2_human_pd_symbol3_sort) <- c(1:nrow(sox2_human_pd_symbol3_sort))

# get the log2(mean(fast)/mean(slow))
class(sox2_human_pd_symbol3_sort)
sox2_human_pd_symbol3_sort$mean_G2_m <- rowMeans(sox2_human_pd_symbol3_sort[,c(8,9,10)])

sox2_human_pd_symbol3_sort$mean_EG1_m <- rowMeans(sox2_human_pd_symbol3_sort[,c(2,3,4)])

sox2_human_pd_symbol3_sort$mean_LG1_m <- rowMeans(sox2_human_pd_symbol3_sort[,c(5,6,7)])

# remove 0 values
nrow(sox2_human_pd_symbol3_sort)
sox2_human_pd_symbol3_sort <- sox2_human_pd_symbol3_sort[! rowSums(sox2_human_pd_symbol3_sort[,c(2,3,4)])==0,]
sox2_human_pd_symbol3_sort <- sox2_human_pd_symbol3_sort[! rowSums(sox2_human_pd_symbol3_sort[,c(5,6,7)])==0,]
sox2_human_pd_symbol3_sort <- sox2_human_pd_symbol3_sort[! rowSums(sox2_human_pd_symbol3_sort[,c(8,9,10)])==0,]


sox2_human_pd_symbol3_sort$log2f2s_E <- log2(sox2_human_pd_symbol3_sort$mean_EG1_m/sox2_human_pd_symbol3_sort$mean_G2_m)

sox2_human_pd_symbol3_sort$log2f2s_L <- log2(sox2_human_pd_symbol3_sort$mean_LG1_m/sox2_human_pd_symbol3_sort$mean_G2_m)

#sox2_human_pd_symbol3_sort$number <- as.numeric(rownames(sox2_human_pd_symbol3_sort)) 

sox2_human_pd_symbol3_sort$median_window_E <- rep(0,nrow(sox2_human_pd_symbol3_sort))

sox2_human_pd_symbol3_sort$median_window_L <- rep(0,nrow(sox2_human_pd_symbol3_sort))

#df <- df[is.finite(rowSums(df)),]


barplot(sox2_human_pd_symbol3_sort$log2f2s_E[is.finite(sox2_human_pd_symbol3_sort$log2f2s_E)],ylim = c(-2,2))

barplot(sox2_human_pd_symbol3_sort$log2f2s_L[is.finite(sox2_human_pd_symbol3_sort$log2f2s_L)],ylim = c(-2,2))

# FOR log2_E-----------------------

# 10 window size
nrow(sox2_human_pd_symbol3_sort)
temp <- c(sox2_human_pd_symbol3_sort$log2f2s_E,sox2_human_pd_symbol3_sort$log2f2s_E)
i=1
for (i in 1:339 ) { sox2_human_pd_symbol3_sort$median_window_E[i]= median(temp[i:i+9])}

barplot(sox2_human_pd_symbol3_sort$median_window_E,ylim=c(-2,2),xlab = ("median_window_10"),ylab = ("log2(EG1/G2)")
)
# 5 window size
i=1
for (i in 1:339 ) { sox2_human_pd_symbol3_sort$median_window_E[i]= median(temp[i:i+4])}
barplot(sox2_human_pd_symbol3_sort$median_window_E,ylim=c(-2,2),xlab = ("median_window_5"),ylab = ("log2(EG1/G2)")
)

#  get medain of same peaktime genes
# split dataframe is good funtion
out <- split(sox2_human_pd_symbol3_sort, f=sox2_human_pd_symbol3_sort$peak_time)
sox2_human_pd_symbol3_sort_peak <- sapply(out, function(x) median(x$log2f2s_E))
sox2_human_pd_symbol3_sort_peak <- c(sox2_human_pd_symbol3_sort_peak,sox2_human_pd_symbol3_sort_peak)

barplot(sox2_human_pd_symbol3_sort_peak,ylim = c(-2,2))


# window size 5  get median of same peak time
sox2_human_pd_symbol3_sort_peak <- sapply(out, function(x) median(x$log2f2s_E))
length(sox2_human_pd_symbol3_sort_peak)
sox2_human_pd_symbol3_sort_peak <- c(sox2_human_pd_symbol3_sort_peak,sox2_human_pd_symbol3_sort_peak)

sox2_human_pd_symbol3_sort_peak[77]
barplot(sox2_human_pd_symbol3_sort_peak[1:77],ylim = c(-2,2))

for (i in 1:77) { sox2_human_pd_symbol3_sort_peak[i]= median(sox2_human_pd_symbol3_sort_peak[i:i+4])}
barplot(sox2_human_pd_symbol3_sort_peak[1:77],ylim=c(-1,1),xlab = ("median_window_5"),ylab = ("log2(EG1/G2)")
)
title("sox2: first get median of same peak time \n then plot median of window size 5")

# window size 10  get median of same peak time
sox2_human_pd_symbol3_sort_peak <- sapply(out, function(x) median(x$log2f2s_E))
sox2_human_pd_symbol3_sort_peak <- c(sox2_human_pd_symbol3_sort_peak,sox2_human_pd_symbol3_sort_peak)
for (i in 1:77) { sox2_human_pd_symbol3_sort_peak[i]= median(sox2_human_pd_symbol3_sort_peak[i:i+9])}
barplot(sox2_human_pd_symbol3_sort_peak[1:77],ylim=c(-1,1),xlab = ("median_window_10"),ylab = ("log2(EG1/G2)")
)
title("sox2: first get median of same peak time \n then plot median of window size 10")



# FOR log2_L-----------------------

# 10 window size
nrow(sox2_human_pd_symbol3_sort)
temp <- c(sox2_human_pd_symbol3_sort$log2f2s_L,sox2_human_pd_symbol3_sort$log2f2s_L)
i=1
for (i in 1:339 ) { sox2_human_pd_symbol3_sort$median_window_E[i]= median(temp[i:i+9])}

barplot(sox2_human_pd_symbol3_sort$median_window_E,ylim=c(-2,2),xlab = ("median_window_10"),ylab = ("log2(LG1/G2)")
)
# 5 window size
i=1
for (i in 1:339 ) { sox2_human_pd_symbol3_sort$median_window_E[i]= median(temp[i:i+4])}
barplot(sox2_human_pd_symbol3_sort$median_window_E,ylim=c(-2,2),xlab = ("median_window_5"),ylab = ("log2(LG1/G2)")
)

#  get medain of same peaktime genes
# split dataframe is good funtion
out <- split(sox2_human_pd_symbol3_sort, f=sox2_human_pd_symbol3_sort$peak_time)
sox2_human_pd_symbol3_sort_peak <- sapply(out, function(x) median(x$log2f2s_L))
sox2_human_pd_symbol3_sort_peak <- c(sox2_human_pd_symbol3_sort_peak,sox2_human_pd_symbol3_sort_peak)

barplot(sox2_human_pd_symbol3_sort_peak,ylim = c(-2,2))


# window size 5  get median of same peak time
sox2_human_pd_symbol3_sort_peak <- sapply(out, function(x) median(x$log2f2s_L))
length(sox2_human_pd_symbol3_sort_peak)
sox2_human_pd_symbol3_sort_peak <- c(sox2_human_pd_symbol3_sort_peak,sox2_human_pd_symbol3_sort_peak)

sox2_human_pd_symbol3_sort_peak[77]
barplot(sox2_human_pd_symbol3_sort_peak[1:77],ylim = c(-2,2))

for (i in 1:77) { sox2_human_pd_symbol3_sort_peak[i]= median(sox2_human_pd_symbol3_sort_peak[i:i+4])}
barplot(sox2_human_pd_symbol3_sort_peak[1:77],ylim=c(-1,1),xlab = ("median_window_5"),ylab = ("log2(LG1/G2)")
)
title("sox2: first get median of same peak time \n then plot median of window size 5")

# window size 10  get median of same peak time
sox2_human_pd_symbol3_sort_peak <- sapply(out, function(x) median(x$log2f2s_L))
sox2_human_pd_symbol3_sort_peak <- c(sox2_human_pd_symbol3_sort_peak,sox2_human_pd_symbol3_sort_peak)
for (i in 1:77) { sox2_human_pd_symbol3_sort_peak[i]= median(sox2_human_pd_symbol3_sort_peak[i:i+9])}
barplot(sox2_human_pd_symbol3_sort_peak[1:77],ylim=c(-1,1),xlab = ("median_window_10"),ylab = ("log2(LG1/G2)")
)
title("sox2: first get median of same peak time \n then plot median of window size 10")













