# use the G1/M data compare with fast/slow
# data GSE75066

GSE <- read.table("./G1_M and Fast_Slow/GSE75066_exonic_counts.tsv",header = T,sep = "\t")
head(GSE)

nrow(GSE)

GSE <- GSE[,c(2,7,8,9,13,14,15,19,20,21)]

GSE[,c(2,3,4,5,6,7,8,9,10)] <- t(t(GSE[,c(2,3,4,5,6,7,8,9,10)])/colSums(GSE[,c(2,3,4,5,6,7,8,9,10)])*10^6)

colSums(GSE[,c(2,3,4,5,6,7,8,9,10)])
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

# get peak time genes
peaktime <- read.table("./G1_M and Fast_Slow/peaktime_symbol.tsv",header = T)
head(peaktime)

GSE <- GSE[match(peaktime$external_gene_name,GSE$human_symbol),]
GSE <- na.omit(GSE)
nrow(GSE) #342
head(GSE)


# sort GSE expression data

GSE$peak_time <- peaktime$peak_time[na.omit(match(GSE$human_symbol,peaktime$external_gene_name))]

head(GSE)

rownames(GSE) <- c(1:nrow(GSE))

# get the log2(mean(fast)/mean(slow))
class(GSE)
GSE$mean_G2 <- rowMeans(GSE[,c(8,9,10)])

GSE$mean_EG1 <- rowMeans(GSE[,c(2,3,4)])

GSE$mean_LG1 <- rowMeans(GSE[,c(5,6,7)])

# remove 0 values
nrow(GSE)
GSE <- GSE[! rowSums(GSE[,c(2,3,4)])==0,]
GSE <- GSE[! rowSums(GSE[,c(5,6,7)])==0,]
GSE <- GSE[! rowSums(GSE[,c(8,9,10)])==0,]


GSE$log2f2s_E <- log2(GSE$mean_EG1/GSE$mean_G2)

GSE$log2f2s_L <- log2(GSE$mean_LG1/GSE$mean_G2)

#GSE$number <- as.numeric(rownames(GSE)) 

GSE$median_window_E <- rep(0,nrow(GSE))

GSE$median_window_L <- rep(0,nrow(GSE))

#df <- df[is.finite(rowSums(df)),]


barplot(GSE$log2f2s_E[is.finite(GSE$log2f2s_E)],ylim = c(-2,2))

barplot(GSE$log2f2s_L[is.finite(GSE$log2f2s_L)],ylim = c(-2,2))

# FOR log2_E-----------------------

# 10 window size
nrow(GSE)
temp <- c(GSE$log2f2s_E,GSE$log2f2s_E)
i=1
for (i in 1:337 ) { GSE$median_window_E[i]= median(temp[i:i+9])}

barplot(GSE$median_window_E,ylim=c(-2,2),xlab = ("median_window_10"),ylab = ("log2(EG1/G2)")
)
# 5 window size
i=1
for (i in 1:337 ) { GSE$median_window_E[i]= median(temp[i:i+4])}
barplot(GSE$median_window_E,ylim=c(-2,2),xlab = ("median_window_5"),ylab = ("log2(EG1/G2)")
)

#  get medain of same peaktime genes
# split dataframe is good funtion
out <- split(GSE, f=GSE$peak_time)
GSE_peak <- sapply(out, function(x) median(x$log2f2s_E))
GSE_peak <- c(GSE_peak,GSE_peak)

barplot(GSE_peak,ylim = c(-2,2))


# window size 5  get median of same peak time
GSE_peak <- sapply(out, function(x) median(x$log2f2s_E))
length(GSE_peak)
GSE_peak <- c(GSE_peak,GSE_peak)

GSE_peak[79]
barplot(GSE_peak[1:79],ylim = c(-1,1))

for (i in 1:79) { GSE_peak[i]= median(GSE_peak[i:i+4])}
barplot(GSE_peak[1:79],ylim=c(-1,1),xlab = ("median_window_5"),ylab = ("log2(EG1/G2)")
)
title("GSE: first get median of same peak time \n then plot median of window size 5")

# window size 10  get median of same peak time
GSE_peak <- sapply(out, function(x) median(x$log2f2s_E))
GSE_peak <- c(GSE_peak,GSE_peak)
for (i in 1:79) { GSE_peak[i]= median(GSE_peak[i:i+9])}
barplot(GSE_peak[1:79],ylim=c(-1,1),xlab = ("median_window_10"),ylab = ("log2(EG1/G2)")
)
title("GSE: first get median of same peak time \n then plot median of window size 10")

# FOR log2_L-----------------------

# 10 window size
nrow(GSE)
temp <- c(GSE$log2f2s_L,GSE$log2f2s_L)
i=1
for (i in 1:337 ) { GSE$median_window_E[i]= median(temp[i:i+9])}

barplot(GSE$median_window_E,ylim=c(-2,2),xlab = ("median_window_10"),ylab = ("log2(LG1/G2)")
)
# 5 window size
i=1
for (i in 1:337 ) { GSE$median_window_E[i]= median(temp[i:i+4])}
barplot(GSE$median_window_E,ylim=c(-2,2),xlab = ("median_window_5"),ylab = ("log2(LG1/G2)")
)

#  get medain of same peaktime genes
# split dataframe is good funtion
out <- split(GSE, f=GSE$peak_time)
GSE_peak <- sapply(out, function(x) median(x$log2f2s_L))
GSE_peak <- c(GSE_peak,GSE_peak)

barplot(GSE_peak,ylim = c(-2,2))


# window size 5  get median of same peak time
GSE_peak <- sapply(out, function(x) median(x$log2f2s_L))
length(GSE_peak)
GSE_peak <- c(GSE_peak,GSE_peak)

GSE_peak[79]
barplot(GSE_peak[1:79],ylim = c(-1,1))

for (i in 1:79) { GSE_peak[i]= median(GSE_peak[i:i+4])}
barplot(GSE_peak[1:79],ylim=c(-1,1),xlab = ("median_window_5"),ylab = ("log2(LG1/G2)")
)
title("GSE: first get median of same peak time \n then plot median of window size 5")

# window size 10  get median of same peak time
GSE_peak <- sapply(out, function(x) median(x$log2f2s_L))
GSE_peak <- c(GSE_peak,GSE_peak)
for (i in 1:79) { GSE_peak[i]= median(GSE_peak[i:i+9])}
barplot(GSE_peak[1:79],ylim=c(-1,1),xlab = ("median_window_10"),ylab = ("log2(LG1/G2)")
)
title("GSE: first get median of same peak time \n then plot median of window size 10")




