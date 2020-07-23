# Jan 10 2019
# get the corresponding periodical genes from fib and esc rna-seq data

human_pd <- read.table("./periodical genes/human_periodic.tsv",header = F)
head(human_pd)
tail(human_pd)

nrow(human_pd)

require("biomaRt")
listMarts()
mart <- useMart("ENSEMBL_MART_ENSEMBL")
listmart<-listDatasets(mart)
View(listmart)

mart <- useDataset("hsapiens_gene_ensembl", mart)


annote <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id","ensembl_peptide_id","external_gene_name"),
  filter="ensembl_peptide_id",
  values=human_pd$V2,
  uniqueRows=T)

write.table(annote,file="./human_periodical_gene_symbol.tsv",col.names = T,row.names = F, sep = "\t",quote = F)

# library(org.Hs.eg.db)
# 
# symb <- mapIds(org.Hs.eg.db, keys=as.character(human_pd$V2) , keytype="ENSEMBLPROT", column="SYMBOL")
# head(symb)



# human_pd$V2[! human_pd$V2 %in% annote$ensembl_peptide_id]
# [1] ENSP00000217429 ENSP00000280193 ENSP00000318089 ENSP00000384703 ENSP00000358153
# [6] ENSP00000358162 ENSP00000361882 ENSP00000366581 ENSP00000308344 ENSP00000324074
# [11] ENSP00000355810 ENSP00000414687 ENSP00000305056 ENSP00000345892 ENSP00000357380

# YOU CAN NOT use biomart to get these gene symbol, so I add them manually.
# now read it into R

human_pd_symbol <- read.table("./human_periodical_gene_symbol.tsv",header = T)
head(human_pd_symbol)

# use symbol to get expression data
fib <- read.table("./GSEA/fib2016+2018_rb_negative.gct",header = T)
head(fib)

fib_human_pd_symbol <- fib[match(human_pd_symbol$external_gene_name,fib$NAME),]
fib_human_pd_symbol <- na.omit(fib_human_pd_symbol)
nrow(fib_human_pd_symbol) #355

head(fib_human_pd_symbol)
rownames(fib_human_pd_symbol) <- c(1:nrow(fib_human_pd_symbol))

# get the log2(mean(fast)/mean(slow))
class(fib_human_pd_symbol)
fib_human_pd_symbol$mean_slow <- rowMeans(fib_human_pd_symbol[,c(3,4,5,6)])

fib_human_pd_symbol$mean_fast <- rowMeans(fib_human_pd_symbol[,c(7,8,9,10)])

fib_human_pd_symbol$log2f2s <- log2(fib_human_pd_symbol$mean_fast/fib_human_pd_symbol$mean_slow)

fib_human_pd_symbol$number <- as.numeric(rownames(fib_human_pd_symbol)) 


# get the mean(log2(fast/slow))--------




# plot 

plot(fib_human_pd_symbol$log2f2s,ylim=c(-2,2),xlab = ("rank"),ylab = ("log2(fast/slow)"),
     pch=20,
     col="black",
     cex=1
     )

abline(lm(fib_human_pd_symbol$log2f2s ~fib_human_pd_symbol$number ))
text(x=150,y=1.8,"Intercept:-0.1885643, slop:-0.0004583")
title("Fibroblasts periodical genes")


# use symbol to get expression data
esc <- read.table("./GSEA/esc2016+2018_rb_negative.gct",header = T)
head(esc)

esc_human_pd_symbol <- esc[match(human_pd_symbol$external_gene_name,esc$NAME),]
esc_human_pd_symbol <- na.omit(esc_human_pd_symbol)
nrow(esc_human_pd_symbol) #357

head(esc_human_pd_symbol)
rownames(esc_human_pd_symbol) <- c(1:nrow(esc_human_pd_symbol))

# get the log2(mean(fast)/mean(slow))
class(esc_human_pd_symbol)
esc_human_pd_symbol$mean_slow <- rowMeans(esc_human_pd_symbol[,c(3,4,5,6)])

esc_human_pd_symbol$mean_fast <- rowMeans(esc_human_pd_symbol[,c(7,8,9,10)])

esc_human_pd_symbol$log2f2s <- log2(esc_human_pd_symbol$mean_fast/esc_human_pd_symbol$mean_slow)

esc_human_pd_symbol$number <- as.numeric(rownames(esc_human_pd_symbol)) 


plot(esc_human_pd_symbol$log2f2s,ylim=c(-2,2),xlab = ("rank"),ylab = ("log2(fast/slow)"),
     pch=20,
     col="black",
     cex=1
)

abline(lm(esc_human_pd_symbol$log2f2s ~esc_human_pd_symbol$number ))
text(x=150,y=1.8,"Intercept:--5.761e-02, slop:-3.733e-05")
title("ESCs periodical genes")



# use the peak order------------------ version 2
# YOU CAN NOT use biomart to get these gene symbol, so I add them manually.
# now read it into R

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



# use symbol to get expression data--------------------------
fib <- read.table("./GSEA/fib2016+2018_rb_negative.gct",header = T)
head(fib)

fib_human_pd_symbol3_sort <- fib[match(human_pd_symbol3_sort$external_gene_name,fib$NAME),]
nrow(fib_human_pd_symbol3_sort)

fib_human_pd_symbol3_sort <- na.omit(fib_human_pd_symbol3_sort)

fib_human_pd_symbol3_sort$peak_time <- human_pd_symbol3_sort$peak_time[na.omit(match(human_pd_symbol3_sort$external_gene_name,fib_human_pd_symbol3_sort$NAME))]

fib_human_pd_symbol3_sort$rank <- human_pd_symbol3_sort$rank[na.omit(match(human_pd_symbol3_sort$external_gene_name,fib_human_pd_symbol3_sort$NAME))]

nrow(fib_human_pd_symbol3_sort) #355

head(fib_human_pd_symbol3_sort)
rownames(fib_human_pd_symbol3_sort) <- c(1:nrow(fib_human_pd_symbol3_sort))

# get the log2(mean(fast)/mean(slow))
class(fib_human_pd_symbol3_sort)
fib_human_pd_symbol3_sort$mean_slow <- rowMeans(fib_human_pd_symbol3_sort[,c(3,4,5,6)])

fib_human_pd_symbol3_sort$mean_fast <- rowMeans(fib_human_pd_symbol3_sort[,c(7,8,9,10)])

fib_human_pd_symbol3_sort$log2f2s <- log2(fib_human_pd_symbol3_sort$mean_fast/fib_human_pd_symbol3_sort$mean_slow)

fib_human_pd_symbol3_sort$number <- as.numeric(rownames(fib_human_pd_symbol3_sort)) 


# get the mean(log2(fast/slow))--------




# plot 
summary(lm(fib_human_pd_symbol3_sort$log2f2s ~fib_human_pd_symbol3_sort$number ))

plot(fib_human_pd_symbol3_sort$peak_time,fib_human_pd_symbol3_sort$log2f2s,ylim=c(-2,2),xlab = ("peak_time"),ylab = ("log2(fast/slow)"),
     pch=20,
     col="black",
     cex=1
)

abline(lm(fib_human_pd_symbol3_sort$log2f2s ~fib_human_pd_symbol3_sort$number ))
text(x=40,y=1.8,"Intercept:-0.1745682, slop:-0.0005369\n Adjusted R-squared:  0.01275,p-value: 0.0188")
title("Fibroblasts periodical genes")


# use symbol to get expression data
esc <- read.table("./GSEA/esc2016+2018_rb_negative.gct",header = T)
head(esc)

esc_human_pd_symbol3_sort <- esc[match(human_pd_symbol3_sort$external_gene_name,esc$NAME),]
nrow(esc_human_pd_symbol3_sort)

esc_human_pd_symbol3_sort <- na.omit(esc_human_pd_symbol3_sort)

esc_human_pd_symbol3_sort$peak_time <- human_pd_symbol3_sort$peak_time[na.omit(match(human_pd_symbol3_sort$external_gene_name,esc_human_pd_symbol3_sort$NAME))]

esc_human_pd_symbol3_sort$rank <- human_pd_symbol3_sort$rank[na.omit(match(human_pd_symbol3_sort$external_gene_name,esc_human_pd_symbol3_sort$NAME))]

nrow(esc_human_pd_symbol3_sort) #357

head(esc_human_pd_symbol3_sort)
rownames(esc_human_pd_symbol3_sort) <- c(1:nrow(esc_human_pd_symbol3_sort))

# get the log2(mean(fast)/mean(slow))
class(esc_human_pd_symbol3_sort)
esc_human_pd_symbol3_sort$mean_slow <- rowMeans(esc_human_pd_symbol3_sort[,c(3,4,5,6)])

esc_human_pd_symbol3_sort$mean_fast <- rowMeans(esc_human_pd_symbol3_sort[,c(7,8,9,10)])

esc_human_pd_symbol3_sort$log2f2s <- log2(esc_human_pd_symbol3_sort$mean_fast/esc_human_pd_symbol3_sort$mean_slow)

esc_human_pd_symbol3_sort$number <- as.numeric(rownames(esc_human_pd_symbol3_sort)) 


# get the mean(log2(fast/slow))--------




# plot 

plot(esc_human_pd_symbol3_sort$peak_time,esc_human_pd_symbol3_sort$log2f2s,ylim=c(-2,2),xlab = ("peak_time"),ylab = ("log2(fast/slow)"),
     pch=20,
     col="black",
     cex=1
)

abline(lm(esc_human_pd_symbol3_sort$log2f2s ~esc_human_pd_symbol3_sort$number ))
text(x=40,y=1.8,"Intercept:0.0714317, slop:-0.0006836\n Adjusted R-squared:  0.06819,p-value: 3.355e-07")
title("ESCs periodical genes")

# use the peak order and slide window------------------ version 3
# YOU CAN NOT use biomart to get these gene symbol, so I add them manually.
# now read it into R

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



# use symbol to get expression data--------------------------
fib <- read.table("./GSEA/fib2016+2018_rb_negative.gct",header = T)
head(fib)

fib_human_pd_symbol3_sort <- fib[na.omit(match(human_pd_symbol3_sort$external_gene_name,fib$NAME)),]
nrow(fib_human_pd_symbol3_sort)



fib_human_pd_symbol3_sort$peak_time <- human_pd_symbol3_sort$peak_time[na.omit(match(fib_human_pd_symbol3_sort$NAME,human_pd_symbol3_sort$external_gene_name))]

head(fib_human_pd_symbol3_sort)

fib_human_pd_symbol3_sort$rank <- human_pd_symbol3_sort$rank[na.omit(match(fib_human_pd_symbol3_sort$NAME,human_pd_symbol3_sort$external_gene_name))]

nrow(fib_human_pd_symbol3_sort) #355

head(fib_human_pd_symbol3_sort)
rownames(fib_human_pd_symbol3_sort) <- c(1:nrow(fib_human_pd_symbol3_sort))

# get the log2(mean(fast)/mean(slow))
class(fib_human_pd_symbol3_sort)
fib_human_pd_symbol3_sort$mean_slow <- rowMeans(fib_human_pd_symbol3_sort[,c(3,4,5,6)])

fib_human_pd_symbol3_sort$mean_fast <- rowMeans(fib_human_pd_symbol3_sort[,c(7,8,9,10)])

fib_human_pd_symbol3_sort$log2f2s <- log2(fib_human_pd_symbol3_sort$mean_fast/fib_human_pd_symbol3_sort$mean_slow)

fib_human_pd_symbol3_sort$number <- as.numeric(rownames(fib_human_pd_symbol3_sort)) 

fib_human_pd_symbol3_sort$median_window_10 <- rep(0,nrow(fib_human_pd_symbol3_sort))

barplot(fib_human_pd_symbol3_sort$log2f2s)

# 10 window size
temp <- c(fib_human_pd_symbol3_sort$log2f2s,fib_human_pd_symbol3_sort$log2f2s)
i=1
for (i in 1:355 ) { fib_human_pd_symbol3_sort$median_window_10[i]= median(temp[i:i+9])}

barplot(fib_human_pd_symbol3_sort$median_window_10,ylim=c(-2,2),xlab = ("median_window_10"),ylab = ("log2(fast/slow)")
)
# 5 window size
i=1
for (i in 1:355) { fib_human_pd_symbol3_sort$median_window_10[i]= median(temp[i:i+4])}
barplot(fib_human_pd_symbol3_sort$median_window_10,ylim=c(-2,2),xlab = ("median_window_5"),ylab = ("log2(fast/slow)")
)

#  get medain of same peaktime genes
# split dataframe is good funtion
out <- split(fib_human_pd_symbol3_sort, f=fib_human_pd_symbol3_sort$peak_time)
fib_human_pd_symbol3_sort_peak <- sapply(out, function(x) median(x$log2f2s))
fib_human_pd_symbol3_sort_peak <- c(fib_human_pd_symbol3_sort_peak,fib_human_pd_symbol3_sort_peak)
barplot(fib_human_pd_symbol3_sort_peak)


# window size 5  get median of same peak time
fib_human_pd_symbol3_sort_peak <- sapply(out, function(x) median(x$log2f2s))
length(fib_human_pd_symbol3_sort_peak)
fib_human_pd_symbol3_sort_peak <- c(fib_human_pd_symbol3_sort_peak,fib_human_pd_symbol3_sort_peak)

fib_human_pd_symbol3_sort_peak[79]
barplot(fib_human_pd_symbol3_sort_peak[1:79])

for (i in 1:79) { fib_human_pd_symbol3_sort_peak[i]= median(fib_human_pd_symbol3_sort_peak[i:i+4])}
barplot(fib_human_pd_symbol3_sort_peak[1:79],ylim=c(-2,2),xlab = ("median_window_5"),ylab = ("log2(fast/slow)")
)
title("fibroblast: first get median of same peak time \n then plot median of window size 5")
# window size 10  get median of same peak time
fib_human_pd_symbol3_sort_peak <- sapply(out, function(x) median(x$log2f2s))
fib_human_pd_symbol3_sort_peak <- c(fib_human_pd_symbol3_sort_peak,fib_human_pd_symbol3_sort_peak)
for (i in 1:79) { fib_human_pd_symbol3_sort_peak[i]= median(fib_human_pd_symbol3_sort_peak[i:i+9])}
barplot(fib_human_pd_symbol3_sort_peak[1:79],ylim=c(-2,2),xlab = ("median_window_10"),ylab = ("log2(fast/slow)")
)
title("fibroblast: first get median of same peak time \n then plot median of window size 10")



# get the mean(log2(fast/slow))--------




# plot 
summary(lm(fib_human_pd_symbol3_sort$median_window_10 ~ fib_human_pd_symbol3_sort$number ))

barplot(fib_human_pd_symbol3_sort$median_window_10)

plot(fib_human_pd_symbol3_sort$peak_time,fib_human_pd_symbol3_sort$median_window_10,ylim=c(-2,2),xlab = ("median_window_10"),ylab = ("log2(fast/slow)"),
     pch=20,
     col="black",
     cex=1
)

abline(lm(fib_human_pd_symbol3_sort$log2f2s ~fib_human_pd_symbol3_sort$number ))
text(x=40,y=1.8,"Intercept:-0.1745682, slop:-0.0005369\n Adjusted R-squared:  0.01275,p-value: 0.0188")
title("Fibroblasts periodical genes")


# use symbol to get expression data--------------
esc <- read.table("./GSEA/esc2016+2018_rb_negative.gct",header = T)
head(esc)

esc_human_pd_symbol3_sort <- esc[na.omit(match(human_pd_symbol3_sort$external_gene_name,esc$NAME)),]
nrow(esc_human_pd_symbol3_sort)



esc_human_pd_symbol3_sort$peak_time <- human_pd_symbol3_sort$peak_time[na.omit(match(esc_human_pd_symbol3_sort$NAME,human_pd_symbol3_sort$external_gene_name))]

head(esc_human_pd_symbol3_sort)

esc_human_pd_symbol3_sort$rank <- human_pd_symbol3_sort$rank[na.omit(match(esc_human_pd_symbol3_sort$NAME,human_pd_symbol3_sort$external_gene_name))]

nrow(esc_human_pd_symbol3_sort) #357

head(esc_human_pd_symbol3_sort)
rownames(esc_human_pd_symbol3_sort) <- c(1:nrow(esc_human_pd_symbol3_sort))

# get the log2(mean(fast)/mean(slow))
class(esc_human_pd_symbol3_sort)
esc_human_pd_symbol3_sort$mean_slow <- rowMeans(esc_human_pd_symbol3_sort[,c(3,4,5,6)])

esc_human_pd_symbol3_sort$mean_fast <- rowMeans(esc_human_pd_symbol3_sort[,c(7,8,9,10)])

esc_human_pd_symbol3_sort$log2f2s <- log2(esc_human_pd_symbol3_sort$mean_fast/esc_human_pd_symbol3_sort$mean_slow)

esc_human_pd_symbol3_sort$number <- as.numeric(rownames(esc_human_pd_symbol3_sort)) 

esc_human_pd_symbol3_sort$median_window_10 <- rep(0,nrow(esc_human_pd_symbol3_sort))



# 10 window size
temp <- c(esc_human_pd_symbol3_sort$log2f2s,esc_human_pd_symbol3_sort$log2f2s)
i=1
for (i in 1:357 ) { esc_human_pd_symbol3_sort$median_window_10[i]= median(temp[i:i+9])}

barplot(esc_human_pd_symbol3_sort$median_window_10,ylim=c(-2,2),xlab = ("median_window_10"),ylab = ("log2(fast/slow)")
)
# 5 window size
i=1
for (i in 1:357) { esc_human_pd_symbol3_sort$median_window_10[i]= median(temp[i:i+4])}
barplot(esc_human_pd_symbol3_sort$median_window_10,ylim=c(-2,2),xlab = ("median_window_5"),ylab = ("log2(fast/slow)")
)

barplot(esc_human_pd_symbol3_sort$log2f2s)
#  get medain of same peaktime genes
# split dataframe is good funtion
out <- split(esc_human_pd_symbol3_sort, f=esc_human_pd_symbol3_sort$peak_time)
esc_human_pd_symbol3_sort_peak <- sapply(out, function(x) median(x$log2f2s))
esc_human_pd_symbol3_sort_peak <- c(esc_human_pd_symbol3_sort_peak,esc_human_pd_symbol3_sort_peak)
barplot(esc_human_pd_symbol3_sort_peak[1:79])


# window size 5  get median of same peak time
esc_human_pd_symbol3_sort_peak <- sapply(out, function(x) median(x$log2f2s))
length(esc_human_pd_symbol3_sort_peak)
esc_human_pd_symbol3_sort_peak <- c(esc_human_pd_symbol3_sort_peak,esc_human_pd_symbol3_sort_peak)

esc_human_pd_symbol3_sort_peak[79]


for (i in 1:79) { esc_human_pd_symbol3_sort_peak[i]= median(esc_human_pd_symbol3_sort_peak[i:i+4])}
barplot(esc_human_pd_symbol3_sort_peak[1:79],ylim=c(-1,1),xlab = ("median_window_5"),ylab = ("log2(fast/slow)")
)
title("EScs: first get median of same peak time \n then plot median of window size 5")
# window size 10  get median of same peak time
esc_human_pd_symbol3_sort_peak <- sapply(out, function(x) median(x$log2f2s))
esc_human_pd_symbol3_sort_peak <- c(esc_human_pd_symbol3_sort_peak,esc_human_pd_symbol3_sort_peak)
for (i in 1:79) { esc_human_pd_symbol3_sort_peak[i]= median(esc_human_pd_symbol3_sort_peak[i:i+9])}
barplot(esc_human_pd_symbol3_sort_peak[1:79],ylim=c(-1,1),xlab = ("median_window_10"),ylab = ("log2(fast/slow)")
)
title("EScs: first get median of same peak time \n then plot median of window size 10")
