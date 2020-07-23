#　Jan 16 2019

# yeast sliding window

yeast_e <- read.table("./G1_M and Fast_Slow/FitFlow__Fast_Slow.tab",header = F,sep = "\t", quote="")
head(yeast_e)
# yeast_e <- read.table("./G1_M and Fast_Slow/yeast_0116.gct",header = T,sep = "\t")
# head(yeast_e)

# use the peak order and slide window for yeast
# order by peak time
yeast_peak <- read.table("./periodical genes/cerevisiae_periodic.tsv",header = T,sep = "\t")
head(yeast_peak)
yeast_peak <- yeast_peak[order(yeast_peak$peaktime,decreasing = F),]


yeast_sort <- yeast_e[match(yeast_peak$gene,yeast_e$V1),]

yeast_sort$peaktime <- yeast_peak$peaktime[match(yeast_sort$V1,yeast_peak$gene)]

head(yeast_sort)

yeast_sort[yeast_sort$V2=="TSL1",]

colnames(yeast_sort) <- c("V1","V2","V3", "fast","slow", "peaktime")

yeast_sort$log2f2s <- log2(yeast_sort$fast/yeast_sort$slow)

yeast_sort$median_window_10 <- rep(0,nrow(yeast_sort))



nrow(yeast_sort)

yeast_sort <- na.omit(yeast_sort)

yeast_sort <- yeast_sort[is.finite(yeast_sort$log2f2s),]

nrow(yeast_sort) #5964

barplot(yeast_sort$log2f2s)

head(yeast_sort)
# 10 window size
temp <- c(yeast_sort$log2f2s,yeast_sort$log2f2s)
i=1
for (i in 1:5964 ) { yeast_sort$median_window_10[i]= median(temp[i:i+9])}

barplot(yeast_sort$median_window_10,ylim=c(-6,6),xlab = ("median_window_10"),ylab = ("log2(fast/slow)")
)
# 100 window size
i=1
for (i in 1:5964) { yeast_sort$median_window_100[i]= median(temp[i:i+99])}
barplot(yeast_sort$median_window_100,ylim=c(-6,6),xlab = ("median_window_100"),ylab = ("log2(fast/slow)")
)
# 1000 window size
i=1
for (i in 1:5964) { yeast_sort$median_window_1000[i]= median(temp[i:i+999])}
barplot(yeast_sort$median_window_1000,ylim=c(-6,6),xlab = ("median_window_1000"),ylab = ("log2(fast/slow)")
)



#  get medain of same peaktime genes
# split dataframe is good funtion
out <- split(yeast_sort, f=yeast_sort$peaktime)
yeast_sort_peak <- sapply(out, function(x) median(x$log2f2s))
yeast_sort_peak <- c(yeast_sort_peak,yeast_sort_peak)
barplot(yeast_sort_peak)


# window size 5  get median of same peak time
yeast_sort_peak <- sapply(out, function(x) median(x$log2f2s))
length(yeast_sort_peak)
yeast_sort_peak <- c(yeast_sort_peak,yeast_sort_peak)
yeast_sort_peak[100]
barplot(yeast_sort_peak[1:100])

for (i in 1:100) { yeast_sort_peak[i]= median(yeast_sort_peak[i:i+4])}
barplot(yeast_sort_peak[1:100],ylim=c(-1,1),xlab = ("median_window_5"),ylab = ("log2(fast/slow)")
)
title("yeast: first get median of same peak time \n then plot median of window size 5")
# window size 10  get median of same peak time
yeast_sort_peak <- sapply(out, function(x) median(x$log2f2s))
yeast_sort_peak <- c(yeast_sort_peak,yeast_sort_peak)
for (i in 1:100) { yeast_sort_peak[i]= median(yeast_sort_peak[i:i+9])}
barplot(yeast_sort_peak[1:100],ylim=c(-1,1),xlab = ("median_window_10"),ylab = ("log2(fast/slow)")
)
title("yeast: first get median of same peak time \n then plot median of window size 10")











