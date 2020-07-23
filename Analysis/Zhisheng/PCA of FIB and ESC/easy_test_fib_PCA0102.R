#########################

#Jan 2 2019
#plot PCA of CFSE FIB and/or ESC
library(DESeq2)
library(gplots)
library(ggplot2)
library(ggfortify)
fib <- read.table("./GSEA/FIB_GSEA_1231.gct", header = T)
head(fib)

require(iris)
head(iris)
head(fib)
fib <- fib[,c(3,4,5,6,7,8,9,10)]
fib <- t(fib)
rownames(fib1)
fib1 <- as.data.frame(fib)

fib1$name <- rownames(fib1)
autoplot(prcomp(fib1[,-19383]),data = fib1, colour = 'name')
iris[,-1]
ncol(fib1)





