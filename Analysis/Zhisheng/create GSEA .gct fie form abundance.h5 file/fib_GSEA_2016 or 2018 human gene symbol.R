#do GSEA with 2016 data and 2018 data separately
#use human symbol level data






#load fib
fib <- read.table("./GSEA/FIB_GSEA_1231.gct",header = T)
head(fib)

#save 2016 data
fib1 <- fib[,c(1,2,3,4,7,8)]
head(fib1)
write.table(fib1,file="./GSEA/fib2016_symbol.gct",row.names = F,col.names = T,quote = F, sep="\t")
#save 2018 data
fib2 <- fib[,c(1,2,5,6,9,10)] 
head(fib2)
write.table(fib2,file="./GSEA/fib2018_symbol.gct",row.names = F,col.names = T,quote = F, sep="\t")

#load esc
esc <- read.table("./GSEA/ESC_GSEA_0104.gct",header = T)
head(esc)

#save 2016 data
esc1 <- esc[,c(1,2,3,4,7,8)]
head(esc1)
write.table(esc1,file="./GSEA/esc2016_symbol.gct",row.names = F,col.names = T,quote = F, sep="\t")
#save 2018 data
esc2 <- esc[,c(1,2,5,6,9,10)]
head(esc2)
write.table(esc2,file="./GSEA/esc2018_symbol.gct",row.names = F,col.names = T,quote = F, sep="\t")






