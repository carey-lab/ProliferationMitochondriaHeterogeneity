
# Jan 08 2019 
# get corresponding categories for all gene sets

# fib2016_2018_fast_rb_sn_ng data-----
fib2016_2018_fast_rb_sn_ng <- read.table("./GSEA/result_from_server2/fib2016+2018_rb_ng_sn_gsea_report_for_FAST_1546779207404.tsv", header = T,sep="\t")
fib2016_2018_slow_rb_sn_ng <- read.table("./GSEA/result_from_server2/fib2016+2018_rb_ng_sn_gsea_report_for_SLOW_1546779207404.tsv", header = T,sep="\t")



fib2016_2018_fast_rb_sn_ng <-as.data.frame(fib2016_2018_fast_rb_sn_ng[,c(1,4,6,7,8)])
fib2016_2018_slow_rb_sn_ng <-as.data.frame(fib2016_2018_slow_rb_sn_ng[,c(1,4,6,7,8)])


fib2016_2018_rb_sn_ng <- rbind(fib2016_2018_slow_rb_sn_ng,fib2016_2018_fast_rb_sn_ng)
fib2016_2018_rb_sn_ng <- na.omit(fib2016_2018_rb_sn_ng)

head(fib2016_2018_rb_sn_ng)
nrow(fib2016_2018_rb_sn_ng)

#-----------load esc2016_2018_rb_sn_ng data----

esc2016_2018_fast_rb_sn_ng <- read.table("./GSEA/result_from_server2/esc2016+2018_rb_negative_sn_gsea_report_for_FAST_1546782445406.tsv", header = T,sep="\t")
esc2016_2018_slow_rb_sn_ng <- read.table("./GSEA/result_from_server2/esc2016+2018_rb_negative_sn_gsea_report_for_SLOW_1546782445406.tsv", header = T,sep="\t")


esc2016_2018_fast_rb_sn_ng <-as.data.frame(esc2016_2018_fast_rb_sn_ng[,c(1,4,6,7,8)])
esc2016_2018_slow_rb_sn_ng <-as.data.frame(esc2016_2018_slow_rb_sn_ng[,c(1,4,6,7,8)])

esc2016_2018_rb_sn_ng <- rbind(esc2016_2018_slow_rb_sn_ng,esc2016_2018_fast_rb_sn_ng)
esc2016_2018_rb_sn_ng <- na.omit(esc2016_2018_rb_sn_ng)

nrow(esc2016_2018_rb_sn_ng)
head(esc2016_2018_rb_sn_ng)


# esc2016_2018_rb_sn_ng[match(fib2016_2018_rb_sn_ng$NAME,esc2016_2018_rb_sn_ng$NAME),]

# esc add fib special gene sets
fib_s <- fib2016_2018_rb_sn_ng[! fib2016_2018_rb_sn_ng$NAME %in% esc2016_2018_rb_sn_ng$NAME,]
fib_s[,c(3,4,5)] <- 0
head(fib_s)
esc2016_2018_rb_sn_ng <- rbind(esc2016_2018_rb_sn_ng,fib_s)
nrow(esc2016_2018_rb_sn_ng)


# fib add fib special gene sets
esc_s <- esc2016_2018_rb_sn_ng[! esc2016_2018_rb_sn_ng$NAME %in% fib2016_2018_rb_sn_ng$NAME,]
esc_s[,c(3,4,5)] <- 0
head(esc_s)
fib2016_2018_rb_sn_ng <- rbind(fib2016_2018_rb_sn_ng,esc_s)
nrow(fib2016_2018_rb_sn_ng)




# merge all fib and esc

fib_esc <- cbind(fib2016_2018_rb_sn_ng,esc2016_2018_rb_sn_ng[match(fib2016_2018_rb_sn_ng$NAME,esc2016_2018_rb_sn_ng$NAME),])

nrow(fib_esc)



# special to esc or fib was removed
fib_esc <- na.omit(fib_esc)



head(fib_esc)
class(fib_esc)
colnames(fib_esc) <- c("NAME.FIB","SIZE_FIB","NES_FIB","NOM.p.val","FDR.q.val","NAME.ESC","SIZE_ESC","NES_ESC","NOM.p.val","FDR.q.val")
# fib_esc$top_pos <- fib_esc$NES_FIB+fib_esc$NES_ESC
# fib_esc$top_neg <- fib_esc$NES_FIB-fib_esc$NES_ESC


# get categories
category <- read.csv("NES_all_categories.csv",header = T)
head(category)
length(category[,3][-1])

fib_esc$category <- category[,3][-1]
head(fib_esc)


fib_esc$fib_add_esc <- fib_esc$NES_FIB +fib_esc$NES_ESC
fib_esc$fib_minus_esc <- fib_esc$NES_FIB-fib_esc$NES_ESC
head(fib_esc)

fib_esc <- fib_esc[,c(1,11,2,3,4,5,6,7,8,9,10,12,13)]

write.table(fib_esc,file="NES_all.tsv",col.names = T,row.names = F,quote = F,sep = "\t")
``