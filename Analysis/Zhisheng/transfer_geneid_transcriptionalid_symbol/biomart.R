
#获取小鼠基因的注释表
#并且做了3个基因表达量的变化-------------18.12.13


require("biomaRt")
listMarts()
datasets <- listDatasets()
mart <- useMart("ENSEMBL_MART_ENSEMBL")
listmart<-listDatasets(mart)
View(listmart)
mart <- useDataset("mmusculus_gene_ensembl", mart)

esc1 <- read.delim("RNASeqESC_Replicate1.tab.tsv",header=T,row.names=1,check.names=FALSE,sep=" ")


esc1ens <- rownames(esc1)

esc1Lookup <- gsub("\\.[0-9]*$", "", esc1ens)

annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_transcript_id","ensembl_gene_id", "gene_biotype", "external_gene_name"),
  filter="ensembl_transcript_id",
  values=esc1Lookup,
  uniqueRows=TRUE)



annotLookup2 <- data.frame(
  esc1ens[match(annotLookup$ensembl_transcript_id, esc1Lookup)],
  annotLookup)

colnames(annotLookup2) <- c(
  "original_id",
  c("ensembl_transcript_id","ensembl_gene_id", "gene_biotype", "external_gene_name"))

write.table(annotLookup2,file = "./ESC1.table")

annotLookup2 <- read.csv("ESC1.table",sep=" ")
head(df)


esc2 <- read.delim("RNASeqESC_Replicate2.tab.tsv",header=T,row.names=1,check.names=FALSE,sep=" ")
fib1 <- read.delim("RNASeqFIB_Replicate1.tab.tsv",header=T,row.names=1,check.names=FALSE,sep=" ")
fib2 <- read.delim("RNASeqFIB_Replicate2.tab.tsv",header=T,row.names=1,check.names=FALSE,sep=" ")


#读取基因对应的那一行的ensembl id
tars <- annotLookup2[annotLookup2$external_gene_name=="Tars",]
pmaip1 <- annotLookup2[annotLookup2$external_gene_name=="Pmaip1",]
tnfrsf21 <- annotLookup2[annotLookup2$external_gene_name=="Tnfrsf21",]

#将对应基因（转录本）RNA-seq data 的数据读入
Tars_esc1 <- esc1[rownames(esc1)==tars$original_id,]
Pmaip1_esc1 <- esc1[rownames(esc1)==pmaip1$original_id,]
Tnfrsf21_esc1 <- esc1[rownames(esc1)==tnfrsf21$original_id,]

colnames(Tars_esc1) <- c("slow","med","fast")
colnames(Pmaip1_esc1) <- c("slow","med","fast")
colnames(Tnfrsf21_esc1) <- c("slow","med","fast")
rownames(Tars_esc1) <- "Tars_esc1"
rownames(Pmaip1_esc1) <- "Pmaip1_esc1"
rownames(Tnfrsf21_esc1) <- "Tnfrsf21_esc1"


Tars_esc2 <- esc2[rownames(esc2)==tars$original_id,]
Pmaip1_esc2 <- esc2[rownames(esc2)==pmaip1$original_id,]
Tnfrsf21_esc2 <- esc2[rownames(esc2)==tnfrsf21$original_id,]

colnames(Tars_esc2) <- c("slow","med","fast")
colnames(Pmaip1_esc2) <- c("slow","med","fast")
colnames(Tnfrsf21_esc2) <- c("slow","med","fast")
rownames(Tars_esc2) <- "Tars_esc2"
rownames(Pmaip1_esc2) <- "Pmaip1_esc2"
rownames(Tnfrsf21_esc2) <- "Tnfrsf21_esc2"

Tars_fib1 <- fib1[rownames(fib1)==tars$original_id,]
Pmaip1_fib1 <- fib1[rownames(fib1)==pmaip1$original_id,]
Tnfrsf21_fib1 <- fib1[rownames(fib1)==tnfrsf21$original_id,]

colnames(Tars_fib1) <- c("slow","med","fast")
colnames(Pmaip1_fib1) <- c("slow","med","fast")
colnames(Tnfrsf21_fib1) <- c("slow","med","fast")
rownames(Tars_fib1) <- "Tars_fib1"
rownames(Pmaip1_fib1) <- "Pmaip1_fib1"
rownames(Tnfrsf21_fib1) <- "Tnfrsf21_fib1"

Tars_fib2 <- fib2[rownames(fib2)==tars$original_id,]
Pmaip1_fib2 <- fib2[rownames(fib2)==pmaip1$original_id,]
Tnfrsf21_fib2 <- fib2[rownames(fib2)==tnfrsf21$original_id,]

colnames(Tars_fib2) <- c("slow","med","fast")
colnames(Pmaip1_fib2) <- c("slow","med","fast")
colnames(Tnfrsf21_fib2) <- c("slow","med","fast")
rownames(Tars_fib2) <- "Tars_fib2"
rownames(Pmaip1_fib2) <- "Pmaip1_fib2"
rownames(Tnfrsf21_fib2) <- "Tnfrsf21_fib2"

require(ggplot2)
library(reshape2)

#-----------------------Tars-------------------------------


Tars_all <- rbind(Tars_esc1,Tars_esc2,Tars_fib1,Tars_fib2)

rownames(Tars_all) <- NULL

Tars_all <- log(Tars_all,base=2)
  
Tars_all <- data.frame(rep=c(1,1,2,2),name=c("Tars_esc1","Tars_esc2","Tars_fib1","Tars_fib2"),Tars_all)


Tars_all.m <- melt(Tars_all,id.vars = c("name","rep"))


p1 <- ggplot(Tars_all.m,aes(x =as.factor(Tars_all.m$variable), y = Tars_all.m$value,group=Tars_all.m$name, colour=as.factor(Tars_all.m$rep))) +
  geom_point(size = 2) +
  
  geom_line(size=1) +
  
  labs(title = "Tars") +
  
  scale_colour_discrete(name  ="Cell Type",
                        breaks=c("1", "2"),
                        labels=c("ESC", "FIB")) +
  #去除图例标题
  guides(fill=guide_legend(title="Cell Type")) +
  #主题
  theme(
    plot.title = element_text(lineheight=.8, size=20,face="bold",hjust = 0.5),
    #delete background
    panel.grid.major =element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    
    #加上坐标轴
    axis.line = element_line(colour = "black",size=0.5),
    #刻度线
    axis.ticks = element_line(size=0.5),
    axis.ticks.length=unit(0.15,"cm"),
    
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    #x轴标签
    axis.text.x = element_text(size=14,color='black'),
    #y轴标签
    axis.text.y = element_text(size=14,color='black'),
    #图例
    legend.title = element_text(colour="bLACK", size=14, face="bold")
    )
   

#---------------------------------------

Pmaip1_all <- rbind(Pmaip1_esc1,Pmaip1_esc2,Pmaip1_fib1,Pmaip1_fib2)

rownames(Pmaip1_all) <- NULL

Pmaip1_all <- log(Pmaip1_all,base=2)

Pmaip1_all <- data.frame(rep=c(1,1,2,2),name=c("Pmaip1_esc1","Pmaip1_esc2","Pmaip1_fib1","Pmaip1_fib2"),Pmaip1_all)


Pmaip1_all.m <- melt(Pmaip1_all,id.vars = c("name","rep"))


p2 <- ggplot(Pmaip1_all.m,aes(x =as.factor(Pmaip1_all.m$variable), y = Pmaip1_all.m$value,group=Pmaip1_all.m$name, colour=as.factor(Pmaip1_all.m$rep))) +
  geom_point(size = 2) +
  
  geom_line(size=1) +
  
  labs(title = "Pmaip1") +
  
  scale_colour_discrete(name  ="Cell Type",
                        breaks=c("1", "2"),
                        labels=c("ESC", "FIB")) +
  #去除图例标题
  guides(fill=guide_legend(title="Cell Type")) +
  #主题
  theme(
    plot.title = element_text(lineheight=.8, size=20,face="bold",hjust = 0.5),
    #delete background
    panel.grid.major =element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    
    #加上坐标轴
    axis.line = element_line(colour = "black",size=0.5),
    #刻度线
    axis.ticks = element_line(size=0.5),
    axis.ticks.length=unit(0.15,"cm"),
    
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    #x轴标签
    axis.text.x = element_text(size=14,color='black'),
    #y轴标签
    axis.text.y = element_text(size=14,color='black'),
    #图例
    legend.title = element_text(colour="bLACK", size=14, face="bold")
  )

p

#--------------------------------

Tnfrsf21_all <- rbind(Tnfrsf21_esc1,Tnfrsf21_esc2,Tnfrsf21_fib1,Tnfrsf21_fib2)

rownames(Tnfrsf21_all) <- NULL

Tnfrsf21_all <- log(Tnfrsf21_all,base=2)

Tnfrsf21_all <- data.frame(rep=c(1,1,2,2),name=c("Tnfrsf21_esc1","Tnfrsf21_esc2","Tnfrsf21_fib1","Tnfrsf21_fib2"),Tnfrsf21_all)


Tnfrsf21_all.m <- melt(Tnfrsf21_all,id.vars = c("name","rep"))


p3 <- ggplot(Tnfrsf21_all.m,aes(x =as.factor(Tnfrsf21_all.m$variable), y = Tnfrsf21_all.m$value,group=Tnfrsf21_all.m$name, colour=as.factor(Tnfrsf21_all.m$rep))) +
  geom_point(size = 2) +
  
  geom_line(size=1) +
  
  labs(title = "Tnfrsf21") +
  
  scale_colour_discrete(name  ="Cell Type",
                        breaks=c("1", "2"),
                        labels=c("ESC", "FIB")) +
  #去除图例标题
  guides(fill=guide_legend(title="Cell Type")) +
  #主题
  theme(
    plot.title = element_text(lineheight=.8, size=20,face="bold",hjust = 0.5),
    #delete background
    panel.grid.major =element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    
    #加上坐标轴
    axis.line = element_line(colour = "black",size=0.5),
    #刻度线
    axis.ticks = element_line(size=0.5),
    axis.ticks.length=unit(0.15,"cm"),
    
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    #x轴标签
    axis.text.x = element_text(size=14,color='black'),
    #y轴标签
    axis.text.y = element_text(size=14,color='black'),
    #图例
    legend.title = element_text(colour="bLACK", size=14, face="bold")
  )


p









#---------------------------------------------------
p <- p + geom_point()

###
library(plotly)
###
p1
ggplotly(p1)
  
p
##########
  
  
  
  
  
  geom_errorbar(aes(ymin = ymin, ymax = ymax))








errbar(x, y, yplus, yminus, cap, xlab, ylab, add=FALSE, lty=1, ylim, lwd=1, Type=rep(1,length(y)), ... )