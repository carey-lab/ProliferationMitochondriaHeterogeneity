
#将基因分成不同的基因sets然后画图--18.12.23

annotLookup2 <- read.table("annotation.table",header = T)

human_gene <- read.table("mouse2human_ensembl_gene_id.table",header = T)

fib1 <- read.delim("RNASeqFIB_Replicate1.tab.tsv",header=T,row.names=1,check.names=FALSE,sep=" ")
fib2 <- read.delim("RNASeqFIB_Replicate2.tab.tsv",header=T,row.names=1,check.names=FALSE,sep=" ")
esc1 <- read.delim("RNASeqESC_Replicate1.tab.tsv",header=T,row.names=1,check.names=FALSE,sep=" ")
esc2 <- read.delim("RNASeqESC_Replicate2.tab.tsv",header=T,row.names=1,check.names=FALSE,sep=" ")


#------------------DNA_repair-----------------------------------------

DNA_repair <- read.table("../GO_data/HALLMARK_DNA_REPAIR.txt",header = T, sep="\t")


#有一些NA值，证明一些human_DNA_repair_gene 没有对应的老鼠基因
mouse_DNA_repair_gid <- human_gene$mouse_ensembl_gene_id[match(DNA_repair[2:nrow(DNA_repair),],human_gene$human_external_gene_name)]


mouse_DNA_repair_oid <- as.character(annotLookup2$original_id[match(mouse_DNA_repair_gid,annotLookup2$ensembl_gene_id)])



#载入fib1
temp=match(mouse_DNA_repair_oid,row.names(fib1))
temp=temp[!is.na(temp)]
mouse_DNA_repair_exp1 <- fib1[temp,]
#载入fib2
temp=match(mouse_DNA_repair_oid,row.names(fib2))
temp=temp[!is.na(temp)]
mouse_DNA_repair_exp2 <- fib2[temp,]
#import esc1
temp=match(mouse_DNA_repair_oid,row.names(esc1))
temp=temp[!is.na(temp)]
mouse_DNA_repair_exp3 <- esc1[temp,]
#载入esc2
temp=match(mouse_DNA_repair_oid,row.names(esc2))
temp=temp[!is.na(temp)]
mouse_DNA_repair_exp4 <- esc2[temp,]




# FIB_slow_TPM  FIB_med_TPM FIB_fast_TPM 
# 39.37212     37.07594     36.16368 
#plot it
library(ggplot2)
library(reshape2)


#mouse_DNA_repair_exp <- log(mouse_DNA_repair_exp,2)
rownames(mouse_DNA_repair_exp1) <- c()

colnames(mouse_DNA_repair_exp2) <- c("FIB_slow_TPM","FIB_med_TPM" ,"FIB_fast_TPM")
rownames(mouse_DNA_repair_exp2) <- c()
rownames(mouse_DNA_repair_exp3) <- c()
rownames(mouse_DNA_repair_exp4) <- c()




mouse_DNA_repair_exp1.m <- melt(mouse_DNA_repair_exp1)
mouse_DNA_repair_exp2.m <- melt(mouse_DNA_repair_exp2)
mouse_DNA_repair_exp3.m <- melt(mouse_DNA_repair_exp3)
mouse_DNA_repair_exp4.m <- melt(mouse_DNA_repair_exp4)


mouse_DNA_repair_exp1.m$value <- log(mouse_DNA_repair_exp1.m$value+0.1,2)
mouse_DNA_repair_exp2.m$value <- log(mouse_DNA_repair_exp2.m$value+0.1,2)
mouse_DNA_repair_exp3.m$value <- log(mouse_DNA_repair_exp3.m$value+0.1,2)
mouse_DNA_repair_exp4.m$value <- log(mouse_DNA_repair_exp4.m$value+0.1,2)




head(mouse_DNA_repair_exp1.m)



#plot fib1 violin plot-------------

ggplot(mouse_DNA_repair_exp1.m,aes(
  as.factor(mouse_DNA_repair_exp1.m$variable),
  mouse_DNA_repair_exp1.m$value,
  color=as.factor(mouse_DNA_repair_exp1.m$variable)
  )) +
  
  geom_violin() +
  
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="---", size= 6, color= "black") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="FIB1") +


  
  labs(title = "DNA_repair", x="", y="Expression log2(TPM)") +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

#plot fib2 violin plot-------------

ggplot(mouse_DNA_repair_exp2.m,aes(
  as.factor(mouse_DNA_repair_exp2.m$variable),
  mouse_DNA_repair_exp2.m$value,
  color=as.factor(mouse_DNA_repair_exp2.m$variable)
)) +
  
  geom_violin() +
  
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="---", size= 6, color= "black") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="FIB2") +
  
  
  
  labs(title = "DNA_repair", x="", y="Expression log2(TPM)") +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))


#plot ESC1 violin plot-------------

ggplot(mouse_DNA_repair_exp3.m,aes(
  as.factor(mouse_DNA_repair_exp3.m$variable),
  mouse_DNA_repair_exp3.m$value,
  color=as.factor(mouse_DNA_repair_exp3.m$variable)
)) +
  
  geom_violin() +
  
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="---", size= 6, color= "black") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="ESC1") +
  
  
  
  labs(title = "DNA_repair", x="", y="Expression log2(TPM)") +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

#plot ESC2 violin plot-------------

ggplot(mouse_DNA_repair_exp4.m,aes(
  as.factor(mouse_DNA_repair_exp4.m$variable),
  mouse_DNA_repair_exp4.m$value,
  color=as.factor(mouse_DNA_repair_exp4.m$variable)
)) +
  
  geom_violin() +
  
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="---", size= 6, color= "black") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="ESC2") +
  
  
  
  labs(title = "DNA_repair", x="", y="Expression log2(TPM)") +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

#-----------------------GO RIBOSOME-----------------------------------------

GO_RIBOSOME <- read.table("../GO_data/GO0005840GO_RIBOSOME.txt",header = T, sep="\t")

#有一些NA值，证明一些human_DNA_repair_gene 没有对应的老鼠基因
mouse_GO_RIBOSOME_gid <- human_gene$mouse_ensembl_gene_id[match(GO_RIBOSOME[2:nrow(GO_RIBOSOME),],human_gene$human_external_gene_name)]


mouse_GO_RIBOSOME_oid <- as.character(annotLookup2$original_id[match(mouse_GO_RIBOSOME_gid,annotLookup2$ensembl_gene_id)])



#载入fib1
temp=match(mouse_GO_RIBOSOME_oid,row.names(fib1))
temp=temp[!is.na(temp)]
mouse_GO_RIBOSOME_exp1 <- fib1[temp,]
#载入fib2
temp=match(mouse_GO_RIBOSOME_oid,row.names(fib2))
temp=temp[!is.na(temp)]
mouse_GO_RIBOSOME_exp2 <- fib2[temp,]
#import esc1
temp=match(mouse_GO_RIBOSOME_oid,row.names(esc1))
temp=temp[!is.na(temp)]
mouse_GO_RIBOSOME_exp3 <- esc1[temp,]
#载入esc2
temp=match(mouse_GO_RIBOSOME_oid,row.names(esc2))
temp=temp[!is.na(temp)]
mouse_GO_RIBOSOME_exp4 <- esc2[temp,]




# FIB_slow_TPM  FIB_med_TPM FIB_fast_TPM 
# 39.37212     37.07594     36.16368 
#plot it


#mouse_GO_RIBOSOME_exp <- log(mouse_GO_RIBOSOME_exp,2)
rownames(mouse_GO_RIBOSOME_exp1) <- c()

colnames(mouse_GO_RIBOSOME_exp2) <- c("FIB_slow_TPM","FIB_med_TPM" ,"FIB_fast_TPM")
rownames(mouse_GO_RIBOSOME_exp2) <- c()
rownames(mouse_GO_RIBOSOME_exp3) <- c()
rownames(mouse_GO_RIBOSOME_exp4) <- c()




mouse_GO_RIBOSOME_exp1.m <- melt(mouse_GO_RIBOSOME_exp1)
mouse_GO_RIBOSOME_exp2.m <- melt(mouse_GO_RIBOSOME_exp2)
mouse_GO_RIBOSOME_exp3.m <- melt(mouse_GO_RIBOSOME_exp3)
mouse_GO_RIBOSOME_exp4.m <- melt(mouse_GO_RIBOSOME_exp4)


mouse_GO_RIBOSOME_exp1.m$value <- log(mouse_GO_RIBOSOME_exp1.m$value+0.1,2)
mouse_GO_RIBOSOME_exp2.m$value <- log(mouse_GO_RIBOSOME_exp2.m$value+0.1,2)
mouse_GO_RIBOSOME_exp3.m$value <- log(mouse_GO_RIBOSOME_exp3.m$value+0.1,2)
mouse_GO_RIBOSOME_exp4.m$value <- log(mouse_GO_RIBOSOME_exp4.m$value+0.1,2)

colMeans(mouse_GO_RIBOSOME_exp2)
#-----plot fib1-------------------
ggplot(mouse_GO_RIBOSOME_exp1.m,aes(
  as.factor(mouse_GO_RIBOSOME_exp1.m$variable),
  mouse_GO_RIBOSOME_exp1.m$value,
  color=as.factor(mouse_GO_RIBOSOME_exp1.m$variable)
)) +
  
  geom_violin() +
  
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  
  stat_summary(fun.y = "mean", geom = "text", label="---", size= 6, color= "black") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="FIB1") +
  
  
  
  labs(title = "GO_RIBOSOME", x="", y="Expression log2(TPM)") +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

#plot fib2 violin plot-------------

ggplot(mouse_GO_RIBOSOME_exp2.m,aes(
  as.factor(mouse_GO_RIBOSOME_exp2.m$variable),
  mouse_GO_RIBOSOME_exp2.m$value,
  color=as.factor(mouse_GO_RIBOSOME_exp2.m$variable)
)) +
  
  geom_violin() +
  
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="---", size= 6, color= "black") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="FIB2") +
  
  
  
  labs(title = "GO_RIBOSOME", x="", y="Expression log2(TPM)") +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))


#plot ESC1 violin plot-------------

ggplot(mouse_GO_RIBOSOME_exp3.m,aes(
  as.factor(mouse_GO_RIBOSOME_exp3.m$variable),
  mouse_GO_RIBOSOME_exp3.m$value,
  color=as.factor(mouse_GO_RIBOSOME_exp3.m$variable)
)) +
  
  geom_violin() +
  
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="---", size= 6, color= "black") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="ESC1") +
  
  
  
  labs(title = "GO_RIBOSOME", x="", y="Expression log2(TPM)") +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

#plot ESC2 violin plot-------------

ggplot(mouse_GO_RIBOSOME_exp4.m,aes(
  as.factor(mouse_GO_RIBOSOME_exp4.m$variable),
  mouse_GO_RIBOSOME_exp4.m$value,
  color=as.factor(mouse_GO_RIBOSOME_exp4.m$variable)
)) +
  
  geom_violin() +
  
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="---", size= 6, color= "black") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="ESC2") +
  
  
  
  labs(title = "GO_RIBOSOME", x="", y="Expression log2(TPM)") +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))


#-----------------------------RESPONSE_TO_DNA_DAMAGE_STIMULUS-------------

RESPONSE_TO_DNA_DAMAGE_STIMULUS <- read.table("../GO_data/GO0005840RESPONSE_TO_DNA_DAMAGE_STIMULUS.txt",header = T, sep="\t")

#有一些NA值，证明一些human_DNA_repair_gene 没有对应的老鼠基因
mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_gid <- human_gene$mouse_ensembl_gene_id[match(DNA_repair[2:nrow(DNA_repair),],human_gene$human_external_gene_name)]


mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_oid <- as.character(annotLookup2$original_id[match(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_gid,annotLookup2$ensembl_gene_id)])



#载入fib1
temp=match(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_oid,row.names(fib1))
temp=temp[!is.na(temp)]
mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1 <- fib1[temp,]
#载入fib2
temp=match(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_oid,row.names(fib2))
temp=temp[!is.na(temp)]
mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2 <- fib2[temp,]
#import esc1
temp=match(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_oid,row.names(esc1))
temp=temp[!is.na(temp)]
mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3 <- esc1[temp,]
#载入esc2
temp=match(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_oid,row.names(esc2))
temp=temp[!is.na(temp)]
mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4 <- esc2[temp,]




# FIB_slow_TPM  FIB_med_TPM FIB_fast_TPM 
# 39.37212     37.07594     36.16368 
#plot it


#mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp <- log(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp,2)
rownames(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1) <- c()

colnames(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2) <- c("FIB_slow_TPM","FIB_med_TPM" ,"FIB_fast_TPM")
rownames(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2) <- c()
rownames(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3) <- c()
rownames(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4) <- c()




mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1.m <- melt(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1)
mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2.m <- melt(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2)
mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3.m <- melt(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3)
mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4.m <- melt(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4)


mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1.m$value <- log(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1.m$value+0.1,2)
mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2.m$value <- log(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2.m$value+0.1,2)
mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3.m$value <- log(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3.m$value+0.1,2)
mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4.m$value <- log(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4.m$value+0.1,2)


#-----plot fib1-------------------
ggplot(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1.m,aes(
  as.factor(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1.m$variable),
  mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1.m$value,
  color=as.factor(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1.m$variable)
)) +
  
  geom_violin() +
  
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="---", size= 6, color= "black") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="FIB1") +
  
  
  
  labs(title = "RESPONSE_TO_DNA_DAMAGE_STIMULUS", x="", y="Expression log2(TPM)") +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

#plot fib2 violin plot-------------

ggplot(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2.m,aes(
  as.factor(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2.m$variable),
  mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2.m$value,
  color=as.factor(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2.m$variable)
)) +
  
  geom_violin() +
  
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="---", size= 6, color= "black") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="FIB2") +
  
  
  
  labs(title = "RESPONSE_TO_DNA_DAMAGE_STIMULUS", x="", y="Expression log2(TPM)") +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))


#plot ESC1 violin plot-------------

ggplot(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3.m,aes(
  as.factor(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3.m$variable),
  mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3.m$value,
  color=as.factor(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3.m$variable)
)) +
  
  geom_violin() +
  
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="---", size= 6, color= "black") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="ESC1") +
  
  
  
  labs(title = "RESPONSE_TO_DNA_DAMAGE_STIMULUS", x="", y="Expression log2(TPM)") +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

#plot ESC2 violin plot-------------

ggplot(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4.m,aes(
  as.factor(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4.m$variable),
  mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4.m$value,
  color=as.factor(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4.m$variable)
)) +
  
  geom_violin() +
  
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="---", size= 6, color= "black") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="ESC2") +
  
  
  
  labs(title = "RESPONSE_TO_DNA_DAMAGE_STIMULUS", x="", y="Expression log2(TPM)") +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))


#----------------------------------------

REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS <- read.table("../GO_data/GO0005840REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS.txt",header = T, sep="\t")

#有一些NA值，证明一些human_DNA_repair_gene 没有对应的老鼠基因
mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_gid <- human_gene$mouse_ensembl_gene_id[match(DNA_repair[2:nrow(DNA_repair),],human_gene$human_external_gene_name)]


mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_oid <- as.character(annotLookup2$original_id[match(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_gid,annotLookup2$ensembl_gene_id)])



#载入fib1
temp=match(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_oid,row.names(fib1))
temp=temp[!is.na(temp)]
mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1 <- fib1[temp,]
#载入fib2
temp=match(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_oid,row.names(fib2))
temp=temp[!is.na(temp)]
mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2 <- fib2[temp,]
#import esc1
temp=match(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_oid,row.names(esc1))
temp=temp[!is.na(temp)]
mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3 <- esc1[temp,]
#载入esc2
temp=match(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_oid,row.names(esc2))
temp=temp[!is.na(temp)]
mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4 <- esc2[temp,]




# FIB_slow_TPM  FIB_med_TPM FIB_fast_TPM 
# 39.37212     37.07594     36.16368 
#plot it


#mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp <- log(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp,2)
rownames(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1) <- c()

colnames(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2) <- c("FIB_slow_TPM","FIB_med_TPM" ,"FIB_fast_TPM")
rownames(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2) <- c()
rownames(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3) <- c()
rownames(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4) <- c()




mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1.m <- melt(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1)
mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2.m <- melt(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2)
mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3.m <- melt(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3)
mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4.m <- melt(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4)


mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1.m$value <- log(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1.m$value+0.1,2)
mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2.m$value <- log(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2.m$value+0.1,2)
mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3.m$value <- log(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3.m$value+0.1,2)
mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4.m$value <- log(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4.m$value+0.1,2)

#-----plot fib1-------------------
ggplot(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1.m,aes(
  as.factor(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1.m$variable),
  mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1.m$value,
  color=as.factor(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1.m$variable)
)) +
  
  geom_violin() +
  
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="---", size= 6, color= "black") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="FIB1") +
  
  
  
  labs(title = "REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS", x="", y="Expression log2(TPM)") +
  
  theme(
    plot.title = element_text(lineheight=.8, size=12,face="bold",hjust = 0.5))

#plot fib2 violin plot-------------

ggplot(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2.m,aes(
  as.factor(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2.m$variable),
  mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2.m$value,
  color=as.factor(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2.m$variable)
)) +
  
  geom_violin() +
  
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="---", size= 6, color= "black") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="FIB2") +
  
  
  
  labs(title = "REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS", x="", y="Expression log2(TPM)") +
  
  theme(
    plot.title = element_text(lineheight=.8, size=12,face="bold",hjust = 0.5))


#plot ESC1 violin plot-------------

ggplot(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3.m,aes(
  as.factor(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3.m$variable),
  mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3.m$value,
  color=as.factor(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3.m$variable)
)) +
  
  geom_violin() +
  
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="---", size= 6, color= "black") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="ESC1") +
  
  
  
  labs(title = "REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS", x="", y="Expression log2(TPM)") +
  
  theme(
    plot.title = element_text(lineheight=.8, size=12,face="bold",hjust = 0.5))

#plot ESC2 violin plot-------------

ggplot(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4.m,aes(
  as.factor(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4.m$variable),
  mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4.m$value,
  color=as.factor(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4.m$variable)
)) +
  
  geom_violin() +
  
  scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="---", size= 6, color= "black") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="ESC2") +
  
  
  
  labs(title = "REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS", x="", y="Expression log2(TPM)") +
  
  theme(
    plot.title = element_text(lineheight=.8, size=12,face="bold",hjust = 0.5))


