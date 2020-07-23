
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



slope_func <- function(exp.table) {
  
  x <- c(1,2,3)
  slope <- c(rep(0,nrow(exp.table)))
  for (i in 1:nrow(exp.table)) {
    y <- as.numeric(exp.table[i,]);
    t <- lm(y~x)
    slope[i] <- t$coefficients[[2]]
  }
  
  return(slope)
  }

#---------------slope------------------------------















#log2(fast/slow)------------------------------------
#fib1
#除去0，都加0.1
#mouse_DNA_repair_exp1_slope <- mouse_DNA_repair_exp1+0.1

#mouse_DNA_repair_exp1_f2s <- log(mouse_DNA_repair_exp1$FIB_fast_TPM/mouse_DNA_repair_exp1$FIB_slow_TPM,2)




mouse_DNA_repair_exp1_c <- log(mouse_DNA_repair_exp1+0.1,2)
mouse_DNA_repair_exp1_slope <- slope_func(mouse_DNA_repair_exp1_c)
mouse_DNA_repair_exp1_slope <- as.data.frame(mouse_DNA_repair_exp1_slope)




library(ggplot2)
library(reshape2)
require(gridExtra)

ggplot(mouse_DNA_repair_exp1_slope,aes(
  x="",
  y=mouse_DNA_repair_exp1_slope
  #color=as.factor(mouse_DNA_repair_exp1.m$variable)
)) +
  
  geom_violin(fill = "grey80",colour = "#3366FF") +
  
  #scale_fill_brewer(palette="Dark2")
  
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="------", size= 6, color= "red") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="FIB1") +
  
  labs(title = "DNA_repair", x="log2FIB1", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1226/DNA_repair_fib1_log_slope.png", width = 2, height = 4)


#FIB2------------------------------------

mouse_DNA_repair_exp2_c <- log(mouse_DNA_repair_exp2+0.1,2)
mouse_DNA_repair_exp2_slope <- slope_func(mouse_DNA_repair_exp2_c)
mouse_DNA_repair_exp2_slope <- as.data.frame(mouse_DNA_repair_exp2_slope)




library(ggplot2)
library(reshape2)
require(gridExtra)

ggplot(mouse_DNA_repair_exp2_slope,aes(
  x="",
  y=mouse_DNA_repair_exp2_slope
  #color=as.factor(mouse_DNA_repair_exp2.m$variable)
)) +
  
  geom_violin(fill = "grey80",colour = "#3366FF") +
  
  #scale_fill_brewer(palette="Dark2")
  
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="------", size= 6, color= "red") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="FIB2") +
  
  labs(title = "DNA_repair", x="log2FIB2", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1226/DNA_repair_FIB2_log_slope.png", width = 2, height = 4)




#----------esc1----------------------

mouse_DNA_repair_exp3_c <- log(mouse_DNA_repair_exp3+0.1,2)
mouse_DNA_repair_exp3_slope <- slope_func(mouse_DNA_repair_exp3_c)
mouse_DNA_repair_exp3_slope <- as.data.frame(mouse_DNA_repair_exp3_slope)




library(ggplot2)
library(reshape2)
require(gridExtra)

ggplot(mouse_DNA_repair_exp3_slope,aes(
  x="",
  y=mouse_DNA_repair_exp3_slope
  #color=as.factor(mouse_DNA_repair_exp3.m$variable)
)) +
  
  geom_violin(fill = "grey80",colour = "#3366FF") +
  
  #scale_fill_brewer(palette="Dark2")
  
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="------", size= 6, color= "red") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="ESC1") +
  
  labs(title = "DNA_repair", x="log2ESC1", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1226/DNA_repair_ESC1_log_slope.png", width = 2, height = 4)



#----------------esc2--------------------
mouse_DNA_repair_exp4_c <- log(mouse_DNA_repair_exp4+0.1,2)
mouse_DNA_repair_exp4_slope <- slope_func(mouse_DNA_repair_exp4_c)
mouse_DNA_repair_exp4_slope <- as.data.frame(mouse_DNA_repair_exp4_slope)




library(ggplot2)
library(reshape2)
require(gridExtra)

ggplot(mouse_DNA_repair_exp4_slope,aes(
  x="",
  y=mouse_DNA_repair_exp4_slope
  #color=as.factor(mouse_DNA_repair_exp4.m$variable)
)) +
  
  geom_violin(fill = "grey80",colour = "#3366FF") +
  
  #scale_fill_brewer(palette="Dark2")
  
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="------", size= 6, color= "red") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="ESC2") +
  
  labs(title = "DNA_repair", x="log2ESC2", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1226/DNA_repair_ESC2_log_slope.png", width = 2, height = 4)



#----------go-ribosome----------------------------

GO_RIBOSOME <- read.table("../GO_data/GO0005840GO_RIBOSOME.txt",header = T, sep="\t")


#有一些NA值，证明一些human_GO_RIBOSOME_gene 没有对应的老鼠基因
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






#fib1
#除去0，都加0.1
#mouse_GO_RIBOSOME_exp1_slope <- mouse_GO_RIBOSOME_exp1+0.1

#mouse_GO_RIBOSOME_exp1_f2s <- log(mouse_GO_RIBOSOME_exp1$FIB_fast_TPM/mouse_GO_RIBOSOME_exp1$FIB_slow_TPM,2)


mouse_GO_RIBOSOME_exp1_c <- log(mouse_GO_RIBOSOME_exp1+0.1,2)
mouse_GO_RIBOSOME_exp1_slope <- slope_func(mouse_GO_RIBOSOME_exp1_c)
mouse_GO_RIBOSOME_exp1_slope <- as.data.frame(mouse_GO_RIBOSOME_exp1_slope)




library(ggplot2)
library(reshape2)
require(gridExtra)

ggplot(mouse_GO_RIBOSOME_exp1_slope,aes(
  x="",
  y=mouse_GO_RIBOSOME_exp1_slope
  #color=as.factor(mouse_GO_RIBOSOME_exp1.m$variable)
)) +
  
  geom_violin(fill = "grey80",colour = "#3366FF") +
  
  #scale_fill_brewer(palette="Dark2")
  
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="------", size= 6, color= "red") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="FIB1") +
  
  labs(title = "GO_RIBOSOME", x="log2FIB1", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1226/GO_RIBOSOME_fib1_log_slope.png", width = 2, height = 4)


#FIB2------------------------------------

mouse_GO_RIBOSOME_exp2_c <- log(mouse_GO_RIBOSOME_exp2+0.1,2)
mouse_GO_RIBOSOME_exp2_slope <- slope_func(mouse_GO_RIBOSOME_exp2_c)
mouse_GO_RIBOSOME_exp2_slope <- as.data.frame(mouse_GO_RIBOSOME_exp2_slope)




library(ggplot2)
library(reshape2)
require(gridExtra)

ggplot(mouse_GO_RIBOSOME_exp2_slope,aes(
  x="",
  y=mouse_GO_RIBOSOME_exp2_slope
  #color=as.factor(mouse_GO_RIBOSOME_exp2.m$variable)
)) +
  
  geom_violin(fill = "grey80",colour = "#3366FF") +
  
  #scale_fill_brewer(palette="Dark2")
  
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="------", size= 6, color= "red") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="FIB2") +
  
  labs(title = "GO_RIBOSOME", x="log2FIB2", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1226/GO_RIBOSOME_FIB2_log_slope.png", width = 2, height = 4)




#----------esc1----------------------

mouse_GO_RIBOSOME_exp3_c <- log(mouse_GO_RIBOSOME_exp3+0.1,2)
mouse_GO_RIBOSOME_exp3_slope <- slope_func(mouse_GO_RIBOSOME_exp3_c)
mouse_GO_RIBOSOME_exp3_slope <- as.data.frame(mouse_GO_RIBOSOME_exp3_slope)




library(ggplot2)
library(reshape2)
require(gridExtra)

ggplot(mouse_GO_RIBOSOME_exp3_slope,aes(
  x="",
  y=mouse_GO_RIBOSOME_exp3_slope
  #color=as.factor(mouse_GO_RIBOSOME_exp3.m$variable)
)) +
  
  geom_violin(fill = "grey80",colour = "#3366FF") +
  
  #scale_fill_brewer(palette="Dark2")
  
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="------", size= 6, color= "red") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="ESC1") +
  
  labs(title = "GO_RIBOSOME", x="log2ESC1", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1226/GO_RIBOSOME_ESC1_log_slope.png", width = 2, height = 4)



#----------------esc2--------------------
mouse_GO_RIBOSOME_exp4_c <- log(mouse_GO_RIBOSOME_exp4+0.1,2)
mouse_GO_RIBOSOME_exp4_slope <- slope_func(mouse_GO_RIBOSOME_exp4_c)
mouse_GO_RIBOSOME_exp4_slope <- as.data.frame(mouse_GO_RIBOSOME_exp4_slope)




library(ggplot2)
library(reshape2)
require(gridExtra)

ggplot(mouse_GO_RIBOSOME_exp4_slope,aes(
  x="",
  y=mouse_GO_RIBOSOME_exp4_slope
  #color=as.factor(mouse_GO_RIBOSOME_exp4.m$variable)
)) +
  
  geom_violin(fill = "grey80",colour = "#3366FF") +
  
  #scale_fill_brewer(palette="Dark2")
  
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="------", size= 6, color= "red") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="ESC2") +
  
  labs(title = "GO_RIBOSOME", x="log2ESC2", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1226/GO_RIBOSOME_ESC2_log_slope.png", width = 2, height = 4)



#response to dna damage---------------------------


RESPONSE_TO_DNA_DAMAGE_STIMULUS <- read.table("../GO_data/GO0006974RESPONSE_TO_DNA_DAMAGE_STIMULUS.txt",header = T, sep="\t")


#有一些NA值，证明一些human_RESPONSE_TO_DNA_DAMAGE_STIMULUS_gene 没有对应的老鼠基因
mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_gid <- human_gene$mouse_ensembl_gene_id[match(RESPONSE_TO_DNA_DAMAGE_STIMULUS[2:nrow(RESPONSE_TO_DNA_DAMAGE_STIMULUS),],human_gene$human_external_gene_name)]


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






#log2(fast/slow)------------------------------------
#fib1
#除去0，都加0.1
#mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1_slope <- mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1+0.1

#mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1_f2s <- log(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1$FIB_fast_TPM/mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1$FIB_slow_TPM,2)


mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1_c <- log(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1+0.1,2)
mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1_slope <- slope_func(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1_c)
mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1_slope <- as.data.frame(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1_slope)




library(ggplot2)
library(reshape2)
require(gridExtra)

ggplot(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1_slope,aes(
  x="",
  y=mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1_slope
  #color=as.factor(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1.m$variable)
)) +
  
  geom_violin(fill = "grey80",colour = "#3366FF") +
  
  #scale_fill_brewer(palette="Dark2")
  
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="------", size= 6, color= "red") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="FIB1") +
  
  labs(title = "RESPONSE_TO_DNA_DAMAGE_STIMULUS", x="log2FIB1", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1226/RESPONSE_TO_DNA_DAMAGE_STIMULUS_fib1_log_slope.png", width = 2, height = 4)


#FIB2------------------------------------

mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2_c <- log(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2+0.1,2)
mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2_slope <- slope_func(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2_c)
mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2_slope <- as.data.frame(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2_slope)




library(ggplot2)
library(reshape2)
require(gridExtra)

ggplot(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2_slope,aes(
  x="",
  y=mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2_slope
  #color=as.factor(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2.m$variable)
)) +
  
  geom_violin(fill = "grey80",colour = "#3366FF") +
  
  #scale_fill_brewer(palette="Dark2")
  
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="------", size= 6, color= "red") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="FIB2") +
  
  labs(title = "RESPONSE_TO_DNA_DAMAGE_STIMULUS", x="log2FIB2", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1226/RESPONSE_TO_DNA_DAMAGE_STIMULUS_FIB2_log_slope.png", width = 2, height = 4)




#----------esc1----------------------

mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3_c <- log(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3+0.1,2)
mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3_slope <- slope_func(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3_c)
mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3_slope <- as.data.frame(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3_slope)




library(ggplot2)
library(reshape2)
require(gridExtra)

ggplot(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3_slope,aes(
  x="",
  y=mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3_slope
  #color=as.factor(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3.m$variable)
)) +
  
  geom_violin(fill = "grey80",colour = "#3366FF") +
  
  #scale_fill_brewer(palette="Dark2")
  
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="------", size= 6, color= "red") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="ESC1") +
  
  labs(title = "RESPONSE_TO_DNA_DAMAGE_STIMULUS", x="log2ESC1", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1226/RESPONSE_TO_DNA_DAMAGE_STIMULUS_ESC1_log_slope.png", width = 2, height = 4)



#----------------esc2--------------------
mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4_c <- log(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4+0.1,2)
mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4_slope <- slope_func(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4_c)
mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4_slope <- as.data.frame(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4_slope)




library(ggplot2)
library(reshape2)
require(gridExtra)

ggplot(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4_slope,aes(
  x="",
  y=mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4_slope
  #color=as.factor(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4.m$variable)
)) +
  
  geom_violin(fill = "grey80",colour = "#3366FF") +
  
  #scale_fill_brewer(palette="Dark2")
  
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="------", size= 6, color= "red") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="ESC2") +
  
  labs(title = "RESPONSE_TO_DNA_DAMAGE_STIMULUS", x="log2ESC2", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1226/RESPONSE_TO_DNA_DAMAGE_STIMULUS_ESC2_log_slope.png", width = 2, height = 4)



#-----------REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS--------------

REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS <- read.table("../GO_data/GO_2001020GO_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS.txt",header = T, sep="\t")


#有一些NA值，证明一些human_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_gene 没有对应的老鼠基因
mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_gid <- human_gene$mouse_ensembl_gene_id[match(REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS[2:nrow(REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS),],human_gene$human_external_gene_name)]


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






#log2(fast/slow)------------------------------------
#fib1
#除去0，都加0.1
#mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1_slope <- mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1+0.1

#mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1_f2s <- log(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1$FIB_fast_TPM/mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1$FIB_slow_TPM,2)


mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1_c <- log(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1+0.1,2)
mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1_slope <- slope_func(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1_c)
mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1_slope <- as.data.frame(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1_slope)




library(ggplot2)
library(reshape2)
require(gridExtra)

ggplot(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1_slope,aes(
  x="",
  y=mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1_slope
  #color=as.factor(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1.m$variable)
)) +
  
  geom_violin(fill = "grey80",colour = "#3366FF") +
  
  #scale_fill_brewer(palette="Dark2")
  
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="------", size= 6, color= "red") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="FIB1") +
  
  labs(title = "REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS", x="log2FIB1", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1226/REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_fib1_log_slope.png", width = 2, height = 4)


#FIB2------------------------------------

mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2_c <- log(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2+0.1,2)
mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2_slope <- slope_func(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2_c)
mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2_slope <- as.data.frame(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2_slope)




library(ggplot2)
library(reshape2)
require(gridExtra)

ggplot(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2_slope,aes(
  x="",
  y=mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2_slope
  #color=as.factor(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2.m$variable)
)) +
  
  geom_violin(fill = "grey80",colour = "#3366FF") +
  
  #scale_fill_brewer(palette="Dark2")
  
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="------", size= 6, color= "red") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="FIB2") +
  
  labs(title = "REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS", x="log2FIB2", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1226/REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_FIB2_log_slope.png", width = 2, height = 4)




#----------esc1----------------------

mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3_c <- log(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3+0.1,2)
mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3_slope <- slope_func(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3_c)
mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3_slope <- as.data.frame(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3_slope)




library(ggplot2)
library(reshape2)
require(gridExtra)

ggplot(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3_slope,aes(
  x="",
  y=mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3_slope
  #color=as.factor(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3.m$variable)
)) +
  
  geom_violin(fill = "grey80",colour = "#3366FF") +
  
  #scale_fill_brewer(palette="Dark2")
  
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="------", size= 6, color= "red") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="ESC1") +
  
  labs(title = "REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS", x="log2ESC1", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1226/REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_ESC1_log_slope.png", width = 2, height = 4)



#----------------esc2--------------------
mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4_c <- log(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4+0.1,2)
mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4_slope <- slope_func(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4_c)
mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4_slope <- as.data.frame(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4_slope)




library(ggplot2)
library(reshape2)
require(gridExtra)

ggplot(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4_slope,aes(
  x="",
  y=mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4_slope
  #color=as.factor(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4.m$variable)
)) +
  
  geom_violin(fill = "grey80",colour = "#3366FF") +
  
  #scale_fill_brewer(palette="Dark2")
  
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="------", size= 6, color= "red") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="ESC2") +
  
  labs(title = "REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS", x="log2ESC2", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1226/REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_ESC2_log_slope.png", width = 2, height = 4)


#PLOT FIB1+FIB2


#grid.arrange(FIB1, FIB2, ncol=2)





#test cor
mouse_DNA_repair_exp1 <- as.data.frame(t(mouse_DNA_repair_exp1)) 

mouse_DNA_repair_exp1$group <- c(1,2,3)

mouse_DNA_repair_exp1t <- t(cor(x=mouse_DNA_repair_exp1$group,
    y=mouse_DNA_repair_exp1[-137]))

mouse_DNA_repair_exp1[is.na(mouse_DNA_repair_exp1)] <- 0







#------new go gene sets----------------------




#------------------REGULATION_OF_DNA_DAMAGE_CHECKPOINT--------------------

REGULATION_OF_DNA_DAMAGE_CHECKPOINT <- read.table("../GO_data/GO_REGULATION_OF_DNA_DAMAGE_CHECKPOINT.txt",header = T, sep="\t")


#有一些NA值，证明一些human_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_gene 没有对应的老鼠基因
mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_gid <- human_gene$mouse_ensembl_gene_id[match(REGULATION_OF_DNA_DAMAGE_CHECKPOINT[2:nrow(REGULATION_OF_DNA_DAMAGE_CHECKPOINT),],human_gene$human_external_gene_name)]


mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_oid <- as.character(annotLookup2$original_id[match(mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_gid,annotLookup2$ensembl_gene_id)])



#载入fib1
temp=match(mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_oid,row.names(fib1))
temp=temp[!is.na(temp)]
mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_exp1 <- fib1[temp,]
#载入fib2
temp=match(mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_oid,row.names(fib2))
temp=temp[!is.na(temp)]
mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_exp2 <- fib2[temp,]
#import esc1
temp=match(mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_oid,row.names(esc1))
temp=temp[!is.na(temp)]
mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_exp3 <- esc1[temp,]
#载入esc2
temp=match(mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_oid,row.names(esc2))
temp=temp[!is.na(temp)]
mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_exp4 <- esc2[temp,]






#fib1
#除去0，都加0.1
#mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_exp1_slope <- mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_exp1+0.1

#mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_exp1_f2s <- log(mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_exp1$FIB_fast_TPM/mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_exp1$FIB_slow_TPM,2)


mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_exp1_c <- log(mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_exp1+0.1,2)
mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_exp1_slope <- slope_func(mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_exp1_c)
mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_exp1_slope <- as.data.frame(mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_exp1_slope)




library(ggplot2)
library(reshape2)
require(gridExtra)

ggplot(mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_exp1_slope,aes(
  x="",
  y=mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_exp1_slope
  #color=as.factor(mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_exp1.m$variable)
)) +
  
  geom_violin(fill = "grey80",colour = "#3366FF") +
  
  #scale_fill_brewer(palette="Dark2")
  
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="------", size= 6, color= "red") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="FIB1") +
  
  labs(title = "REGULATION_OF_DNA_DAMAGE_CHECKPOINT", x="log2FIB1", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1226/REGULATION_OF_DNA_DAMAGE_CHECKPOINT_fib1_log_slope.png", width = 2, height = 4)


#FIB2------------------------------------

mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_exp2_c <- log(mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_exp2+0.1,2)
mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_exp2_slope <- slope_func(mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_exp2_c)
mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_exp2_slope <- as.data.frame(mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_exp2_slope)




library(ggplot2)
library(reshape2)
require(gridExtra)

ggplot(mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_exp2_slope,aes(
  x="",
  y=mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_exp2_slope
  #color=as.factor(mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_exp2.m$variable)
)) +
  
  geom_violin(fill = "grey80",colour = "#3366FF") +
  
  #scale_fill_brewer(palette="Dark2")
  
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="------", size= 6, color= "red") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="FIB2") +
  
  labs(title = "REGULATION_OF_DNA_DAMAGE_CHECKPOINT", x="log2FIB2", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1226/REGULATION_OF_DNA_DAMAGE_CHECKPOINT_FIB2_log_slope.png", width = 2, height = 4)




#----------esc1----------------------

mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_exp3_c <- log(mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_exp3+0.1,2)
mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_exp3_slope <- slope_func(mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_exp3_c)
mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_exp3_slope <- as.data.frame(mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_exp3_slope)




library(ggplot2)
library(reshape2)
require(gridExtra)

ggplot(mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_exp3_slope,aes(
  x="",
  y=mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_exp3_slope
  #color=as.factor(mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_exp3.m$variable)
)) +
  
  geom_violin(fill = "grey80",colour = "#3366FF") +
  
  #scale_fill_brewer(palette="Dark2")
  
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="------", size= 6, color= "red") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="ESC1") +
  
  labs(title = "REGULATION_OF_DNA_DAMAGE_CHECKPOINT", x="log2ESC1", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1226/REGULATION_OF_DNA_DAMAGE_CHECKPOINT_ESC1_log_slope.png", width = 2, height = 4)



#----------------esc2--------------------
mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_exp4_c <- log(mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_exp4+0.1,2)
mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_exp4_slope <- slope_func(mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_exp4_c)
mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_exp4_slope <- as.data.frame(mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_exp4_slope)




library(ggplot2)
library(reshape2)
require(gridExtra)

ggplot(mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_exp4_slope,aes(
  x="",
  y=mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_exp4_slope
  #color=as.factor(mouse_REGULATION_OF_DNA_DAMAGE_CHECKPOINT_exp4.m$variable)
)) +
  
  geom_violin(fill = "grey80",colour = "#3366FF") +
  
  #scale_fill_brewer(palette="Dark2")
  
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="------", size= 6, color= "red") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="ESC2") +
  
  labs(title = "REGULATION_OF_DNA_DAMAGE_CHECKPOINT", x="log2ESC2", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1226/REGULATION_OF_DNA_DAMAGE_CHECKPOINT_ESC2_log_slope.png", width = 2, height = 4)





#------------------GO_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS








POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS <- read.table("../GO_data/GO_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS.txt",header = T, sep="\t")


#有一些NA值，证明一些human_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_gene 没有对应的老鼠基因
mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_gid <- human_gene$mouse_ensembl_gene_id[match(POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS[2:nrow(POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS),],human_gene$human_external_gene_name)]


mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_oid <- as.character(annotLookup2$original_id[match(mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_gid,annotLookup2$ensembl_gene_id)])



#载入fib1
temp=match(mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_oid,row.names(fib1))
temp=temp[!is.na(temp)]
mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1 <- fib1[temp,]
#载入fib2
temp=match(mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_oid,row.names(fib2))
temp=temp[!is.na(temp)]
mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2 <- fib2[temp,]
#import esc1
temp=match(mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_oid,row.names(esc1))
temp=temp[!is.na(temp)]
mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3 <- esc1[temp,]
#载入esc2
temp=match(mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_oid,row.names(esc2))
temp=temp[!is.na(temp)]
mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4 <- esc2[temp,]






#fib1
#除去0，都加0.1
#mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1_slope <- mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1+0.1

#mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1_f2s <- log(mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1$FIB_fast_TPM/mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1$FIB_slow_TPM,2)


mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1_c <- log(mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1+0.1,2)
mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1_slope <- slope_func(mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1_c)
mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1_slope <- as.data.frame(mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1_slope)




library(ggplot2)
library(reshape2)
require(gridExtra)

ggplot(mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1_slope,aes(
  x="",
  y=mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1_slope
  #color=as.factor(mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1.m$variable)
)) +
  
  geom_violin(fill = "grey80",colour = "#3366FF") +
  
  #scale_fill_brewer(palette="Dark2")
  
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="------", size= 6, color= "red") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="FIB1") +
  
  labs(title = "POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS", x="log2FIB1", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1226/POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_fib1_log_slope.png", width = 2, height = 4)


#FIB2------------------------------------

mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2_c <- log(mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2+0.1,2)
mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2_slope <- slope_func(mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2_c)
mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2_slope <- as.data.frame(mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2_slope)




library(ggplot2)
library(reshape2)
require(gridExtra)

ggplot(mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2_slope,aes(
  x="",
  y=mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2_slope
  #color=as.factor(mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2.m$variable)
)) +
  
  geom_violin(fill = "grey80",colour = "#3366FF") +
  
  #scale_fill_brewer(palette="Dark2")
  
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="------", size= 6, color= "red") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="FIB2") +
  
  labs(title = "POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS", x="log2FIB2", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1226/POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_FIB2_log_slope.png", width = 2, height = 4)




#----------esc1----------------------

mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3_c <- log(mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3+0.1,2)
mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3_slope <- slope_func(mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3_c)
mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3_slope <- as.data.frame(mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3_slope)




library(ggplot2)
library(reshape2)
require(gridExtra)

ggplot(mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3_slope,aes(
  x="",
  y=mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3_slope
  #color=as.factor(mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3.m$variable)
)) +
  
  geom_violin(fill = "grey80",colour = "#3366FF") +
  
  #scale_fill_brewer(palette="Dark2")
  
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="------", size= 6, color= "red") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="ESC1") +
  
  labs(title = "POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS", x="log2ESC1", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1226/POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_ESC1_log_slope.png", width = 2, height = 4)



#----------------esc2--------------------
mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4_c <- log(mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4+0.1,2)
mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4_slope <- slope_func(mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4_c)
mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4_slope <- as.data.frame(mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4_slope)




library(ggplot2)
library(reshape2)
require(gridExtra)

ggplot(mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4_slope,aes(
  x="",
  y=mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4_slope
  #color=as.factor(mouse_POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4.m$variable)
)) +
  
  geom_violin(fill = "grey80",colour = "#3366FF") +
  
  #scale_fill_brewer(palette="Dark2")
  
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="------", size= 6, color= "red") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="ESC2") +
  
  labs(title = "POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS", x="log2ESC2", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1226/POSITIVE_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_ESC2_log_slope.png", width = 2, height = 4)


#------------------WONG_MITOCHONDRIA_GENE_MODULE--------------------





WONG_MITOCHONDRIA_GENE_MODULE <- read.table("../GO_data/WONG_MITOCHONDRIA_GENE_MODULE.txt",header = T, sep="\t")


#有一些NA值，证明一些human_WONG_MITOCHONDRIA_GENE_MODULE_gene 没有对应的老鼠基因
mouse_WONG_MITOCHONDRIA_GENE_MODULE_gid <- human_gene$mouse_ensembl_gene_id[match(WONG_MITOCHONDRIA_GENE_MODULE[2:nrow(WONG_MITOCHONDRIA_GENE_MODULE),],human_gene$human_external_gene_name)]


mouse_WONG_MITOCHONDRIA_GENE_MODULE_oid <- as.character(annotLookup2$original_id[match(mouse_WONG_MITOCHONDRIA_GENE_MODULE_gid,annotLookup2$ensembl_gene_id)])



#载入fib1
temp=match(mouse_WONG_MITOCHONDRIA_GENE_MODULE_oid,row.names(fib1))
temp=temp[!is.na(temp)]
mouse_WONG_MITOCHONDRIA_GENE_MODULE_exp1 <- fib1[temp,]
#载入fib2
temp=match(mouse_WONG_MITOCHONDRIA_GENE_MODULE_oid,row.names(fib2))
temp=temp[!is.na(temp)]
mouse_WONG_MITOCHONDRIA_GENE_MODULE_exp2 <- fib2[temp,]
#import esc1
temp=match(mouse_WONG_MITOCHONDRIA_GENE_MODULE_oid,row.names(esc1))
temp=temp[!is.na(temp)]
mouse_WONG_MITOCHONDRIA_GENE_MODULE_exp3 <- esc1[temp,]
#载入esc2
temp=match(mouse_WONG_MITOCHONDRIA_GENE_MODULE_oid,row.names(esc2))
temp=temp[!is.na(temp)]
mouse_WONG_MITOCHONDRIA_GENE_MODULE_exp4 <- esc2[temp,]






#fib1
#除去0，都加0.1
#mouse_WONG_MITOCHONDRIA_GENE_MODULE_exp1_slope <- mouse_WONG_MITOCHONDRIA_GENE_MODULE_exp1+0.1

#mouse_WONG_MITOCHONDRIA_GENE_MODULE_exp1_f2s <- log(mouse_WONG_MITOCHONDRIA_GENE_MODULE_exp1$FIB_fast_TPM/mouse_WONG_MITOCHONDRIA_GENE_MODULE_exp1$FIB_slow_TPM,2)


mouse_WONG_MITOCHONDRIA_GENE_MODULE_exp1_c <- log(mouse_WONG_MITOCHONDRIA_GENE_MODULE_exp1+0.1,2)
mouse_WONG_MITOCHONDRIA_GENE_MODULE_exp1_slope <- slope_func(mouse_WONG_MITOCHONDRIA_GENE_MODULE_exp1_c)
mouse_WONG_MITOCHONDRIA_GENE_MODULE_exp1_slope <- as.data.frame(mouse_WONG_MITOCHONDRIA_GENE_MODULE_exp1_slope)




library(ggplot2)
library(reshape2)
require(gridExtra)

ggplot(mouse_WONG_MITOCHONDRIA_GENE_MODULE_exp1_slope,aes(
  x="",
  y=mouse_WONG_MITOCHONDRIA_GENE_MODULE_exp1_slope
  #color=as.factor(mouse_WONG_MITOCHONDRIA_GENE_MODULE_exp1.m$variable)
)) +
  
  geom_violin(fill = "grey80",colour = "#3366FF") +
  
  #scale_fill_brewer(palette="Dark2")
  
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="------", size= 6, color= "red") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="FIB1") +
  
  labs(title = "WONG_MITOCHONDRIA_GENE_MODULE", x="log2FIB1", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1226/WONG_MITOCHONDRIA_GENE_MODULE_fib1_log_slope.png", width = 2, height = 4)


#FIB2------------------------------------

mouse_WONG_MITOCHONDRIA_GENE_MODULE_exp2_c <- log(mouse_WONG_MITOCHONDRIA_GENE_MODULE_exp2+0.1,2)
mouse_WONG_MITOCHONDRIA_GENE_MODULE_exp2_slope <- slope_func(mouse_WONG_MITOCHONDRIA_GENE_MODULE_exp2_c)
mouse_WONG_MITOCHONDRIA_GENE_MODULE_exp2_slope <- as.data.frame(mouse_WONG_MITOCHONDRIA_GENE_MODULE_exp2_slope)




library(ggplot2)
library(reshape2)
require(gridExtra)

ggplot(mouse_WONG_MITOCHONDRIA_GENE_MODULE_exp2_slope,aes(
  x="",
  y=mouse_WONG_MITOCHONDRIA_GENE_MODULE_exp2_slope
  #color=as.factor(mouse_WONG_MITOCHONDRIA_GENE_MODULE_exp2.m$variable)
)) +
  
  geom_violin(fill = "grey80",colour = "#3366FF") +
  
  #scale_fill_brewer(palette="Dark2")
  
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="------", size= 6, color= "red") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="FIB2") +
  
  labs(title = "WONG_MITOCHONDRIA_GENE_MODULE", x="log2FIB2", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1226/WONG_MITOCHONDRIA_GENE_MODULE_FIB2_log_slope.png", width = 2, height = 4)




#----------esc1----------------------

mouse_WONG_MITOCHONDRIA_GENE_MODULE_exp3_c <- log(mouse_WONG_MITOCHONDRIA_GENE_MODULE_exp3+0.1,2)
mouse_WONG_MITOCHONDRIA_GENE_MODULE_exp3_slope <- slope_func(mouse_WONG_MITOCHONDRIA_GENE_MODULE_exp3_c)
mouse_WONG_MITOCHONDRIA_GENE_MODULE_exp3_slope <- as.data.frame(mouse_WONG_MITOCHONDRIA_GENE_MODULE_exp3_slope)




library(ggplot2)
library(reshape2)
require(gridExtra)

ggplot(mouse_WONG_MITOCHONDRIA_GENE_MODULE_exp3_slope,aes(
  x="",
  y=mouse_WONG_MITOCHONDRIA_GENE_MODULE_exp3_slope
  #color=as.factor(mouse_WONG_MITOCHONDRIA_GENE_MODULE_exp3.m$variable)
)) +
  
  geom_violin(fill = "grey80",colour = "#3366FF") +
  
  #scale_fill_brewer(palette="Dark2")
  
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="------", size= 6, color= "red") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="ESC1") +
  
  labs(title = "WONG_MITOCHONDRIA_GENE_MODULE", x="log2ESC1", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1226/WONG_MITOCHONDRIA_GENE_MODULE_ESC1_log_slope.png", width = 2, height = 4)



#----------------esc2--------------------
mouse_WONG_MITOCHONDRIA_GENE_MODULE_exp4_c <- log(mouse_WONG_MITOCHONDRIA_GENE_MODULE_exp4+0.1,2)
mouse_WONG_MITOCHONDRIA_GENE_MODULE_exp4_slope <- slope_func(mouse_WONG_MITOCHONDRIA_GENE_MODULE_exp4_c)
mouse_WONG_MITOCHONDRIA_GENE_MODULE_exp4_slope <- as.data.frame(mouse_WONG_MITOCHONDRIA_GENE_MODULE_exp4_slope)




library(ggplot2)
library(reshape2)
require(gridExtra)

ggplot(mouse_WONG_MITOCHONDRIA_GENE_MODULE_exp4_slope,aes(
  x="",
  y=mouse_WONG_MITOCHONDRIA_GENE_MODULE_exp4_slope
  #color=as.factor(mouse_WONG_MITOCHONDRIA_GENE_MODULE_exp4.m$variable)
)) +
  
  geom_violin(fill = "grey80",colour = "#3366FF") +
  
  #scale_fill_brewer(palette="Dark2")
  
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="------", size= 6, color= "red") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="ESC2") +
  
  labs(title = "WONG_MITOCHONDRIA_GENE_MODULE", x="log2ESC2", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1226/WONG_MITOCHONDRIA_GENE_MODULE_ESC2_log_slope.png", width = 2, height = 4)


#------------------REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION--------------------





REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION <- read.table("../GO_data/GO_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION.txt",header = T, sep="\t")


#有一些NA值，证明一些human_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_gene 没有对应的老鼠基因
mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_gid <- human_gene$mouse_ensembl_gene_id[match(REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION[2:nrow(REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION),],human_gene$human_external_gene_name)]


mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_oid <- as.character(annotLookup2$original_id[match(mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_gid,annotLookup2$ensembl_gene_id)])



#载入fib1
temp=match(mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_oid,row.names(fib1))
temp=temp[!is.na(temp)]
mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_exp1 <- fib1[temp,]
#载入fib2
temp=match(mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_oid,row.names(fib2))
temp=temp[!is.na(temp)]
mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_exp2 <- fib2[temp,]
#import esc1
temp=match(mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_oid,row.names(esc1))
temp=temp[!is.na(temp)]
mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_exp3 <- esc1[temp,]
#载入esc2
temp=match(mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_oid,row.names(esc2))
temp=temp[!is.na(temp)]
mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_exp4 <- esc2[temp,]






#fib1
#除去0，都加0.1
#mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_exp1_slope <- mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_exp1+0.1

#mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_exp1_f2s <- log(mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_exp1$FIB_fast_TPM/mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_exp1$FIB_slow_TPM,2)


mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_exp1_c <- log(mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_exp1+0.1,2)
mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_exp1_slope <- slope_func(mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_exp1_c)
mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_exp1_slope <- as.data.frame(mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_exp1_slope)




library(ggplot2)
library(reshape2)
require(gridExtra)

ggplot(mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_exp1_slope,aes(
  x="",
  y=mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_exp1_slope
  #color=as.factor(mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_exp1.m$variable)
)) +
  
  geom_violin(fill = "grey80",colour = "#3366FF") +
  
  #scale_fill_brewer(palette="Dark2")
  
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="------", size= 6, color= "red") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="FIB1") +
  
  labs(title = "REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION", x="log2FIB1", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1226/REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_fib1_log_slope.png", width = 2, height = 4)


#FIB2------------------------------------

mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_exp2_c <- log(mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_exp2+0.1,2)
mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_exp2_slope <- slope_func(mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_exp2_c)
mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_exp2_slope <- as.data.frame(mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_exp2_slope)




library(ggplot2)
library(reshape2)
require(gridExtra)

ggplot(mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_exp2_slope,aes(
  x="",
  y=mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_exp2_slope
  #color=as.factor(mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_exp2.m$variable)
)) +
  
  geom_violin(fill = "grey80",colour = "#3366FF") +
  
  #scale_fill_brewer(palette="Dark2")
  
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="------", size= 6, color= "red") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="FIB2") +
  
  labs(title = "REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION", x="log2FIB2", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1226/REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_FIB2_log_slope.png", width = 2, height = 4)




#----------esc1----------------------

mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_exp3_c <- log(mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_exp3+0.1,2)
mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_exp3_slope <- slope_func(mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_exp3_c)
mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_exp3_slope <- as.data.frame(mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_exp3_slope)




library(ggplot2)
library(reshape2)
require(gridExtra)

ggplot(mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_exp3_slope,aes(
  x="",
  y=mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_exp3_slope
  #color=as.factor(mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_exp3.m$variable)
)) +
  
  geom_violin(fill = "grey80",colour = "#3366FF") +
  
  #scale_fill_brewer(palette="Dark2")
  
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="------", size= 6, color= "red") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="ESC1") +
  
  labs(title = "REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION", x="log2ESC1", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1226/REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_ESC1_log_slope.png", width = 2, height = 4)



#----------------esc2--------------------
mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_exp4_c <- log(mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_exp4+0.1,2)
mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_exp4_slope <- slope_func(mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_exp4_c)
mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_exp4_slope <- as.data.frame(mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_exp4_slope)




library(ggplot2)
library(reshape2)
require(gridExtra)

ggplot(mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_exp4_slope,aes(
  x="",
  y=mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_exp4_slope
  #color=as.factor(mouse_REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_exp4.m$variable)
)) +
  
  geom_violin(fill = "grey80",colour = "#3366FF") +
  
  #scale_fill_brewer(palette="Dark2")
  
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="------", size= 6, color= "red") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="ESC2") +
  
  labs(title = "REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION", x="log2ESC2", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1226/REGULATION_OF_MITOCHONDRIAL_DEPOLARIZATION_ESC2_log_slope.png", width = 2, height = 4)


#------------------RIBOSOME_BIOGENESIS_AND_ASSEMBLY--------------------





RIBOSOME_BIOGENESIS_AND_ASSEMBLY <- read.table("../GO_data/GO_0042254RIBOSOME_BIOGENESIS_AND_ASSEMBLY.txt",header = T, sep="\t")


#有一些NA值，证明一些human_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_gene 没有对应的老鼠基因
mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_gid <- human_gene$mouse_ensembl_gene_id[match(RIBOSOME_BIOGENESIS_AND_ASSEMBLY[2:nrow(RIBOSOME_BIOGENESIS_AND_ASSEMBLY),],human_gene$human_external_gene_name)]


mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_oid <- as.character(annotLookup2$original_id[match(mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_gid,annotLookup2$ensembl_gene_id)])



#载入fib1
temp=match(mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_oid,row.names(fib1))
temp=temp[!is.na(temp)]
mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_exp1 <- fib1[temp,]
#载入fib2
temp=match(mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_oid,row.names(fib2))
temp=temp[!is.na(temp)]
mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_exp2 <- fib2[temp,]
#import esc1
temp=match(mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_oid,row.names(esc1))
temp=temp[!is.na(temp)]
mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_exp3 <- esc1[temp,]
#载入esc2
temp=match(mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_oid,row.names(esc2))
temp=temp[!is.na(temp)]
mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_exp4 <- esc2[temp,]






#fib1
#除去0，都加0.1
#mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_exp1_slope <- mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_exp1+0.1

#mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_exp1_f2s <- log(mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_exp1$FIB_fast_TPM/mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_exp1$FIB_slow_TPM,2)


mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_exp1_c <- log(mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_exp1+0.1,2)
mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_exp1_slope <- slope_func(mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_exp1_c)
mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_exp1_slope <- as.data.frame(mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_exp1_slope)




library(ggplot2)
library(reshape2)
require(gridExtra)

ggplot(mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_exp1_slope,aes(
  x="",
  y=mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_exp1_slope
  #color=as.factor(mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_exp1.m$variable)
)) +
  
  geom_violin(fill = "grey80",colour = "#3366FF") +
  
  #scale_fill_brewer(palette="Dark2")
  
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="------", size= 6, color= "red") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="FIB1") +
  
  labs(title = "RIBOSOME_BIOGENESIS_AND_ASSEMBLY", x="log2FIB1", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1226/RIBOSOME_BIOGENESIS_AND_ASSEMBLY_fib1_log_slope.png", width = 2, height = 4)


#FIB2------------------------------------

mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_exp2_c <- log(mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_exp2+0.1,2)
mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_exp2_slope <- slope_func(mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_exp2_c)
mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_exp2_slope <- as.data.frame(mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_exp2_slope)




library(ggplot2)
library(reshape2)
require(gridExtra)

ggplot(mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_exp2_slope,aes(
  x="",
  y=mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_exp2_slope
  #color=as.factor(mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_exp2.m$variable)
)) +
  
  geom_violin(fill = "grey80",colour = "#3366FF") +
  
  #scale_fill_brewer(palette="Dark2")
  
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="------", size= 6, color= "red") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="FIB2") +
  
  labs(title = "RIBOSOME_BIOGENESIS_AND_ASSEMBLY", x="log2FIB2", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1226/RIBOSOME_BIOGENESIS_AND_ASSEMBLY_FIB2_log_slope.png", width = 2, height = 4)




#----------esc1----------------------

mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_exp3_c <- log(mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_exp3+0.1,2)
mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_exp3_slope <- slope_func(mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_exp3_c)
mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_exp3_slope <- as.data.frame(mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_exp3_slope)




library(ggplot2)
library(reshape2)
require(gridExtra)

ggplot(mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_exp3_slope,aes(
  x="",
  y=mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_exp3_slope
  #color=as.factor(mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_exp3.m$variable)
)) +
  
  geom_violin(fill = "grey80",colour = "#3366FF") +
  
  #scale_fill_brewer(palette="Dark2")
  
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="------", size= 6, color= "red") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="ESC1") +
  
  labs(title = "RIBOSOME_BIOGENESIS_AND_ASSEMBLY", x="log2ESC1", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1226/RIBOSOME_BIOGENESIS_AND_ASSEMBLY_ESC1_log_slope.png", width = 2, height = 4)



#----------------esc2--------------------
mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_exp4_c <- log(mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_exp4+0.1,2)
mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_exp4_slope <- slope_func(mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_exp4_c)
mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_exp4_slope <- as.data.frame(mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_exp4_slope)




library(ggplot2)
library(reshape2)
require(gridExtra)

ggplot(mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_exp4_slope,aes(
  x="",
  y=mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_exp4_slope
  #color=as.factor(mouse_RIBOSOME_BIOGENESIS_AND_ASSEMBLY_exp4.m$variable)
)) +
  
  geom_violin(fill = "grey80",colour = "#3366FF") +
  
  #scale_fill_brewer(palette="Dark2")
  
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="------", size= 6, color= "red") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="ESC2") +
  
  labs(title = "RIBOSOME_BIOGENESIS_AND_ASSEMBLY", x="log2ESC2", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1226/RIBOSOME_BIOGENESIS_AND_ASSEMBLY_ESC2_log_slope.png", width = 2, height = 4)


#------------------HALLMARK_PEROXISOME--------------------





HALLMARK_PEROXISOME <- read.table("../GO_data/HALLMARK_PEROXISOME.txt",header = T, sep="\t")


#有一些NA值，证明一些human_HALLMARK_PEROXISOME_gene 没有对应的老鼠基因
mouse_HALLMARK_PEROXISOME_gid <- human_gene$mouse_ensembl_gene_id[match(HALLMARK_PEROXISOME[2:nrow(HALLMARK_PEROXISOME),],human_gene$human_external_gene_name)]


mouse_HALLMARK_PEROXISOME_oid <- as.character(annotLookup2$original_id[match(mouse_HALLMARK_PEROXISOME_gid,annotLookup2$ensembl_gene_id)])



#载入fib1
temp=match(mouse_HALLMARK_PEROXISOME_oid,row.names(fib1))
temp=temp[!is.na(temp)]
mouse_HALLMARK_PEROXISOME_exp1 <- fib1[temp,]
#载入fib2
temp=match(mouse_HALLMARK_PEROXISOME_oid,row.names(fib2))
temp=temp[!is.na(temp)]
mouse_HALLMARK_PEROXISOME_exp2 <- fib2[temp,]
#import esc1
temp=match(mouse_HALLMARK_PEROXISOME_oid,row.names(esc1))
temp=temp[!is.na(temp)]
mouse_HALLMARK_PEROXISOME_exp3 <- esc1[temp,]
#载入esc2
temp=match(mouse_HALLMARK_PEROXISOME_oid,row.names(esc2))
temp=temp[!is.na(temp)]
mouse_HALLMARK_PEROXISOME_exp4 <- esc2[temp,]






#fib1
#除去0，都加0.1
#mouse_HALLMARK_PEROXISOME_exp1_slope <- mouse_HALLMARK_PEROXISOME_exp1+0.1

#mouse_HALLMARK_PEROXISOME_exp1_f2s <- log(mouse_HALLMARK_PEROXISOME_exp1$FIB_fast_TPM/mouse_HALLMARK_PEROXISOME_exp1$FIB_slow_TPM,2)


mouse_HALLMARK_PEROXISOME_exp1_c <- log(mouse_HALLMARK_PEROXISOME_exp1+0.1,2)
mouse_HALLMARK_PEROXISOME_exp1_slope <- slope_func(mouse_HALLMARK_PEROXISOME_exp1_c)
mouse_HALLMARK_PEROXISOME_exp1_slope <- as.data.frame(mouse_HALLMARK_PEROXISOME_exp1_slope)




library(ggplot2)
library(reshape2)
require(gridExtra)

ggplot(mouse_HALLMARK_PEROXISOME_exp1_slope,aes(
  x="",
  y=mouse_HALLMARK_PEROXISOME_exp1_slope
  #color=as.factor(mouse_HALLMARK_PEROXISOME_exp1.m$variable)
)) +
  
  geom_violin(fill = "grey80",colour = "#3366FF") +
  
  #scale_fill_brewer(palette="Dark2")
  
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="------", size= 6, color= "red") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="FIB1") +
  
  labs(title = "HALLMARK_PEROXISOME", x="log2FIB1", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1226/HALLMARK_PEROXISOME_fib1_log_slope.png", width = 2, height = 4)


#FIB2------------------------------------

mouse_HALLMARK_PEROXISOME_exp2_c <- log(mouse_HALLMARK_PEROXISOME_exp2+0.1,2)
mouse_HALLMARK_PEROXISOME_exp2_slope <- slope_func(mouse_HALLMARK_PEROXISOME_exp2_c)
mouse_HALLMARK_PEROXISOME_exp2_slope <- as.data.frame(mouse_HALLMARK_PEROXISOME_exp2_slope)




library(ggplot2)
library(reshape2)
require(gridExtra)

ggplot(mouse_HALLMARK_PEROXISOME_exp2_slope,aes(
  x="",
  y=mouse_HALLMARK_PEROXISOME_exp2_slope
  #color=as.factor(mouse_HALLMARK_PEROXISOME_exp2.m$variable)
)) +
  
  geom_violin(fill = "grey80",colour = "#3366FF") +
  
  #scale_fill_brewer(palette="Dark2")
  
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="------", size= 6, color= "red") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="FIB2") +
  
  labs(title = "HALLMARK_PEROXISOME", x="log2FIB2", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1226/HALLMARK_PEROXISOME_FIB2_log_slope.png", width = 2, height = 4)




#----------esc1----------------------

mouse_HALLMARK_PEROXISOME_exp3_c <- log(mouse_HALLMARK_PEROXISOME_exp3+0.1,2)
mouse_HALLMARK_PEROXISOME_exp3_slope <- slope_func(mouse_HALLMARK_PEROXISOME_exp3_c)
mouse_HALLMARK_PEROXISOME_exp3_slope <- as.data.frame(mouse_HALLMARK_PEROXISOME_exp3_slope)




library(ggplot2)
library(reshape2)
require(gridExtra)

ggplot(mouse_HALLMARK_PEROXISOME_exp3_slope,aes(
  x="",
  y=mouse_HALLMARK_PEROXISOME_exp3_slope
  #color=as.factor(mouse_HALLMARK_PEROXISOME_exp3.m$variable)
)) +
  
  geom_violin(fill = "grey80",colour = "#3366FF") +
  
  #scale_fill_brewer(palette="Dark2")
  
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="------", size= 6, color= "red") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="ESC1") +
  
  labs(title = "HALLMARK_PEROXISOME", x="log2ESC1", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1226/HALLMARK_PEROXISOME_ESC1_log_slope.png", width = 2, height = 4)



#----------------esc2--------------------
mouse_HALLMARK_PEROXISOME_exp4_c <- log(mouse_HALLMARK_PEROXISOME_exp4+0.1,2)
mouse_HALLMARK_PEROXISOME_exp4_slope <- slope_func(mouse_HALLMARK_PEROXISOME_exp4_c)
mouse_HALLMARK_PEROXISOME_exp4_slope <- as.data.frame(mouse_HALLMARK_PEROXISOME_exp4_slope)




library(ggplot2)
library(reshape2)
require(gridExtra)

ggplot(mouse_HALLMARK_PEROXISOME_exp4_slope,aes(
  x="",
  y=mouse_HALLMARK_PEROXISOME_exp4_slope
  #color=as.factor(mouse_HALLMARK_PEROXISOME_exp4.m$variable)
)) +
  
  geom_violin(fill = "grey80",colour = "#3366FF") +
  
  #scale_fill_brewer(palette="Dark2")
  
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="------", size= 6, color= "red") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="ESC2") +
  
  labs(title = "HALLMARK_PEROXISOME", x="log2ESC2", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1226/HALLMARK_PEROXISOME_ESC2_log_slope.png", width = 2, height = 4)



#------------------HALLMARK_P53_PATHWAY--------------------





HALLMARK_P53_PATHWAY <- read.table("../GO_data/HALLMARK_P53_PATHWAY.txt",header = T, sep="\t")


#有一些NA值，证明一些human_HALLMARK_P53_PATHWAY_gene 没有对应的老鼠基因
mouse_HALLMARK_P53_PATHWAY_gid <- human_gene$mouse_ensembl_gene_id[match(HALLMARK_P53_PATHWAY[2:nrow(HALLMARK_P53_PATHWAY),],human_gene$human_external_gene_name)]


mouse_HALLMARK_P53_PATHWAY_oid <- as.character(annotLookup2$original_id[match(mouse_HALLMARK_P53_PATHWAY_gid,annotLookup2$ensembl_gene_id)])



#载入fib1
temp=match(mouse_HALLMARK_P53_PATHWAY_oid,row.names(fib1))
temp=temp[!is.na(temp)]
mouse_HALLMARK_P53_PATHWAY_exp1 <- fib1[temp,]
#载入fib2
temp=match(mouse_HALLMARK_P53_PATHWAY_oid,row.names(fib2))
temp=temp[!is.na(temp)]
mouse_HALLMARK_P53_PATHWAY_exp2 <- fib2[temp,]
#import esc1
temp=match(mouse_HALLMARK_P53_PATHWAY_oid,row.names(esc1))
temp=temp[!is.na(temp)]
mouse_HALLMARK_P53_PATHWAY_exp3 <- esc1[temp,]
#载入esc2
temp=match(mouse_HALLMARK_P53_PATHWAY_oid,row.names(esc2))
temp=temp[!is.na(temp)]
mouse_HALLMARK_P53_PATHWAY_exp4 <- esc2[temp,]






#fib1
#除去0，都加0.1
#mouse_HALLMARK_P53_PATHWAY_exp1_slope <- mouse_HALLMARK_P53_PATHWAY_exp1+0.1

#mouse_HALLMARK_P53_PATHWAY_exp1_f2s <- log(mouse_HALLMARK_P53_PATHWAY_exp1$FIB_fast_TPM/mouse_HALLMARK_P53_PATHWAY_exp1$FIB_slow_TPM,2)


mouse_HALLMARK_P53_PATHWAY_exp1_c <- log(mouse_HALLMARK_P53_PATHWAY_exp1+0.1,2)
mouse_HALLMARK_P53_PATHWAY_exp1_slope <- slope_func(mouse_HALLMARK_P53_PATHWAY_exp1_c)
mouse_HALLMARK_P53_PATHWAY_exp1_slope <- as.data.frame(mouse_HALLMARK_P53_PATHWAY_exp1_slope)




library(ggplot2)
library(reshape2)
require(gridExtra)

ggplot(mouse_HALLMARK_P53_PATHWAY_exp1_slope,aes(
  x="",
  y=mouse_HALLMARK_P53_PATHWAY_exp1_slope
  #color=as.factor(mouse_HALLMARK_P53_PATHWAY_exp1.m$variable)
)) +
  
  geom_violin(fill = "grey80",colour = "#3366FF") +
  
  #scale_fill_brewer(palette="Dark2")
  
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="------", size= 6, color= "red") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="FIB1") +
  
  labs(title = "HALLMARK_P53_PATHWAY", x="log2FIB1", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1226/HALLMARK_P53_PATHWAY_fib1_log_slope.png", width = 2, height = 4)


#FIB2------------------------------------

mouse_HALLMARK_P53_PATHWAY_exp2_c <- log(mouse_HALLMARK_P53_PATHWAY_exp2+0.1,2)
mouse_HALLMARK_P53_PATHWAY_exp2_slope <- slope_func(mouse_HALLMARK_P53_PATHWAY_exp2_c)
mouse_HALLMARK_P53_PATHWAY_exp2_slope <- as.data.frame(mouse_HALLMARK_P53_PATHWAY_exp2_slope)




library(ggplot2)
library(reshape2)
require(gridExtra)

ggplot(mouse_HALLMARK_P53_PATHWAY_exp2_slope,aes(
  x="",
  y=mouse_HALLMARK_P53_PATHWAY_exp2_slope
  #color=as.factor(mouse_HALLMARK_P53_PATHWAY_exp2.m$variable)
)) +
  
  geom_violin(fill = "grey80",colour = "#3366FF") +
  
  #scale_fill_brewer(palette="Dark2")
  
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="------", size= 6, color= "red") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="FIB2") +
  
  labs(title = "HALLMARK_P53_PATHWAY", x="log2FIB2", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1226/HALLMARK_P53_PATHWAY_FIB2_log_slope.png", width = 2, height = 4)




#----------esc1----------------------

mouse_HALLMARK_P53_PATHWAY_exp3_c <- log(mouse_HALLMARK_P53_PATHWAY_exp3+0.1,2)
mouse_HALLMARK_P53_PATHWAY_exp3_slope <- slope_func(mouse_HALLMARK_P53_PATHWAY_exp3_c)
mouse_HALLMARK_P53_PATHWAY_exp3_slope <- as.data.frame(mouse_HALLMARK_P53_PATHWAY_exp3_slope)




library(ggplot2)
library(reshape2)
require(gridExtra)

ggplot(mouse_HALLMARK_P53_PATHWAY_exp3_slope,aes(
  x="",
  y=mouse_HALLMARK_P53_PATHWAY_exp3_slope
  #color=as.factor(mouse_HALLMARK_P53_PATHWAY_exp3.m$variable)
)) +
  
  geom_violin(fill = "grey80",colour = "#3366FF") +
  
  #scale_fill_brewer(palette="Dark2")
  
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="------", size= 6, color= "red") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="ESC1") +
  
  labs(title = "HALLMARK_P53_PATHWAY", x="log2ESC1", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1226/HALLMARK_P53_PATHWAY_ESC1_log_slope.png", width = 2, height = 4)



#----------------esc2--------------------
mouse_HALLMARK_P53_PATHWAY_exp4_c <- log(mouse_HALLMARK_P53_PATHWAY_exp4+0.1,2)
mouse_HALLMARK_P53_PATHWAY_exp4_slope <- slope_func(mouse_HALLMARK_P53_PATHWAY_exp4_c)
mouse_HALLMARK_P53_PATHWAY_exp4_slope <- as.data.frame(mouse_HALLMARK_P53_PATHWAY_exp4_slope)




library(ggplot2)
library(reshape2)
require(gridExtra)

ggplot(mouse_HALLMARK_P53_PATHWAY_exp4_slope,aes(
  x="",
  y=mouse_HALLMARK_P53_PATHWAY_exp4_slope
  #color=as.factor(mouse_HALLMARK_P53_PATHWAY_exp4.m$variable)
)) +
  
  geom_violin(fill = "grey80",colour = "#3366FF") +
  
  #scale_fill_brewer(palette="Dark2")
  
  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
  
  #geom_jitter(height = 0, width = 0.1) +
  
  geom_boxplot(width=0.3) +
  stat_summary(fun.y = "mean", geom = "text", label="------", size= 6, color= "red") +
  #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.4)+
  
  #stat_summary(geom = "errorbar", fun.data = mean_se, position = "dodge")
  
  scale_colour_discrete(name  ="ESC2") +
  
  labs(title = "HALLMARK_P53_PATHWAY", x="log2ESC2", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1226/HALLMARK_P53_PATHWAY_ESC2_log_slope.png", width = 2, height = 4)