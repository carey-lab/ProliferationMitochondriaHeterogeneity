
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




mouse_DNA_repair_exp1 <- mouse_DNA_repair_exp1+0.1

mouse_DNA_repair_exp1_f2s <- log(mouse_DNA_repair_exp1$FIB_fast_TPM/mouse_DNA_repair_exp1$FIB_slow_TPM,2)

mouse_DNA_repair_exp1_f2s <- as.data.frame(mouse_DNA_repair_exp1_f2s)
mouse_DNA_repair_exp1_slope <- slope_func(mouse_DNA_repair_exp1)
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
  
  labs(title = "DNA_repair", x="FIB1", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1225/DNA_repair_fib1_slope.png", width = 2, height = 4)

#FIB2------------------------------------

mouse_DNA_repair_exp2_slope <- slope_func(mouse_DNA_repair_exp2)
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
  
  labs(title = "DNA_repair", x="FIB2", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1225/DNA_repair_FIB2_slope.png", width = 2, height = 4)


#----------esc1----------------------

mouse_DNA_repair_exp3_slope <- slope_func(mouse_DNA_repair_exp3)
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
  
  labs(title = "DNA_repair", x="ESC1", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1225/DNA_repair_ESC1_slope.png", width = 2, height = 4)

#----------------esc2--------------------
mouse_DNA_repair_exp4_slope <- slope_func(mouse_DNA_repair_exp4)
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
  
  labs(title = "DNA_repair", x="ESC2", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1225/DNA_repair_ESC2_slope.png", width = 2, height = 4)



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


mouse_GO_RIBOSOME_exp1_slope <- slope_func(mouse_GO_RIBOSOME_exp1)
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
  
  labs(title = "GO_RIBOSOME", x="FIB1", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1225/GO_RIBOSOME_fib1_slope.png", width = 2, height = 4)

#FIB2------------------------------------

mouse_GO_RIBOSOME_exp2_slope <- slope_func(mouse_GO_RIBOSOME_exp2)
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
  
  labs(title = "GO_RIBOSOME", x="FIB2", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1225/GO_RIBOSOME_FIB2_slope.png", width = 2, height = 4)


#----------esc1----------------------

mouse_GO_RIBOSOME_exp3_slope <- slope_func(mouse_GO_RIBOSOME_exp3)
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
  
  labs(title = "GO_RIBOSOME", x="ESC1", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1225/GO_RIBOSOME_ESC1_slope.png", width = 2, height = 4)

#----------------esc2--------------------
mouse_GO_RIBOSOME_exp4_slope <- slope_func(mouse_GO_RIBOSOME_exp4)
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
  
  labs(title = "GO_RIBOSOME", x="ESC2", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1225/GO_RIBOSOME_ESC2_slope.png", width = 2, height = 4)



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


mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1_slope <- slope_func(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1)
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
  
  labs(title = "RESPONSE_TO_DNA_DAMAGE_STIMULUS", x="FIB1", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1225/RESPONSE_TO_DNA_DAMAGE_STIMULUS_fib1_slope.png", width = 2, height = 4)

#FIB2------------------------------------

mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2_slope <- slope_func(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2)
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
  
  labs(title = "RESPONSE_TO_DNA_DAMAGE_STIMULUS", x="FIB2", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1225/RESPONSE_TO_DNA_DAMAGE_STIMULUS_FIB2_slope.png", width = 2, height = 4)


#----------esc1----------------------

mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3_slope <- slope_func(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3)
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
  
  labs(title = "RESPONSE_TO_DNA_DAMAGE_STIMULUS", x="ESC1", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1225/RESPONSE_TO_DNA_DAMAGE_STIMULUS_ESC1_slope.png", width = 2, height = 4)

#----------------esc2--------------------
mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4_slope <- slope_func(mouse_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4)
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
  
  labs(title = "RESPONSE_TO_DNA_DAMAGE_STIMULUS", x="ESC2", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1225/RESPONSE_TO_DNA_DAMAGE_STIMULUS_ESC2_slope.png", width = 2, height = 4)

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


mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1_slope <- slope_func(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp1)
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
  
  labs(title = "REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS", x="FIB1", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1225/REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_fib1_slope.png", width = 2, height = 4)

#FIB2------------------------------------

mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2_slope <- slope_func(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp2)
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
  
  labs(title = "REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS", x="FIB2", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1225/REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_FIB2_slope.png", width = 2, height = 4)


#----------esc1----------------------

mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3_slope <- slope_func(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp3)
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
  
  labs(title = "REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS", x="ESC1", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1225/REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_ESC1_slope.png", width = 2, height = 4)

#----------------esc2--------------------
mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4_slope <- slope_func(mouse_REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_exp4)
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
  
  labs(title = "REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS", x="ESC2", y="slope") +
  
  #coord_fixed(0.2) +
  
  geom_hline(yintercept=0, linetype="dashed", 
             color = "blue", size=0.5) +
  
  theme(
    plot.title = element_text(lineheight=.8, size=15,face="bold",hjust = 0.5))

ggsave("./figures/1225/REGULATION_OF_RESPONSE_TO_DNA_DAMAGE_STIMULUS_ESC2_slope.png", width = 2, height = 4)


#PLOT FIB1+FIB2


#grid.arrange(FIB1, FIB2, ncol=2)





























#test cor
mouse_DNA_repair_exp1 <- as.data.frame(t(mouse_DNA_repair_exp1)) 

mouse_DNA_repair_exp1$group <- c(1,2,3)

mouse_DNA_repair_exp1t <- t(cor(x=mouse_DNA_repair_exp1$group,
    y=mouse_DNA_repair_exp1[-137]))

mouse_DNA_repair_exp1[is.na(mouse_DNA_repair_exp1)] <- 0



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


