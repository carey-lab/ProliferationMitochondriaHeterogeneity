# Jan 17

# get the yeast GSEA result-----

yeast_fast <- read.table("./G1_M and Fast_Slow/result from server0116/yeast_0116_gsea_report_for_fast_1547608544528.tsv",header = T,sep = "\t")
head(yeast_fast)

yeast_slow <- read.table("./G1_M and Fast_Slow/result from server0116/yeast_0116_gsea_report_for_slow_1547608544528.tsv",header = T,sep = "\t")

yeast <- rbind(yeast_slow,yeast_fast)
nrow(yeast)
head(yeast)
yeast <- yeast[,c(1,4,5,6,7,8)]


# create fibroblasts, ESCs and yeast table

# read fibroblasts and ESCs
fib_esc <- read.table("./NES_all.tsv",header = T,sep = "\t")
head(fib_esc)

# merge
yeast <- yeast[match(fib_esc$NAME.FIB,yeast$NAME),]

colnames(yeast) <- c("NAME.yeast","SIZE_yeast","ES_yeast","NES_yeast","NOM.p.val","FDR.q.val")

NES_FIB_ESCs_yeast <- cbind(fib_esc,yeast)




write.table(NES_FIB_ESCs_yeast,file = "NES_FIB_ESCs_yeast.tsv",col.names = T,row.names = F,sep = "\t",quote = F)

# all yeast GSEA vs. FIB------

NES_FIB_ESCs_yeast <- read.table("./NES_FIB_ESCs_yeast.tsv",header = T,sep = "\t")
head(NES_FIB_ESCs_yeast)
colnames(NES_FIB_ESCs_yeast)
yeast_fib <- NES_FIB_ESCs_yeast[,c(1,4,17)]
head(yeast_fib)

nrow(yeast_fib)
yeast_fib <- na.omit(yeast_fib)

yeast_fib <- yeast_fib[! (yeast_fib$NES_FIB == 0 | yeast_fib$NES_yeast ==0),]
#plot-------------------

library(grid)
library(gridExtra)
grob_p = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(yeast_fib$NES_yeast, yeast_fib$NES_FIB,method = "pearson"), 4) ), x = 0.5, y = 0.97, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))
grob_s = grobTree(textGrob(paste("Spearman Correlation : ", round(cor(yeast_fib$NES_yeast, yeast_fib$NES_FIB,method = "spearman"), 4) ), x = 0.5, y = 0.93, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))



library(ggplot2)
ggplot(yeast_fib, aes(NES_yeast, NES_FIB))+ 
  geom_point() + 
  
  #stat_summary(fun.data=mean_cl_normal) + 
  #add a regression line
  geom_smooth(method='lm',formula=y~x,se=TRUE, fullrange=FALSE, level=0.95) +
  #stat_cor(method = "pearson", label.x = 0, label.y = 0) +
  
  annotation_custom(grob_p) +
  
  annotation_custom(grob_s) +
  
  
  
  labs(title = "yeast vs. fibroblasts", x="NES_yeast", y="NES_fib") +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", size=1) +
  geom_vline(xintercept=0, linetype="dashed", 
             color = "black", size=1) +
  
  #geom_point(size = 2) +
  
  #geom_line(size=1) +
  
  #labs(title = "Tars") +
  
  # #scale_colour_discrete(name  ="Cell Type",
  #                       breaks=c("1", "2"),
  #                       labels=c("ESC", "FIB")) +
  #去除图例标题
#guides(fill=guide_legend(title="Cell Type")) +
#主题
theme(
  plot.title = element_text(lineheight=.8, size=20,face="bold",hjust = 0.5),
  
  axis.title.x = element_text(color="blue", size=14, face="bold"),
  axis.title.y = element_text(color="blue", size=14, face="bold"),
  #delete background
  panel.grid.major =element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  
  #加上坐标轴
  axis.line = element_line(colour = "black",size=1.0),
  #刻度线
  axis.ticks = element_line(size=0.5),
  axis.ticks.length=unit(0.15,"cm"),
  
  #x轴标签
  axis.text.x = element_text(size=14,color='black'),
  #y轴标签
  axis.text.y = element_text(size=14,color='black'),
  #图例
  legend.title = element_text(colour="black", size=14, face="bold")
)



# all ESCs GSEA vs. yeast------
colnames(NES_FIB_ESCs_yeast)
ESCs_yeast <- NES_FIB_ESCs_yeast[,c(1,9,17)]
head(ESCs_yeast)

nrow(ESCs_yeast)
ESCs_yeast <- na.omit(ESCs_yeast)

ESCs_yeast <- ESCs_yeast[! (ESCs_yeast$NES_ESC == 0 | ESCs_yeast$NES_yeast ==0),]

#plot-------------------

library(grid)
library(gridExtra)
grob_p = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(ESCs_yeast$NES_ESC, ESCs_yeast$NES_yeast,method = "pearson"), 4) ), x = 0.5, y = 0.97, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))
grob_s = grobTree(textGrob(paste("Spearman Correlation : ", round(cor(ESCs_yeast$NES_ESC, ESCs_yeast$NES_yeast,method = "spearman"), 4) ), x = 0.5, y = 0.93, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))



library(ggplot2)
ggplot(ESCs_yeast, aes(NES_ESC, NES_yeast))+ 
  geom_point() + 
  
  #stat_summary(fun.data=mean_cl_normal) + 
  #add a regression line
  geom_smooth(method='lm',formula=y~x,se=TRUE, fullrange=FALSE, level=0.95) +
  #stat_cor(method = "pearson", label.x = 0, label.y = 0) +
  
  annotation_custom(grob_p) +
  
  annotation_custom(grob_s) +
  
  
  
  labs(title = "ESCs vs. yeast", x="NES_ESC", y="NES_yeast") +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", size=1) +
  geom_vline(xintercept=0, linetype="dashed", 
             color = "black", size=1) +
  
  #geom_point(size = 2) +
  
  #geom_line(size=1) +
  
  #labs(title = "Tars") +
  
  # #scale_colour_discrete(name  ="Cell Type",
  #                       breaks=c("1", "2"),
  #                       labels=c("ESC", "FIB")) +
  #去除图例标题
#guides(fill=guide_legend(title="Cell Type")) +
#主题
theme(
  plot.title = element_text(lineheight=.8, size=20,face="bold",hjust = 0.5),
  axis.title.x = element_text(color="blue", size=14, face="bold"),
  axis.title.y = element_text(color="blue", size=14, face="bold"),
  #delete background
  panel.grid.major =element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  
  #加上坐标轴
  axis.line = element_line(colour = "black",size=1.0),
  #刻度线
  axis.ticks = element_line(size=0.5),
  axis.ticks.length=unit(0.15,"cm"),
  
  #x轴标签
  axis.text.x = element_text(size=14,color='black'),
  #y轴标签
  axis.text.y = element_text(size=14,color='black'),
  #图例
  legend.title = element_text(colour="black", size=14, face="bold")
)



# all FIB GSEA vs. ESCs------
colnames(NES_FIB_ESCs_yeast)
fib_yeast <- NES_FIB_ESCs_yeast[,c(1,9,17)]
head(fib_yeast)

nrow(fib_yeast)
fib_yeast <- na.omit(fib_yeast)



#plot-------------------

library(grid)
library(gridExtra)
grob_p = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(fib_yeast$NES_ESC, fib_yeast$NES_yeast,method = "pearson"), 4) ), x = 0.5, y = 0.97, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))
grob_s = grobTree(textGrob(paste("Spearman Correlation : ", round(cor(fib_yeast$NES_ESC, fib_yeast$NES_yeast,method = "spearman"), 4) ), x = 0.5, y = 0.93, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))



library(ggplot2)
ggplot(fib_yeast, aes(NES_ESC, NES_yeast))+ 
  geom_point() + 
  
  #stat_summary(fun.data=mean_cl_normal) + 
  #add a regression line
  geom_smooth(method='lm',formula=y~x,se=TRUE, fullrange=FALSE, level=0.95) +
  #stat_cor(method = "pearson", label.x = 0, label.y = 0) +
  
  annotation_custom(grob_p) +
  
  annotation_custom(grob_s) +
  
  
  
  labs(title = "ESCs vs. yeast", x="NES_ESC", y="NES_yeast") +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", size=1) +
  geom_vline(xintercept=0, linetype="dashed", 
             color = "black", size=1) +
  
  #geom_point(size = 2) +
  
  #geom_line(size=1) +
  
  #labs(title = "Tars") +
  
  # #scale_colour_discrete(name  ="Cell Type",
  #                       breaks=c("1", "2"),
  #                       labels=c("ESC", "FIB")) +
  #去除图例标题
#guides(fill=guide_legend(title="Cell Type")) +
#主题
theme(
  plot.title = element_text(lineheight=.8, size=20,face="bold",hjust = 0.5),
  axis.title.x = element_text(color="blue", size=14, face="bold"),
  axis.title.y = element_text(color="blue", size=14, face="bold"),
  #delete background
  panel.grid.major =element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  
  #加上坐标轴
  axis.line = element_line(colour = "black",size=1.0),
  #刻度线
  axis.ticks = element_line(size=0.5),
  axis.ticks.length=unit(0.15,"cm"),
  
  #x轴标签
  axis.text.x = element_text(size=14,color='black'),
  #y轴标签
  axis.text.y = element_text(size=14,color='black'),
  #图例
  legend.title = element_text(colour="black", size=14, face="bold")
)



# all FIB GSEA vs. ESCs------
colnames(NES_FIB_ESCs_yeast)
fib_yeast <- NES_FIB_ESCs_yeast[,c(1,4,9)]
head(fib_yeast)

nrow(fib_yeast)
fib_yeast <- na.omit(fib_yeast)

# row_sub = apply(fib_yeast, 1, function(row) all(row !=0 ))
# ##Subset as usual
# fib_yeast <- fib_yeast[row_sub,]

fib_yeast <- fib_yeast[! (fib_yeast$NES_FIB == 0 | fib_yeast$NES_ESC ==0),]

#plot-------------------

library(grid)
library(gridExtra)
grob_p = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(fib_yeast$NES_FIB, fib_yeast$NES_ESC,method = "pearson"), 4) ), x = 0.5, y = 0.97, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))
grob_s = grobTree(textGrob(paste("Spearman Correlation : ", round(cor(fib_yeast$NES_FIB, fib_yeast$NES_ESC,method = "spearman"), 4) ), x = 0.5, y = 0.93, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))



library(ggplot2)
ggplot(fib_yeast, aes(NES_FIB, NES_ESC))+ 
  geom_point() + 
  
  #stat_summary(fun.data=mean_cl_normal) + 
  #add a regression line
  geom_smooth(method='lm',formula=y~x,se=TRUE, fullrange=FALSE, level=0.95) +
  #stat_cor(method = "pearson", label.x = 0, label.y = 0) +
  
  annotation_custom(grob_p) +
  
  annotation_custom(grob_s) +
  
  
  
  labs(title = "fibroblasts vs. ESCs", x="NES_FIB", y="NES_ESC") +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", size=1) +
  geom_vline(xintercept=0, linetype="dashed", 
             color = "black", size=1) +
  
  #geom_point(size = 2) +
  
  #geom_line(size=1) +
  
  #labs(title = "Tars") +
  
  # #scale_colour_discrete(name  ="Cell Type",
  #                       breaks=c("1", "2"),
  #                       labels=c("ESC", "FIB")) +
  #去除图例标题
#guides(fill=guide_legend(title="Cell Type")) +
#主题
theme(
  plot.title = element_text(lineheight=.8, size=20,face="bold",hjust = 0.5),
  axis.title.x = element_text(color="blue", size=14, face="bold"),
  axis.title.y = element_text(color="blue", size=14, face="bold"),
  #delete background
  panel.grid.major =element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  
  #加上坐标轴
  axis.line = element_line(colour = "black",size=1.0),
  #刻度线
  axis.ticks = element_line(size=0.5),
  axis.ticks.length=unit(0.15,"cm"),
  
  #x轴标签
  axis.text.x = element_text(size=14,color='black'),
  #y轴标签
  axis.text.y = element_text(size=14,color='black'),
  #图例
  legend.title = element_text(colour="black", size=14, face="bold")
)


















