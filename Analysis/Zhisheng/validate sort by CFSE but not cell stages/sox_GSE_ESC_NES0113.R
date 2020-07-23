# Jan 14 2019
# plot sox2 vs. ESC, GSE vs. ESC, sox2 vs. GSE


########################

# load ESC rb ng sn
esc2016_2018_fast_rb_sn_ng <- read.table("./GSEA/result_from_server2/esc2016+2018_rb_negative_sn_gsea_report_for_FAST_1546782445406.tsv", header = T,sep="\t")
esc2016_2018_slow_rb_sn_ng <- read.table("./GSEA/result_from_server2/esc2016+2018_rb_negative_sn_gsea_report_for_SLOW_1546782445406.tsv", header = T,sep="\t")

head(esc2016_2018_fast_rb_sn_ng)

esc2016_2018_fast_rb_sn_ng <-as.data.frame(esc2016_2018_fast_rb_sn_ng[,c(1,6)])
esc2016_2018_slow_rb_sn_ng <-as.data.frame(esc2016_2018_slow_rb_sn_ng[,c(1,6)])

esc2016_2018_rb_sn_ng <- rbind(esc2016_2018_slow_rb_sn_ng,esc2016_2018_fast_rb_sn_ng)

nrow(esc2016_2018_rb_sn_ng) #14949
# load sox2 data
sox2_EG1_m <- read.table("./G1_M and Fast_Slow/result from server/sox2_EG1_G2_gsea_report_for_EG1_m_1547358456255.tsv", header = T,sep="\t")
sox2_G2_m <- read.table("./G1_M and Fast_Slow/result from server/sox2_EG1_G2_gsea_report_for_G2_m_1547358456255.tsv", header = T,sep="\t")

head(sox2_EG1_m)

sox2_EG1_m <-as.data.frame(sox2_EG1_m[,c(1,6)])
sox2_G2_m <-as.data.frame(sox2_G2_m[,c(1,6)])

sox2 <- rbind(sox2_G2_m,sox2_EG1_m)
nrow(sox2)

#-----plot sox2 vs esc2016+2018 sn ng re-------
sox2_ESC <- cbind(sox2,esc2016_2018_rb_sn_ng[match(sox2$NAME,esc2016_2018_rb_sn_ng$NAME),])
sox2_ESC <-  sox2_ESC[,c(1,2,4)]
head(sox2_ESC)

sox2_ESC <- na.omit(sox2_ESC)
nrow(sox2_ESC)
sox2_ESC$NES.1 <- sox2_ESC$NES.1*(-1)

library(grid)
library(gridExtra)
grob_p = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(sox2_ESC$NES, sox2_ESC$NES.1,method = "pearson"), 4) ), x = 0.5, y = 0.97, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))
grob_s = grobTree(textGrob(paste("Spearman Correlation : ", round(cor(sox2_ESC$NES, sox2_ESC$NES.1,method = "spearman"), 4) ), x = 0.5, y = 0.93, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))



library(ggplot2)
ggplot(sox2_ESC, aes(NES, NES.1))+ 
  geom_point() + 
  
  geom_point(data=sox2_ESC[sox2_ESC$NAME=="HALLMARK_MYC_TARGETS_V1",], aes(NES, NES.1), colour="red", size=3,show.legend = T) +
  
  geom_point(data=sox2_ESC[grep("FISCHER_G.*",sox2_ESC$NAME),], aes(NES, NES.1), colour="blue", size=3) +
  
  geom_point(data=sox2_ESC[grep("WHITFIELD_.*",sox2_ESC$NAME),], aes(NES, NES.1), colour="green", size=3) +
  
  scale_color_gradient(low = "blue", high = "yellow") +
  
  
  #stat_summary(fun.data=mean_cl_normal) + 
  #add a regression line
  geom_smooth(method='lm',formula=y~x,se=TRUE, fullrange=FALSE, level=0.95) +
  #stat_cor(method = "pearson", label.x = 0, label.y = 0) +
  
  annotation_custom(grob_p) +
  
  annotation_custom(grob_s) +
  labs(title = "sox2 vs. ESCs") +
  xlab("sox2.NES") +
  ylab("ESCs.NES") +
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
  axis.title.y = element_text(color="#993333", size=14, face="bold"),
  # axis.title.y = element_blank(),
  # axis.title.x = element_blank(),
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




# load GSE data
GSE_EG1_m <- read.table("./G1_M and Fast_Slow/result from server/GES_EG1_G2_gsea_report_for_EG1_1547361222650.tsv", header = T,sep="\t")
GSE_G2_m <- read.table("./G1_M and Fast_Slow/result from server/GES_EG1_G2_gsea_report_for_G2_1547361222650.tsv", header = T,sep="\t")

head(GSE_EG1_m)

GSE_EG1_m <-as.data.frame(GSE_EG1_m[,c(1,6)])
GSE_G2_m <-as.data.frame(GSE_G2_m[,c(1,6)])

GSE <- rbind(GSE_G2_m,GSE_EG1_m)
nrow(GSE)

#-----plot GSE vs esc2016+2018 sn ng re-------
GSE_ESC <- cbind(GSE,esc2016_2018_rb_sn_ng[match(GSE$NAME,esc2016_2018_rb_sn_ng$NAME),])
GSE_ESC <-  GSE_ESC[,c(2,4)]
head(GSE_ESC)

GSE_ESC <- na.omit(GSE_ESC)
nrow(GSE_ESC)
GSE_ESC$NES.1 <- GSE_ESC$NES.1*(-1)

library(grid)
library(gridExtra)
grob_p = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(GSE_ESC$NES, GSE_ESC$NES.1,method = "pearson"), 4) ), x = 0.5, y = 0.97, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))
grob_s = grobTree(textGrob(paste("Spearman Correlation : ", round(cor(GSE_ESC$NES, GSE_ESC$NES.1,method = "spearman"), 4) ), x = 0.5, y = 0.93, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))



library(ggplot2)
ggplot(GSE_ESC, aes(NES, NES.1))+ 
  geom_point() + 
  
  #stat_summary(fun.data=mean_cl_normal) + 
  #add a regression line
  geom_smooth(method='lm',formula=y~x,se=TRUE, fullrange=FALSE, level=0.95) +
  #stat_cor(method = "pearson", label.x = 0, label.y = 0) +
  
  annotation_custom(grob_p) +
  
  annotation_custom(grob_s) +
  labs(title = "GSE vs. ESCs") +
  xlab("GSE.NES") +
  ylab("ESCs.NES") +
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
  axis.title.y = element_text(color="#993333", size=14, face="bold"),
  # axis.title.y = element_blank(),
  # axis.title.x = element_blank(),
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







#-----plot GSE vs esc2016+2018 sn ng re-------
GSE_sox2 <- cbind(GSE,sox2[match(GSE$NAME,sox2$NAME),])
GSE_sox2 <-  GSE_sox2[,c(2,4)]
head(GSE_sox2)

GSE_sox2 <- na.omit(GSE_sox2)
nrow(GSE_sox2)

library(grid)
library(gridExtra)
grob_p = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(GSE_sox2$NES, GSE_sox2$NES.1,method = "pearson"), 4) ), x = 0.5, y = 0.97, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))
grob_s = grobTree(textGrob(paste("Spearman Correlation : ", round(cor(GSE_sox2$NES, GSE_sox2$NES.1,method = "spearman"), 4) ), x = 0.5, y = 0.93, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))



library(ggplot2)
ggplot(GSE_sox2, aes(NES, NES.1))+ 
  geom_point() + 
  
  geom_point(data=GSE_sox2[10:13, ], aes(NES, NES.1), colour="red", size=5)
  
  #stat_summary(fun.data=mean_cl_normal) + 
  #add a regression line
  geom_smooth(method='lm',formula=y~x,se=TRUE, fullrange=FALSE, level=0.95) +
  #stat_cor(method = "pearson", label.x = 0, label.y = 0) +
  
  annotation_custom(grob_p) +
  
  annotation_custom(grob_s) +
  labs(title = "GSE vs. sox2") +
  xlab("GSE.NES") +
  ylab("sox2.NES") +
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
  axis.title.y = element_text(color="#993333", size=14, face="bold"),
  # axis.title.y = element_blank(),
  # axis.title.x = element_blank(),
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