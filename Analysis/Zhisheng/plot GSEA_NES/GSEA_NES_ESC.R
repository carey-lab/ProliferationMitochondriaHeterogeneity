# plot NES for all the ESC data
# Jan 07 2019
library(readxl)

# load slow data of all log

esc2016_slow_log <- read_xlsx("./GSEA/result_from_server/esc2016_symbol_log_gsea_report_for_SLOW_1546670324514.xls.xlsx", col_names = T)
esc2018_slow_log <- read_xlsx("./GSEA/result_from_server/esc2018_symbol_log_gsea_report_for_SLOW_1546672523895.xls.xlsx", col_names = T)
esc2016_2018_slow_log <- read_xlsx("./GSEA/result_from_server/esc_2016+2018_symbol_log_gsea_report_for_SLOW_1546676884773.xls.xlsx", col_names = T)



head(esc2016_slow_log)
# load fast data of all log

esc2016_fast_log <- read.table("./GSEA/result_from_server/esc2016_symbol_log_gsea_report_for_FAST_1546670324514.tsv", header = T,sep="\t")
esc2018_fast_log <- read.table("./GSEA/result_from_server/esc2018_symbol_log_gsea_report_for_FAST_1546672523895.tsv", header = T,sep="\t")
esc2016_2018_fast_log <- read.table("./GSEA/result_from_server/esc2016+2018_symbol_log_gsea_report_for_FAST_1546676884773.tsv", header = T,sep="\t")



esc2016_slow_log <-as.data.frame(esc2016_slow_log[,c(1,6)])
esc2018_slow_log <-as.data.frame(esc2018_slow_log[,c(1,6)])
esc2016_2018_slow_log <-as.data.frame(esc2016_2018_slow_log[,c(1,6)])


esc2016_fast_log <-as.data.frame(esc2016_fast_log[,c(1,6)])
esc2018_fast_log <-as.data.frame(esc2018_fast_log[,c(1,6)])
esc2016_2018_fast_log <-as.data.frame(esc2016_2018_fast_log[,c(1,6)])


#esc 2016 vs esc 2018 log
head(esc2016_slow_log)
head(esc2018_slow_log)
nrow(esc2016_slow_log) #13091
nrow(esc2018_slow_log) #14152

head(esc2016_fast_log)
head(esc2018_fast_log)
nrow(esc2016_fast_log) #1969
nrow(esc2018_fast_log) #908


esc2016_log <- rbind(esc2016_slow_log,esc2016_fast_log)
esc2018_log <- rbind(esc2018_slow_log,esc2018_fast_log)
esc2016_2018_log <- rbind(esc2016_2018_slow_log,esc2016_2018_fast_log)



#--------plot esc 2016 vs esc 2018 log-----------

esc2016vs2018_log <- cbind(esc2016_log,esc2018_log[match(esc2016_log$NAME,esc2018_log$NAME),])
esc2016vs2018_log <-  esc2016vs2018_log[,c(2,4)]
nrow(esc2016vs2018_log)
esc2016vs2018_log <- na.omit(esc2016vs2018_log)

library(grid)
library(gridExtra)
grob_p = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(esc2016vs2018_log$NES, esc2016vs2018_log$NES.1,method = "pearson"), 4) ), x = 0.5, y = 0.97, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))
grob_s = grobTree(textGrob(paste("Spearman Correlation : ", round(cor(esc2016vs2018_log$NES, esc2016vs2018_log$NES.1,method = "spearman"), 4) ), x = 0.5, y = 0.93, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))



library(ggplot2)
ggplot(esc2016vs2018_log, aes(NES, NES.1))+ 
  geom_point() + 
  
  #stat_summary(fun.data=mean_cl_normal) + 
  #add a regression line
  geom_smooth(method='lm',formula=y~x,se=TRUE, fullrange=FALSE, level=0.95) +
  #stat_cor(method = "pearson", label.x = 0, label.y = 0) +
  
  annotation_custom(grob_p) +
  
  annotation_custom(grob_s) +
  
  
  
  labs(title = "esc2016 vs. esc2018 (log2 not remove) ", x="esc2016 NES", y="esc 2018 NES") +
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
  #delete background
  panel.grid.major =element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  
  #加上坐标轴
  axis.line = element_line(colour = "black",size=1.0),
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
  legend.title = element_text(colour="black", size=14, face="bold")
)









#------plot esc 2016 vs esc 2016+2018 log-------------
esc2016vs2016_2018_log <- cbind(esc2016_log,esc2016_2018_log[match(esc2016_log$NAME,esc2016_2018_log$NAME),])
esc2016vs2016_2018_log <-  esc2016vs2016_2018_log[,c(2,4)]
head(esc2016vs2016_2018_log)
esc2016vs2016_2018_log <- na.omit(esc2016vs2016_2018_log)

library(grid)
library(gridExtra)
grob_p = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(esc2016vs2016_2018_log$NES, esc2016vs2016_2018_log$NES.1,method = "pearson"), 4) ), x = 0.5, y = 0.97, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))
grob_s = grobTree(textGrob(paste("Spearman Correlation : ", round(cor(esc2016vs2016_2018_log$NES, esc2016vs2016_2018_log$NES.1,method = "spearman"), 4) ), x = 0.5, y = 0.93, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))



library(ggplot2)
ggplot(esc2016vs2016_2018_log, aes(NES, NES.1))+ 
  geom_point() + 
  
  #stat_summary(fun.data=mean_cl_normal) + 
  #add a regression line
  geom_smooth(method='lm',formula=y~x,se=TRUE, fullrange=FALSE, level=0.95) +
  #stat_cor(method = "pearson", label.x = 0, label.y = 0) +
  
  annotation_custom(grob_p) +
  
  annotation_custom(grob_s) +
  
  
  
  labs(title = "esc2016 vs. esc2016_2018 (log2 not remove) ", x="esc2016 NES", y="esc 2018 NES") +
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
  #delete background
  panel.grid.major =element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  
  #加上坐标轴
  axis.line = element_line(colour = "black",size=1.0),
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
  legend.title = element_text(colour="black", size=14, face="bold")
)

#------plot esc 2018 vs esc 2016+2018 log-----------

esc2018vs2016_2018_log <- cbind(esc2018_log,esc2016_2018_log[match(esc2018_log$NAME,esc2016_2018_log$NAME),])
esc2018vs2016_2018_log <-  esc2018vs2016_2018_log[,c(2,4)]

esc2018vs2016_2018_log <- na.omit(esc2018vs2016_2018_log)

library(grid)
library(gridExtra)
grob_p = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(esc2018vs2016_2018_log$NES, esc2018vs2016_2018_log$NES.1,method = "pearson"), 4) ), x = 0.5, y = 0.97, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))
grob_s = grobTree(textGrob(paste("Spearman Correlation : ", round(cor(esc2018vs2016_2018_log$NES, esc2018vs2016_2018_log$NES.1,method = "spearman"), 4) ), x = 0.5, y = 0.93, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))



library(ggplot2)
ggplot(esc2018vs2016_2018_log, aes(NES, NES.1))+ 
  geom_point() + 
  
  #stat_summary(fun.data=mean_cl_normal) + 
  #add a regression line
  geom_smooth(method='lm',formula=y~x,se=TRUE, fullrange=FALSE, level=0.95) +
  #stat_cor(method = "pearson", label.x = 0, label.y = 0) +
  
  annotation_custom(grob_p) +
  
  annotation_custom(grob_s) +
  labs(title = "esc2018 vs. esc2016_2018 (log2 not remove) ", x="esc2018 NES", y="esc 2018 NES") +
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
  #delete background
  panel.grid.major =element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  
  #加上坐标轴
  axis.line = element_line(colour = "black",size=1.0),
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
  legend.title = element_text(colour="black", size=14, face="bold")
)









# load not removed 2016+2018 sn data

esc2016_2018_fast_sn <- read.table("./GSEA/result_from_server/esc2016+2018_notre_sn_gsea_report_for_FAST_1546749854327.tsv", header = T,sep="\t")
esc2016_2018_slow_sn <- read.table("./GSEA/result_from_server/esc2016+2018_notre_sn_gsea_report_for_SLOW_1546749854327.tsv", header = T,sep="\t")

esc2016_2018_fast_sn <-as.data.frame(esc2016_2018_fast_sn[,c(1,6)])
esc2016_2018_slow_sn <-as.data.frame(esc2016_2018_slow_sn[,c(1,6)])

esc2016_2018_sn <- rbind(esc2016_2018_slow_sn,esc2016_2018_fast_sn)




#-----plot esc2016 vs esc2016+2018 sn not re-------
esc2016vs2016_2018_sn <- cbind(esc2016_log,esc2016_2018_sn[match(esc2016_log$NAME,esc2016_2018_sn$NAME),])
esc2016vs2016_2018_sn <-  esc2016vs2016_2018_sn[,c(2,4)]
head(esc2016vs2016_2018_sn)

esc2016vs2016_2018_sn <- na.omit(esc2016vs2016_2018_sn)

library(grid)
library(gridExtra)
grob_p = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(esc2016vs2016_2018_sn$NES, esc2016vs2016_2018_sn$NES.1,method = "pearson"), 4) ), x = 0.5, y = 0.97, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))
grob_s = grobTree(textGrob(paste("Spearman Correlation : ", round(cor(esc2016vs2016_2018_sn$NES, esc2016vs2016_2018_sn$NES.1,method = "spearman"), 4) ), x = 0.5, y = 0.93, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))



library(ggplot2)
ggplot(esc2016vs2016_2018_sn, aes(NES, NES.1))+ 
  geom_point() + 
  
  #stat_summary(fun.data=mean_cl_normal) + 
  #add a regression line
  geom_smooth(method='lm',formula=y~x,se=TRUE, fullrange=FALSE, level=0.95) +
  #stat_cor(method = "pearson", label.x = 0, label.y = 0) +
  
  annotation_custom(grob_p) +
  
  annotation_custom(grob_s) +
  labs(title = "esc2016 vs. esc2016_2018 (sn not remove) ", x="esc2016 NES", y="esc 2018 NES") +
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
  #delete background
  panel.grid.major =element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  
  #加上坐标轴
  axis.line = element_line(colour = "black",size=1.0),
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
  legend.title = element_text(colour="black", size=14, face="bold")
)
#-----plot esc2018 vs esc2016+2018 sn not re-------

esc2018vs2016_2018_sn <- cbind(esc2018_log,esc2016_2018_sn[match(esc2018_log$NAME,esc2016_2018_sn$NAME),])
esc2018vs2016_2018_sn <-  esc2018vs2016_2018_sn[,c(2,4)]
head(esc2018vs2016_2018_sn)

esc2018vs2016_2018_sn <- na.omit(esc2018vs2016_2018_sn)

library(grid)
library(gridExtra)
grob_p = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(esc2018vs2016_2018_sn$NES, esc2018vs2016_2018_sn$NES.1,method = "pearson"), 4) ), x = 0.5, y = 0.97, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))
grob_s = grobTree(textGrob(paste("Spearman Correlation : ", round(cor(esc2018vs2016_2018_sn$NES, esc2018vs2016_2018_sn$NES.1,method = "spearman"), 4) ), x = 0.5, y = 0.93, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))



library(ggplot2)
ggplot(esc2018vs2016_2018_sn, aes(NES, NES.1))+ 
  geom_point() + 
  
  #stat_summary(fun.data=mean_cl_normal) + 
  #add a regression line
  geom_smooth(method='lm',formula=y~x,se=TRUE, fullrange=FALSE, level=0.95) +
  #stat_cor(method = "pearson", label.x = 0, label.y = 0) +
  
  annotation_custom(grob_p) +
  
  annotation_custom(grob_s) +
  labs(title = "esc2018 vs. esc2016_2018 (sn not remove) ", x="esc2016 NES", y="esc 2018 NES") +
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
  #delete background
  panel.grid.major =element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  
  #加上坐标轴
  axis.line = element_line(colour = "black",size=1.0),
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
  legend.title = element_text(colour="black", size=14, face="bold")
)















#----------------Jan 07 2019 modify------

# load removed 2016+2018 rb log data-------------


esc2016_2018_fast_rb_log <- read.table("./GSEA/result_from_server2/esc2016+2018_rb_log_gsea_report_for_FAST_1546776448028.tsv", header = T,sep="\t")
esc2016_2018_slow_rb_log <- read.table("./GSEA/result_from_server2/esc2016+2018_rb_log_gsea_report_for_SLOW_1546776448028.tsv", header = T,sep="\t")


esc2016_2018_fast_rb_log <-as.data.frame(esc2016_2018_fast_rb_log[,c(1,6)])
esc2016_2018_slow_rb_log <-as.data.frame(esc2016_2018_slow_rb_log[,c(1,6)])



esc2016_2018_rb_log <- rbind(esc2016_2018_slow_rb_log,esc2016_2018_fast_rb_log)


#----plot esc2016 vs esc2016+2018 log re-----------
esc2016vs2016_2018_rb_log <- cbind(esc2016_log,esc2016_2018_rb_log[match(esc2016_log$NAME,esc2016_2018_rb_log$NAME),])
esc2016vs2016_2018_rb_log <-  esc2016vs2016_2018_rb_log[,c(2,4)]
head(esc2016vs2016_2018_rb_log)
esc2016vs2016_2018_rb_log <- na.omit(esc2016vs2016_2018_rb_log)

library(grid)
library(gridExtra)
grob_p = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(esc2016vs2016_2018_rb_log$NES, esc2016vs2016_2018_rb_log$NES.1,method = "pearson"), 4) ), x = 0.5, y = 0.97, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))
grob_s = grobTree(textGrob(paste("Spearman Correlation : ", round(cor(esc2016vs2016_2018_rb_log$NES, esc2016vs2016_2018_rb_log$NES.1,method = "spearman"), 4) ), x = 0.5, y = 0.93, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))



library(ggplot2)
ggplot(esc2016vs2016_2018_rb_log, aes(NES, NES.1))+ 
  geom_point() + 
  
  #stat_summary(fun.data=mean_cl_normal) + 
  #add a regression line
  geom_smooth(method='lm',formula=y~x,se=TRUE, fullrange=FALSE, level=0.95) +
  #stat_cor(method = "pearson", label.x = 0, label.y = 0) +
  
  annotation_custom(grob_p) +
  
  annotation_custom(grob_s) +
  labs(title = "esc2016 vs. esc2016_2018 (log remove) ", x="esc2016 NES", y="esc 2018 NES") +
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
  #delete background
  panel.grid.major =element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  
  #加上坐标轴
  axis.line = element_line(colour = "black",size=1.0),
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
  legend.title = element_text(colour="black", size=14, face="bold")
)
#----plot esc2018 vs esc2016+2018 log re-----

esc2018vs2016_2018_rb_log <- cbind(esc2018_log,esc2016_2018_rb_log[match(esc2018_log$NAME,esc2016_2018_rb_log$NAME),])
esc2018vs2016_2018_rb_log <-  esc2018vs2016_2018_rb_log[,c(2,4)]
head(esc2018vs2016_2018_rb_log)
esc2018vs2016_2018_rb_log <- na.omit(esc2018vs2016_2018_rb_log)

library(grid)
library(gridExtra)
grob_p = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(esc2018vs2016_2018_rb_log$NES, esc2018vs2016_2018_rb_log$NES.1,method = "pearson"), 4) ), x = 0.5, y = 0.97, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))
grob_s = grobTree(textGrob(paste("Spearman Correlation : ", round(cor(esc2018vs2016_2018_rb_log$NES, esc2018vs2016_2018_rb_log$NES.1,method = "spearman"), 4) ), x = 0.5, y = 0.93, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))



library(ggplot2)
ggplot(esc2018vs2016_2018_rb_log, aes(NES, NES.1))+ 
  geom_point() + 
  
  #stat_summary(fun.data=mean_cl_normal) + 
  #add a regression line
  geom_smooth(method='lm',formula=y~x,se=TRUE, fullrange=FALSE, level=0.95) +
  #stat_cor(method = "pearson", label.x = 0, label.y = 0) +
  
  annotation_custom(grob_p) +
  
  annotation_custom(grob_s) +
  labs(title = "esc2018 vs. esc2016_2018 (log remove) ", x="esc2016 NES", y="esc 2018 NES") +
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
  #delete background
  panel.grid.major =element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  
  #加上坐标轴
  axis.line = element_line(colour = "black",size=1.0),
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
  legend.title = element_text(colour="black", size=14, face="bold")
)





# load removed 2016+2018 rb sn data-------------

esc2016_2018_fast_rb_sn <- read.table("./GSEA/result_from_server2/esc2016+2018_rb_sn_gsea_report_for_FAST_1546776344765.tsv", header = T,sep="\t")
esc2016_2018_slow_rb_sn <- read.table("./GSEA/result_from_server2/esc2016+2018_rb_sn_gsea_report_for_SLOW_1546776344765.tsv", header = T,sep="\t")


esc2016_2018_fast_rb_sn <-as.data.frame(esc2016_2018_fast_rb_sn[,c(1,6)])
esc2016_2018_slow_rb_sn <-as.data.frame(esc2016_2018_slow_rb_sn[,c(1,6)])



esc2016_2018_rb_sn <- rbind(esc2016_2018_slow_rb_sn,esc2016_2018_fast_rb_sn)

#----------plot esc 2016 vs esc 2016+2018 re sn-----
esc2016vs2016_2018_rb_sn <- cbind(esc2016_log,esc2016_2018_rb_sn[match(esc2016_log$NAME,esc2016_2018_rb_sn$NAME),])
esc2016vs2016_2018_rb_sn <-  esc2016vs2016_2018_rb_sn[,c(2,4)]
head(esc2016vs2016_2018_rb_sn)

esc2016vs2016_2018_rb_sn <- na.omit(esc2016vs2016_2018_rb_sn)

library(grid)
library(gridExtra)
grob_p = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(esc2016vs2016_2018_rb_sn$NES, esc2016vs2016_2018_rb_sn$NES.1,method = "pearson"), 4) ), x = 0.5, y = 0.97, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))
grob_s = grobTree(textGrob(paste("Spearman Correlation : ", round(cor(esc2016vs2016_2018_rb_sn$NES, esc2016vs2016_2018_rb_sn$NES.1,method = "spearman"), 4) ), x = 0.5, y = 0.93, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))



library(ggplot2)
ggplot(esc2016vs2016_2018_rb_sn, aes(NES, NES.1))+ 
  geom_point() + 
  
  #stat_summary(fun.data=mean_cl_normal) + 
  #add a regression line
  geom_smooth(method='lm',formula=y~x,se=TRUE, fullrange=FALSE, level=0.95) +
  #stat_cor(method = "pearson", label.x = 0, label.y = 0) +
  
  annotation_custom(grob_p) +
  
  annotation_custom(grob_s) +
  labs(title = "esc2016 vs. esc2016_2018 (sn remove) ", x="esc2016 NES", y="esc 2018 NES") +
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
  #delete background
  panel.grid.major =element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  
  #加上坐标轴
  axis.line = element_line(colour = "black",size=1.0),
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
  legend.title = element_text(colour="black", size=14, face="bold")
)


#----------plot esc 2018 vs esc 2016+2018 re sn--------------
esc2018vs2016_2018_rb_sn <- cbind(esc2018_log,esc2016_2018_rb_sn[match(esc2018_log$NAME,esc2016_2018_rb_sn$NAME),])
esc2018vs2016_2018_rb_sn <-  esc2018vs2016_2018_rb_sn[,c(2,4)]
head(esc2018vs2016_2018_rb_sn)
esc2018vs2016_2018_rb_sn <- na.omit(esc2018vs2016_2018_rb_sn)

library(grid)
library(gridExtra)
grob_p = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(esc2018vs2016_2018_rb_sn$NES, esc2018vs2016_2018_rb_sn$NES.1,method = "pearson"), 4) ), x = 0.5, y = 0.97, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))
grob_s = grobTree(textGrob(paste("Spearman Correlation : ", round(cor(esc2018vs2016_2018_rb_sn$NES, esc2018vs2016_2018_rb_sn$NES.1,method = "spearman"), 4) ), x = 0.5, y = 0.93, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))



library(ggplot2)
ggplot(esc2018vs2016_2018_rb_sn, aes(NES, NES.1))+ 
  geom_point() + 
  
  #stat_summary(fun.data=mean_cl_normal) + 
  #add a regression line
  geom_smooth(method='lm',formula=y~x,se=TRUE, fullrange=FALSE, level=0.95) +
  #stat_cor(method = "pearson", label.x = 0, label.y = 0) +
  
  annotation_custom(grob_p) +
  
  annotation_custom(grob_s) +
  labs(title = "esc2018 vs. esc2016_2018 (sn remove) ", x="esc2016 NES", y="esc 2018 NES") +
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
  #delete background
  panel.grid.major =element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  
  #加上坐标轴
  axis.line = element_line(colour = "black",size=1.0),
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
  legend.title = element_text(colour="black", size=14, face="bold")
)



# as we can see the correlation between esc 2016 or 2018 with esc 2016+2018(re-batch) will be lower, so comapareesc 2016 or 2018(re batch ) with 2016+2018(re-batch)

# load esc2016 and 2018 re batch
esc2016_fast_rb_log <- read.table("./GSEA/result_from_server2/esc2016_rb_negative_log_gsea_report_for_FAST_1546787546324.tsv", header = T,sep="\t")

#A TEST ES VS NES-----------------------
head(esc2016_fast_rb_log_test)
esc2016_fast_rb_log_test <-as.data.frame(esc2016_fast_rb_log[,c(5,6)])
esc2016_fast_rb_log_test[is.na(esc2016_fast_rb_log_test)] <- 0


grob_p = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(esc2016_fast_rb_log_test$ES, esc2016_fast_rb_log_test$NES,method = "pearson"), 4) ), x = 0, y = 0.97, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))
grob_s = grobTree(textGrob(paste("Spearman Correlation : ", round(cor(esc2016_fast_rb_log_test$ES, esc2016_fast_rb_log_test$NES,method = "spearman"), 4) ), x = 0, y = 0.93, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))



library(ggplot2)
ggplot(esc2016_fast_rb_log_test, aes(ES, NES))+ 
  geom_point() + 
  
  #stat_summary(fun.data=mean_cl_normal) + 
  #add a regression line
  geom_smooth(method='lm',formula=y~x,se=TRUE, fullrange=FALSE, level=0.95) +
  #stat_cor(method = "pearson", label.x = 0, label.y = 0) +
  
  annotation_custom(grob_p) +
  
  annotation_custom(grob_s) +
  
  
  
  labs(title = "ES vs. NES", x="esc2016 NES", y="esc 2018 NES") +
 
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
  #delete background
  panel.grid.major =element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  
  #加上坐标轴
  axis.line = element_line(colour = "black",size=1.0),
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
  legend.title = element_text(colour="black", size=14, face="bold")
)






















esc2016_slow_rb_log <- read.table("./GSEA/result_from_server2/esc2016_rb_negative_log_gsea_report_for_SLOW_1546787546324.tsv", header = T,sep="\t")

esc2018_fast_rb_log <- read.table("./GSEA/result_from_server2/esc2018_rb_negativegsea_report_for_FAST_1546785095675.tsv", header = T,sep="\t")
esc2018_slow_rb_log <- read.table("./GSEA/result_from_server2/esc2018_rb_negative_log_gsea_report_for_SLOW_1546785095675.tsv", header = T,sep="\t")

esc2016_fast_rb_log <-as.data.frame(esc2016_fast_rb_log[,c(1,6)])
esc2016_slow_rb_log <-as.data.frame(esc2016_slow_rb_log[,c(1,6)])
esc2018_fast_rb_log <-as.data.frame(esc2018_fast_rb_log[,c(1,6)])
esc2018_slow_rb_log <-as.data.frame(esc2018_slow_rb_log[,c(1,6)])


esc2016_rb_log <- rbind(esc2016_slow_rb_log,esc2016_fast_rb_log)
esc2018_rb_log <- rbind(esc2018_slow_rb_log,esc2018_fast_rb_log)

#----------plot esc 2016 vs 2108 log re batch--------------


esc2016vs2018_rb_log <- cbind(esc2016_rb_log,esc2018_rb_log[match(esc2016_rb_log$NAME,esc2018_rb_log$NAME),])
esc2016vs2018_rb_log <-  esc2016vs2018_rb_log[,c(2,4)]
nrow(esc2016vs2018_rb_log)
esc2016vs2018_rb_log <- na.omit(esc2016vs2018_rb_log)

library(grid)
library(gridExtra)
grob_p = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(esc2016vs2018_rb_log$NES, esc2016vs2018_rb_log$NES.1,method = "pearson"), 4) ), x = 0.5, y = 0.97, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))
grob_s = grobTree(textGrob(paste("Spearman Correlation : ", round(cor(esc2016vs2018_rb_log$NES, esc2016vs2018_rb_log$NES.1,method = "spearman"), 4) ), x = 0.5, y = 0.93, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))



library(ggplot2)
ggplot(esc2016vs2018_rb_log, aes(NES, NES.1))+ 
  geom_point() + 
  
  #stat_summary(fun.data=mean_cl_normal) + 
  #add a regression line
  geom_smooth(method='lm',formula=y~x,se=TRUE, fullrange=FALSE, level=0.95) +
  #stat_cor(method = "pearson", label.x = 0, label.y = 0) +
  
  annotation_custom(grob_p) +
  
  annotation_custom(grob_s) +
  
  
  
  labs(title = "esc2016 vs. esc2018 (log2 removed) ", x="esc2016 NES", y="esc 2018 NES") +
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
  #delete background
  panel.grid.major =element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  
  #加上坐标轴
  axis.line = element_line(colour = "black",size=1.0),
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
  legend.title = element_text(colour="black", size=14, face="bold")
)







#-----------plot esc 2016 re  vs esc 2016+2018 re log batch---------

esc2016vs2016_2018_rb_log <- cbind(esc2016_rb_log,esc2016_2018_rb_log[match(esc2016_rb_log$NAME,esc2016_2018_rb_log$NAME),])
esc2016vs2016_2018_rb_log <-  esc2016vs2016_2018_rb_log[,c(2,4)]
nrow(esc2016vs2016_2018_rb_log)
esc2016vs2016_2018_rb_log <- na.omit(esc2016vs2016_2018_rb_log)

library(grid)
library(gridExtra)
grob_p = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(esc2016vs2016_2018_rb_log$NES, esc2016vs2016_2018_rb_log$NES.1,method = "pearson"), 4) ), x = 0.5, y = 0.97, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))
grob_s = grobTree(textGrob(paste("Spearman Correlation : ", round(cor(esc2016vs2016_2018_rb_log$NES, esc2016vs2016_2018_rb_log$NES.1,method = "spearman"), 4) ), x = 0.5, y = 0.93, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))



library(ggplot2)
ggplot(esc2016vs2016_2018_rb_log, aes(NES, NES.1))+ 
  geom_point() + 
  
  #stat_summary(fun.data=mean_cl_normal) + 
  #add a regression line
  geom_smooth(method='lm',formula=y~x,se=TRUE, fullrange=FALSE, level=0.95) +
  #stat_cor(method = "pearson", label.x = 0, label.y = 0) +
  
  annotation_custom(grob_p) +
  
  annotation_custom(grob_s) +
  
  
  
  labs(title = "esc2016 vs. esc2016+2018 (log2 removed) ", x="esc2016 NES", y="esc 2018 NES") +
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
  #delete background
  panel.grid.major =element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  
  #加上坐标轴
  axis.line = element_line(colour = "black",size=1.0),
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
  legend.title = element_text(colour="black", size=14, face="bold")
)

#-----------plot esc 2018 re  vs esc 2016+2018 re log batch-----------

esc2018vs2016_2018_rb_log <- cbind(esc2018_rb_log,esc2016_2018_rb_log[match(esc2018_rb_log$NAME,esc2016_2018_rb_log$NAME),])
esc2018vs2016_2018_rb_log <-  esc2018vs2016_2018_rb_log[,c(2,4)]
nrow(esc2018vs2016_2018_rb_log)
esc2018vs2016_2018_rb_log <- na.omit(esc2018vs2016_2018_rb_log)

library(grid)
library(gridExtra)
grob_p = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(esc2018vs2016_2018_rb_log$NES, esc2018vs2016_2018_rb_log$NES.1,method = "pearson"), 4) ), x = 0.5, y = 0.97, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))
grob_s = grobTree(textGrob(paste("Spearman Correlation : ", round(cor(esc2018vs2016_2018_rb_log$NES, esc2018vs2016_2018_rb_log$NES.1,method = "spearman"), 4) ), x = 0.5, y = 0.93, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))



library(ggplot2)
ggplot(esc2018vs2016_2018_rb_log, aes(NES, NES.1))+ 
  geom_point() + 
  
  #stat_summary(fun.data=mean_cl_normal) + 
  #add a regression line
  geom_smooth(method='lm',formula=y~x,se=TRUE, fullrange=FALSE, level=0.95) +
  #stat_cor(method = "pearson", label.x = 0, label.y = 0) +
  
  annotation_custom(grob_p) +
  
  annotation_custom(grob_s) +
  
  
  
  labs(title = "esc2018 vs. esc2016+2018 (log2 removed) ", x="esc2016 NES", y="esc 2018 NES") +
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
  #delete background
  panel.grid.major =element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  
  #加上坐标轴
  axis.line = element_line(colour = "black",size=1.0),
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
  legend.title = element_text(colour="black", size=14, face="bold")
)


















# load removed 2016+2018 rb sn data-------------


esc2016_2018_fast_rb_sn <- read.table("./GSEA/result_from_server2/esc2016+2018_rb_sn_gsea_report_for_FAST_1546776344765.tsv", header = T,sep="\t")
esc2016_2018_slow_rb_sn <- read.table("./GSEA/result_from_server2/esc2016+2018_rb_sn_gsea_report_for_SLOW_1546776344765.tsv", header = T,sep="\t")


esc2016_2018_fast_rb_sn <-as.data.frame(esc2016_2018_fast_rb_sn[,c(1,6)])
esc2016_2018_slow_rb_sn <-as.data.frame(esc2016_2018_slow_rb_sn[,c(1,6)])


esc2016_2018_rb_sn <- rbind(esc2016_2018_slow_rb_sn,esc2016_2018_fast_rb_sn)

#------------plot esc 2016 re  vs esc 2016+2018 re sn batch-----
esc2016vs2016_2018_rb_sn <- cbind(esc2016_rb_log,esc2016_2018_rb_sn[match(esc2016_rb_log$NAME,esc2016_2018_rb_sn$NAME),])
esc2016vs2016_2018_rb_sn <-  esc2016vs2016_2018_rb_sn[,c(2,4)]
nrow(esc2016vs2016_2018_rb_sn)
esc2016vs2016_2018_rb_sn <- na.omit(esc2016vs2016_2018_rb_sn)

library(grid)
library(gridExtra)
grob_p = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(esc2016vs2016_2018_rb_sn$NES, esc2016vs2016_2018_rb_sn$NES.1,method = "pearson"), 4) ), x = 0.5, y = 0.97, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))
grob_s = grobTree(textGrob(paste("Spearman Correlation : ", round(cor(esc2016vs2016_2018_rb_sn$NES, esc2016vs2016_2018_rb_sn$NES.1,method = "spearman"), 4) ), x = 0.5, y = 0.93, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))



library(ggplot2)
ggplot(esc2016vs2016_2018_rb_sn, aes(NES, NES.1))+ 
  geom_point() + 
  
  #stat_summary(fun.data=mean_cl_normal) + 
  #add a regression line
  geom_smooth(method='lm',formula=y~x,se=TRUE, fullrange=FALSE, level=0.95) +
  #stat_cor(method = "pearson", label.x = 0, label.y = 0) +
  
  annotation_custom(grob_p) +
  
  annotation_custom(grob_s) +
  
  
  
  labs(title = "esc2016 vs. esc2016+2018 (sn removed) ", x="esc2016 NES", y="esc 2018 NES") +
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
  #delete background
  panel.grid.major =element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  
  #加上坐标轴
  axis.line = element_line(colour = "black",size=1.0),
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
  legend.title = element_text(colour="black", size=14, face="bold")
)


#---plot 2018 vs 2016+2018 sn re----------------

esc2018vs2016_2018_rb_sn <- cbind(esc2018_rb_log,esc2016_2018_rb_sn[match(esc2018_rb_log$NAME,esc2016_2018_rb_sn$NAME),])
esc2018vs2016_2018_rb_sn <-  esc2018vs2016_2018_rb_sn[,c(2,4)]
nrow(esc2018vs2016_2018_rb_sn)
esc2018vs2016_2018_rb_sn <- na.omit(esc2018vs2016_2018_rb_sn)

library(grid)
library(gridExtra)
grob_p = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(esc2018vs2016_2018_rb_sn$NES, esc2018vs2016_2018_rb_sn$NES.1,method = "pearson"), 4) ), x = 0.5, y = 0.97, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))
grob_s = grobTree(textGrob(paste("Spearman Correlation : ", round(cor(esc2018vs2016_2018_rb_sn$NES, esc2018vs2016_2018_rb_sn$NES.1,method = "spearman"), 4) ), x = 0.5, y = 0.93, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))



library(ggplot2)
ggplot(esc2018vs2016_2018_rb_sn, aes(NES, NES.1))+ 
  geom_point() + 
  
  #stat_summary(fun.data=mean_cl_normal) + 
  #add a regression line
  geom_smooth(method='lm',formula=y~x,se=TRUE, fullrange=FALSE, level=0.95) +
  #stat_cor(method = "pearson", label.x = 0, label.y = 0) +
  
  annotation_custom(grob_p) +
  
  annotation_custom(grob_s) +
  
  
  
  labs(title = "esc2018 vs. esc2016+2018 (log2 removed) ", x="esc2016 NES", y="esc 2018 NES") +
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
  #delete background
  panel.grid.major =element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  
  #加上坐标轴
  axis.line = element_line(colour = "black",size=1.0),
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
  legend.title = element_text(colour="black", size=14, face="bold")
)



#--load esc rb sn ng------

esc2016_2018_fast_rb_sn_ng <- read.table("./GSEA/result_from_server2/esc2016+2018_rb_negative_sn_gsea_report_for_FAST_1546782445406.tsv", header = T,sep="\t")
esc2016_2018_slow_rb_sn_ng <- read.table("./GSEA/result_from_server2/esc2016+2018_rb_negative_sn_gsea_report_for_SLOW_1546782445406.tsv", header = T,sep="\t")

head(esc2016_2018_fast_rb_sn_ng)

esc2016_2018_fast_rb_sn_ng <-as.data.frame(esc2016_2018_fast_rb_sn_ng[,c(1,6)])
esc2016_2018_slow_rb_sn_ng <-as.data.frame(esc2016_2018_slow_rb_sn_ng[,c(1,6)])


esc2016_2018_rb_sn_ng <- rbind(esc2016_2018_slow_rb_sn_ng,esc2016_2018_fast_rb_sn_ng)



# plot re sn vs re sn ng-------------------
esc_re_snvsesc_re_sn_ng <- cbind(esc2016_2018_rb_sn,esc2016_2018_rb_sn_ng[match(esc2016_2018_rb_sn$NAME,esc2016_2018_rb_sn_ng$NAME),])
esc_re_snvsesc_re_sn_ng <-  esc_re_snvsesc_re_sn_ng[,c(2,4)]
nrow(esc_re_snvsesc_re_sn_ng)
esc_re_snvsesc_re_sn_ng <- na.omit(esc_re_snvsesc_re_sn_ng)

library(grid)
library(gridExtra)
grob_p = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(esc_re_snvsesc_re_sn_ng$NES, esc_re_snvsesc_re_sn_ng$NES.1,method = "pearson"), 4) ), x = 0.5, y = 0.97, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))
grob_s = grobTree(textGrob(paste("Spearman Correlation : ", round(cor(esc_re_snvsesc_re_sn_ng$NES, esc_re_snvsesc_re_sn_ng$NES.1,method = "spearman"), 4) ), x = 0.5, y = 0.93, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))



library(ggplot2)
ggplot(esc_re_snvsesc_re_sn_ng, aes(NES, NES.1))+ 
  geom_point() + 
  
  #stat_summary(fun.data=mean_cl_normal) + 
  #add a regression line
  geom_smooth(method='lm',formula=y~x,se=TRUE, fullrange=FALSE, level=0.95) +
  #stat_cor(method = "pearson", label.x = 0, label.y = 0) +
  
  annotation_custom(grob_p) +
  
  annotation_custom(grob_s) +
  
  
  
  labs(title = "esc2016+2018 vs. esc2016+2018_ng (sn removed) ", x="esc2016 NES", y="esc 2018 NES") +
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
  #delete background
  panel.grid.major =element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  
  #加上坐标轴
  axis.line = element_line(colour = "black",size=1.0),
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
  legend.title = element_text(colour="black", size=14, face="bold")
)








# compare fib with esc---------------------------


fib_rb_sn_ngvsesc_rb_sn_ng <- cbind(fib2016_2018_rb_sn_ng,esc2016_2018_rb_sn_ng[match(fib2016_2018_rb_sn_ng$NAME,esc2016_2018_rb_sn_ng$NAME),])
fib_rb_sn_ngvsesc_rb_sn_ng <-  fib_rb_sn_ngvsesc_rb_sn_ng[,c(2,4)]
nrow(fib_rb_sn_ngvsesc_rb_sn_ng)
fib_rb_sn_ngvsesc_rb_sn_ng <- na.omit(fib_rb_sn_ngvsesc_rb_sn_ng)

library(grid)
library(gridExtra)
grob_p = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(fib_rb_sn_ngvsesc_rb_sn_ng$NES, fib_rb_sn_ngvsesc_rb_sn_ng$NES.1,method = "pearson"), 4) ), x = 0.5, y = 0.97, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))
grob_s = grobTree(textGrob(paste("Spearman Correlation : ", round(cor(fib_rb_sn_ngvsesc_rb_sn_ng$NES, fib_rb_sn_ngvsesc_rb_sn_ng$NES.1,method = "spearman"), 4) ), x = 0.5, y = 0.93, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))



library(ggplot2)
ggplot(fib_rb_sn_ngvsesc_rb_sn_ng, aes(NES, NES.1))+ 
  geom_point() + 
  
  #stat_summary(fun.data=mean_cl_normal) + 
  #add a regression line
  geom_smooth(method='lm',formula=y~x,se=TRUE, fullrange=FALSE, level=0.95) +
  #stat_cor(method = "pearson", label.x = 0, label.y = 0) +
  
  annotation_custom(grob_p) +
  
  annotation_custom(grob_s) +
  
  
  
  labs(title = "fib_rb_sn_ng vs. esc_rn_sn_ng (sn removed) ", x="esc2016 NES", y="esc 2018 NES") +
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
  #delete background
  panel.grid.major =element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  
  #加上坐标轴
  axis.line = element_line(colour = "black",size=1.0),
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
  legend.title = element_text(colour="black", size=14, face="bold")
)







#------find 4 top 50-----------------
#load esc2016_2018_rb_sn_ng
esc2016_2018_fast_rb_sn_ng <- read.table("./GSEA/result_from_server2/esc2016+2018_rb_negative_sn_gsea_report_for_FAST_1546782445406.tsv", header = T,sep="\t")
esc2016_2018_slow_rb_sn_ng <- read.table("./GSEA/result_from_server2/esc2016+2018_rb_negative_sn_gsea_report_for_SLOW_1546782445406.tsv", header = T,sep="\t")

head(esc2016_2018_fast_rb_sn_ng)

esc2016_2018_fast_rb_sn_ng <-as.data.frame(esc2016_2018_fast_rb_sn_ng[,c(1,6)])
esc2016_2018_slow_rb_sn_ng <-as.data.frame(esc2016_2018_slow_rb_sn_ng[,c(1,6)])


esc2016_2018_rb_sn_ng <- rbind(esc2016_2018_slow_rb_sn_ng,esc2016_2018_fast_rb_sn_ng)
esc2016_2018_rb_sn_ng <- na.omit(esc2016_2018_rb_sn_ng)
nrow(esc2016_2018_rb_sn_ng)
#load fib2016_2018_rb_sn_ng

fib2016_2018_fast_rb_sn_ng <- read.table("./GSEA/result_from_server2/fib2016+2018_rb_ng_sn_gsea_report_for_FAST_1546779207404.tsv", header = T,sep="\t")
fib2016_2018_slow_rb_sn_ng <- read.table("./GSEA/result_from_server2/fib2016+2018_rb_ng_sn_gsea_report_for_SLOW_1546779207404.tsv", header = T,sep="\t")



fib2016_2018_fast_rb_sn_ng <-as.data.frame(fib2016_2018_fast_rb_sn_ng[,c(1,6)])
fib2016_2018_slow_rb_sn_ng <-as.data.frame(fib2016_2018_slow_rb_sn_ng[,c(1,6)])


fib2016_2018_rb_sn_ng <- rbind(fib2016_2018_slow_rb_sn_ng,fib2016_2018_fast_rb_sn_ng)
fib2016_2018_rb_sn_ng <- na.omit(fib2016_2018_rb_sn_ng)
nrow(fib2016_2018_rb_sn_ng)

fib_esc <- cbind(fib2016_2018_rb_sn_ng,esc2016_2018_rb_sn_ng$NES[match(fib2016_2018_rb_sn_ng$NAME,esc2016_2018_rb_sn_ng$NAME)])

nrow(fib_esc)

fib_esc <- na.omit(fib_esc)


# top positive 
head(fib_esc)
colnames(fib_esc) <- c("NAME","NES_FIB","NES_ESC")
fib_esc$top_pos <- fib_esc$NES_FIB+fib_esc$NES_ESC
fib_esc$top_neg <- fib_esc$NES_FIB-fib_esc$NES_ESC


top_pos <- fib_esc[,c(1,4)][order(fib_esc$top_pos,decreasing = T),]

top_neg <- fib_esc[,c(1,5)][order(fib_esc$top_neg,decreasing = T),]

top_pos_plus_100 <- head(top_pos,100)

top_pos_minus_100 <- tail(top_pos,100)

top_neg_plus_100 <- head(top_neg,100)

top_neg_minus_100 <- tail(top_neg,100)

top_100 <- cbind(top_pos_plus_100,top_pos_minus_100,top_neg_plus_100,top_neg_minus_100)

colnames(top_100) <- c("name","top_pos_plus_100","name","top_pos_minus_100","name","top_neg_plus_100","name","top_neg_minus_100")

head(top_100)

write.table(top_100,file="NES_top_100.tsv",col.names = T,row.names = F,quote = F,sep = "\t")
