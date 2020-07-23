# plot NES for all the FIB data
# Jan 06 2019
library(readxl)

# load slow data of all log
fib2016_slow_log <- read_xlsx("./GSEA/result_from_server/fib2016_symbol_log_gsea_report_for_SLOW_1546667114972.xls.xlsx", col_names = T)
fib2018_slow_log <- read_xlsx("./GSEA/result_from_server/fib2018_symbol_log_gsea_report_for_SLOW_1546667957369.xls.xlsx", col_names = T)
fib2016_2018_slow_log <- read_xlsx("./GSEA/result_from_server/fib2016+2018_symbol_log_gsea_report_for_SLOW_1546674692563.xls.xlsx", col_names = T)
esc2016_slow_log <- read_xlsx("./GSEA/result_from_server/esc2016_symbol_log_gsea_report_for_SLOW_1546670324514.xls.xlsx", col_names = T)
esc2018_slow_log <- read_xlsx("./GSEA/result_from_server/esc2018_symbol_log_gsea_report_for_SLOW_1546672523895.xls.xlsx", col_names = T)
esc2016_2018_slow_log <- read_xlsx("./GSEA/result_from_server/esc_2016+2018_symbol_log_gsea_report_for_SLOW_1546676884773.xls.xlsx", col_names = T)



head(fib2016_slow_log)
# 

# load fast data of all log
fib2016_fast_log <- read.table("./GSEA/result_from_server/fib2016_symbol_log_gsea_report_for_FAST_1546667114972.tsv", header = T,sep="\t")
head(fib2016_fast_log)
fib2018_fast_log <- read.table("./GSEA/result_from_server/fib2018_symbol_log_gsea_report_for_FAST_1546667957369.tsv", header = T,sep="\t")
fib2016_2018_fast_log <- read.table("./GSEA/result_from_server/fib2016+2018_symbol_log_gsea_report_for_FAST_1546674692563.tsv", header = T,sep="\t")
esc2016_fast_log <- read.table("./GSEA/result_from_server/esc2016_symbol_log_gsea_report_for_FAST_1546670324514.tsv", header = T,sep="\t")
esc2018_fast_log <- read.table("./GSEA/result_from_server/esc2018_symbol_log_gsea_report_for_FAST_1546672523895.tsv", header = T,sep="\t")
esc2016_2018_fast_log <- read.table("./GSEA/result_from_server/esc2016+2018_symbol_log_gsea_report_for_FAST_1546676884773.tsv", header = T,sep="\t")




fib2016_slow_log <-as.data.frame(fib2016_slow_log[,c(1,6)])
length(fib2016_slow_log$NAME)
length(fib2016_slow_log$NES)
fib2018_slow_log <-as.data.frame(fib2018_slow_log[,c(1,6)])
fib2016_2018_slow_log <-as.data.frame(fib2016_2018_slow_log[,c(1,6)])
esc2016_slow_log <-as.data.frame(esc2016_slow_log[,c(1,6)])
esc2018_slow_log <-as.data.frame(esc2018_slow_log[,c(1,6)])
esc2016_2018_slow_log <-as.data.frame(esc2016_2018_slow_log[,c(1,6)])



fib2016_fast_log <-as.data.frame(fib2016_fast_log[,c(1,6)])
length(fib2016_fast_log$NAME)
length(fib2016_fast_log$NES)
fib2018_fast_log <-as.data.frame(fib2018_fast_log[,c(1,6)])
fib2016_2018_fast_log <-as.data.frame(fib2016_2018_fast_log[,c(1,6)])
esc2016_fast_log <-as.data.frame(esc2016_fast_log[,c(1,6)])
esc2018_fast_log <-as.data.frame(esc2018_fast_log[,c(1,6)])
esc2016_2018_fast_log <-as.data.frame(esc2016_2018_fast_log[,c(1,6)])


#fib 2016 vs fib 2018 log
head(fib2016_slow_log)
head(fib2018_slow_log)
nrow(fib2016_slow_log) #13091
nrow(fib2018_slow_log) #14152

head(fib2016_fast_log)
head(fib2018_fast_log)
nrow(fib2016_fast_log) #1969
nrow(fib2018_fast_log) #908


fib2016_log <- rbind(fib2016_slow_log,fib2016_fast_log)
fib2018_log <- rbind(fib2018_slow_log,fib2018_fast_log)
fib2016_2018_log <- rbind(fib2016_2018_slow_log,fib2016_2018_fast_log)



#--------plot fib 2016 vs fib 2018 log-----------

fib2016vs2018_log <- cbind(fib2016_log,fib2018_log[match(fib2016_log$NAME,fib2018_log$NAME),])
fib2016vs2018_log <-  fib2016vs2018_log[,c(2,4)]
nrow(fib2016vs2018_log)
fib2016vs2018_log <- na.omit(fib2016vs2018_log)

library(grid)
library(gridExtra)
grob_p = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(fib2016vs2018_log$NES, fib2016vs2018_log$NES.1,method = "pearson"), 4) ), x = 0.5, y = 0.97, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))
grob_s = grobTree(textGrob(paste("Spearman Correlation : ", round(cor(fib2016vs2018_log$NES, fib2016vs2018_log$NES.1,method = "spearman"), 4) ), x = 0.5, y = 0.93, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))



library(ggplot2)
ggplot(fib2016vs2018_log, aes(NES, NES.1))+ 
  geom_point() + 
  
  #stat_summary(fun.data=mean_cl_normal) + 
  #add a regression line
  geom_smooth(method='lm',formula=y~x,se=TRUE, fullrange=FALSE, level=0.95) +
  #stat_cor(method = "pearson", label.x = 0, label.y = 0) +
  
  annotation_custom(grob_p) +
  
  annotation_custom(grob_s) +
  
  
  
  labs(title = "fib2016 vs. fib2018 (log2 not remove) ", x="fib2016 NES", y="fib 2018 NES") +
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
  








#------plot fib 2016 vs fib 2016+2018 log-------------
fib2016vs2016_2018_log <- cbind(fib2016_log,fib2016_2018_log[match(fib2016_log$NAME,fib2016_2018_log$NAME),])
fib2016vs2016_2018_log <-  fib2016vs2016_2018_log[,c(2,4)]
head(fib2016vs2016_2018_log)
fib2016vs2016_2018_log <- na.omit(fib2016vs2016_2018_log)

library(grid)
library(gridExtra)
grob_p = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(fib2016vs2016_2018_log$NES, fib2016vs2016_2018_log$NES.1,method = "pearson"), 4) ), x = 0.5, y = 0.97, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))
grob_s = grobTree(textGrob(paste("Spearman Correlation : ", round(cor(fib2016vs2016_2018_log$NES, fib2016vs2016_2018_log$NES.1,method = "spearman"), 4) ), x = 0.5, y = 0.93, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))



library(ggplot2)
ggplot(fib2016vs2016_2018_log, aes(NES, NES.1))+ 
  geom_point() + 
  
  #stat_summary(fun.data=mean_cl_normal) + 
  #add a regression line
  geom_smooth(method='lm',formula=y~x,se=TRUE, fullrange=FALSE, level=0.95) +
  #stat_cor(method = "pearson", label.x = 0, label.y = 0) +
  
  annotation_custom(grob_p) +
  
  annotation_custom(grob_s) +
  
  
  
  labs(title = "fib2016 vs. fib2016_2018 (log2 not remove) ", x="fib2016 NES", y="fib 2018 NES") +
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

#------plot fib 2018 vs fib 2016+2018 log-----------

fib2018vs2016_2018_log <- cbind(fib2018_log,fib2016_2018_log[match(fib2018_log$NAME,fib2016_2018_log$NAME),])
fib2018vs2016_2018_log <-  fib2018vs2016_2018_log[,c(2,4)]

fib2018vs2016_2018_log <- na.omit(fib2018vs2016_2018_log)

library(grid)
library(gridExtra)
grob_p = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(fib2018vs2016_2018_log$NES, fib2018vs2016_2018_log$NES.1,method = "pearson"), 4) ), x = 0.5, y = 0.97, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))
grob_s = grobTree(textGrob(paste("Spearman Correlation : ", round(cor(fib2018vs2016_2018_log$NES, fib2018vs2016_2018_log$NES.1,method = "spearman"), 4) ), x = 0.5, y = 0.93, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))



library(ggplot2)
ggplot(fib2018vs2016_2018_log, aes(NES, NES.1))+ 
  geom_point() + 
  
  #stat_summary(fun.data=mean_cl_normal) + 
  #add a regression line
  geom_smooth(method='lm',formula=y~x,se=TRUE, fullrange=FALSE, level=0.95) +
  #stat_cor(method = "pearson", label.x = 0, label.y = 0) +
  
  annotation_custom(grob_p) +
  
  annotation_custom(grob_s) +
  labs(title = "fib2018 vs. fib2016_2018 (log2 not remove) ", x="fib2018 NES", y="fib 2018 NES") +
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
fib2016_2018_fast_sn <- read.table("./GSEA/result_from_server/fib2016+2018_notre_sn_gsea_report_for_FAST_1546747118845.tsv", header = T,sep="\t")
fib2016_2018_slow_sn <- read.table("./GSEA/result_from_server/fib2016+2018_notre_sn_gsea_report_for_SLOW_1546747118845.tsv", header = T,sep="\t")

esc2016_2018_fast_sn <- read.table("./GSEA/result_from_server/esc2016+2018_notre_sn_gsea_report_for_FAST_1546749854327.tsv", header = T,sep="\t")
esc2016_2018_slow_sn <- read.table("./GSEA/result_from_server/esc2016+2018_notre_sn_gsea_report_for_SLOW_1546749854327.tsv", header = T,sep="\t")

fib2016_2018_fast_sn <-as.data.frame(fib2016_2018_fast_sn[,c(1,6)])
fib2016_2018_slow_sn <-as.data.frame(fib2016_2018_slow_sn[,c(1,6)])
esc2016_2018_fast_sn <-as.data.frame(esc2016_2018_fast_sn[,c(1,6)])
esc2016_2018_slow_sn <-as.data.frame(esc2016_2018_slow_sn[,c(1,6)])


fib2016_2018_sn <- rbind(fib2016_2018_slow_sn,fib2016_2018_fast_sn)
esc2016_2018_sn <- rbind(esc2016_2018_slow_sn,esc2016_2018_fast_sn)




#-----plot fib2016 vs fib2016+2018 sn not re-------
fib2016vs2016_2018_sn <- cbind(fib2016_log,fib2016_2018_sn[match(fib2016_log$NAME,fib2016_2018_sn$NAME),])
fib2016vs2016_2018_sn <-  fib2016vs2016_2018_sn[,c(2,4)]
head(fib2016vs2016_2018_sn)

fib2016vs2016_2018_sn <- na.omit(fib2016vs2016_2018_sn)

library(grid)
library(gridExtra)
grob_p = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(fib2016vs2016_2018_sn$NES, fib2016vs2016_2018_sn$NES.1,method = "pearson"), 4) ), x = 0.5, y = 0.97, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))
grob_s = grobTree(textGrob(paste("Spearman Correlation : ", round(cor(fib2016vs2016_2018_sn$NES, fib2016vs2016_2018_sn$NES.1,method = "spearman"), 4) ), x = 0.5, y = 0.93, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))



library(ggplot2)
ggplot(fib2016vs2016_2018_sn, aes(NES, NES.1))+ 
  geom_point() + 
  
  #stat_summary(fun.data=mean_cl_normal) + 
  #add a regression line
  geom_smooth(method='lm',formula=y~x,se=TRUE, fullrange=FALSE, level=0.95) +
  #stat_cor(method = "pearson", label.x = 0, label.y = 0) +
  
  annotation_custom(grob_p) +
  
  annotation_custom(grob_s) +
  labs(title = "fib2016 vs. fib2016_2018 (sn not remove) ", x="fib2016 NES", y="fib 2018 NES") +
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
#-----plot fib2018 vs fib2016+2018 sn not re-------

fib2018vs2016_2018_sn <- cbind(fib2018_log,fib2016_2018_sn[match(fib2018_log$NAME,fib2016_2018_sn$NAME),])
fib2018vs2016_2018_sn <-  fib2018vs2016_2018_sn[,c(2,4)]
head(fib2018vs2016_2018_sn)

fib2018vs2016_2018_sn <- na.omit(fib2018vs2016_2018_sn)

library(grid)
library(gridExtra)
grob_p = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(fib2018vs2016_2018_sn$NES, fib2018vs2016_2018_sn$NES.1,method = "pearson"), 4) ), x = 0.5, y = 0.97, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))
grob_s = grobTree(textGrob(paste("Spearman Correlation : ", round(cor(fib2018vs2016_2018_sn$NES, fib2018vs2016_2018_sn$NES.1,method = "spearman"), 4) ), x = 0.5, y = 0.93, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))



library(ggplot2)
ggplot(fib2018vs2016_2018_sn, aes(NES, NES.1))+ 
  geom_point() + 
  
  #stat_summary(fun.data=mean_cl_normal) + 
  #add a regression line
  geom_smooth(method='lm',formula=y~x,se=TRUE, fullrange=FALSE, level=0.95) +
  #stat_cor(method = "pearson", label.x = 0, label.y = 0) +
  
  annotation_custom(grob_p) +
  
  annotation_custom(grob_s) +
  labs(title = "fib2018 vs. fib2016_2018 (sn not remove) ", x="fib2016 NES", y="fib 2018 NES") +
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
fib2016_2018_fast_rb_log <- read.table("./GSEA/result_from_server2/fib2016+2018_rb_log_gsea_report_for_FAST_1546773678001.tsv", header = T,sep="\t")
fib2016_2018_slow_rb_log <- read.table("./GSEA/result_from_server2/fib2016+2018_rb_log_gsea_report_for_SLOW_1546773678001.tsv", header = T,sep="\t")

esc2016_2018_fast_rb_log <- read.table("./GSEA/result_from_server2/esc2016+2018_rb_log_gsea_report_for_FAST_1546776448028.tsv", header = T,sep="\t")
esc2016_2018_slow_rb_log <- read.table("./GSEA/result_from_server2/esc2016+2018_rb_log_gsea_report_for_SLOW_1546776448028.tsv", header = T,sep="\t")

fib2016_2018_fast_rb_log <-as.data.frame(fib2016_2018_fast_rb_log[,c(1,6)])
fib2016_2018_slow_rb_log <-as.data.frame(fib2016_2018_slow_rb_log[,c(1,6)])
esc2016_2018_fast_rb_log <-as.data.frame(esc2016_2018_fast_rb_log[,c(1,6)])
esc2016_2018_slow_rb_log <-as.data.frame(esc2016_2018_slow_rb_log[,c(1,6)])


fib2016_2018_rb_log <- rbind(fib2016_2018_slow_rb_log,fib2016_2018_fast_rb_log)
esc2016_2018_rb_log <- rbind(esc2016_2018_slow_rb_log,esc2016_2018_fast_rb_log)


#----plot fib2016 vs fib2016+2018 log re-----------
fib2016vs2016_2018_rb_log <- cbind(fib2016_log,fib2016_2018_rb_log[match(fib2016_log$NAME,fib2016_2018_rb_log$NAME),])
fib2016vs2016_2018_rb_log <-  fib2016vs2016_2018_rb_log[,c(2,4)]
head(fib2016vs2016_2018_rb_log)
fib2016vs2016_2018_rb_log <- na.omit(fib2016vs2016_2018_rb_log)

library(grid)
library(gridExtra)
grob_p = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(fib2016vs2016_2018_rb_log$NES, fib2016vs2016_2018_rb_log$NES.1,method = "pearson"), 4) ), x = 0.5, y = 0.97, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))
grob_s = grobTree(textGrob(paste("Spearman Correlation : ", round(cor(fib2016vs2016_2018_rb_log$NES, fib2016vs2016_2018_rb_log$NES.1,method = "spearman"), 4) ), x = 0.5, y = 0.93, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))



library(ggplot2)
ggplot(fib2016vs2016_2018_rb_log, aes(NES, NES.1))+ 
  geom_point() + 
  
  #stat_summary(fun.data=mean_cl_normal) + 
  #add a regression line
  geom_smooth(method='lm',formula=y~x,se=TRUE, fullrange=FALSE, level=0.95) +
  #stat_cor(method = "pearson", label.x = 0, label.y = 0) +
  
  annotation_custom(grob_p) +
  
  annotation_custom(grob_s) +
  labs(title = "fib2016 vs. fib2016_2018 (log remove) ", x="fib2016 NES", y="fib 2018 NES") +
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
#----plot fib2018 vs fib2016+2018 log re-----

fib2018vs2016_2018_rb_log <- cbind(fib2018_log,fib2016_2018_rb_log[match(fib2018_log$NAME,fib2016_2018_rb_log$NAME),])
fib2018vs2016_2018_rb_log <-  fib2018vs2016_2018_rb_log[,c(2,4)]
head(fib2018vs2016_2018_rb_log)
fib2018vs2016_2018_rb_log <- na.omit(fib2018vs2016_2018_rb_log)

library(grid)
library(gridExtra)
grob_p = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(fib2018vs2016_2018_rb_log$NES, fib2018vs2016_2018_rb_log$NES.1,method = "pearson"), 4) ), x = 0.5, y = 0.97, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))
grob_s = grobTree(textGrob(paste("Spearman Correlation : ", round(cor(fib2018vs2016_2018_rb_log$NES, fib2018vs2016_2018_rb_log$NES.1,method = "spearman"), 4) ), x = 0.5, y = 0.93, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))



library(ggplot2)
ggplot(fib2018vs2016_2018_rb_log, aes(NES, NES.1))+ 
  geom_point() + 
  
  #stat_summary(fun.data=mean_cl_normal) + 
  #add a regression line
  geom_smooth(method='lm',formula=y~x,se=TRUE, fullrange=FALSE, level=0.95) +
  #stat_cor(method = "pearson", label.x = 0, label.y = 0) +
  
  annotation_custom(grob_p) +
  
  annotation_custom(grob_s) +
  labs(title = "fib2018 vs. fib2016_2018 (log remove) ", x="fib2016 NES", y="fib 2018 NES") +
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
fib2016_2018_fast_rb_sn <- read.table("./GSEA/result_from_server2/fib2016+2018_rb_sn_gsea_report_for_FAST_1546773653298.tsv", header = T,sep="\t")
fib2016_2018_slow_rb_sn <- read.table("./GSEA/result_from_server2/fib2016+2018_rb_sn_gsea_report_for_SLOW_1546773653298.tsv", header = T,sep="\t")

esc2016_2018_fast_rb_sn <- read.table("./GSEA/result_from_server2/esc2016+2018_rb_sn_gsea_report_for_FAST_1546749349230.tsv", header = T,sep="\t")
esc2016_2018_slow_rb_sn <- read.table("./GSEA/result_from_server2/esc2016+2018_rb_sn_gsea_report_for_SLOW_1546749349230.tsv", header = T,sep="\t")


fib2016_2018_fast_rb_sn <-as.data.frame(fib2016_2018_fast_rb_sn[,c(1,6)])
fib2016_2018_slow_rb_sn <-as.data.frame(fib2016_2018_slow_rb_sn[,c(1,6)])
esc2016_2018_fast_rb_sn <-as.data.frame(esc2016_2018_fast_rb_sn[,c(1,6)])
esc2016_2018_slow_rb_sn <-as.data.frame(esc2016_2018_slow_rb_sn[,c(1,6)])


fib2016_2018_rb_sn <- rbind(fib2016_2018_slow_rb_sn,fib2016_2018_fast_rb_sn)
esc2016_2018_rb_sn <- rbind(esc2016_2018_slow_rb_sn,esc2016_2018_fast_rb_sn)

#----------plot fib 2016 vs fib 2016+2018 re sn-----
fib2016vs2016_2018_rb_sn <- cbind(fib2016_log,fib2016_2018_rb_sn[match(fib2016_log$NAME,fib2016_2018_rb_sn$NAME),])
fib2016vs2016_2018_rb_sn <-  fib2016vs2016_2018_rb_sn[,c(2,4)]
head(fib2016vs2016_2018_rb_sn)

fib2016vs2016_2018_rb_sn <- na.omit(fib2016vs2016_2018_rb_sn)

library(grid)
library(gridExtra)
grob_p = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(fib2016vs2016_2018_rb_sn$NES, fib2016vs2016_2018_rb_sn$NES.1,method = "pearson"), 4) ), x = 0.5, y = 0.97, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))
grob_s = grobTree(textGrob(paste("Spearman Correlation : ", round(cor(fib2016vs2016_2018_rb_sn$NES, fib2016vs2016_2018_rb_sn$NES.1,method = "spearman"), 4) ), x = 0.5, y = 0.93, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))



library(ggplot2)
ggplot(fib2016vs2016_2018_rb_sn, aes(NES, NES.1))+ 
  geom_point() + 
  
  #stat_summary(fun.data=mean_cl_normal) + 
  #add a regression line
  geom_smooth(method='lm',formula=y~x,se=TRUE, fullrange=FALSE, level=0.95) +
  #stat_cor(method = "pearson", label.x = 0, label.y = 0) +
  
  annotation_custom(grob_p) +
  
  annotation_custom(grob_s) +
  labs(title = "fib2016 vs. fib2016_2018 (sn remove) ", x="fib2016 NES", y="fib 2018 NES") +
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


#----------plot fib 2018 vs fib 2016+2018 re sn--------------
fib2018vs2016_2018_rb_sn <- cbind(fib2018_log,fib2016_2018_rb_sn[match(fib2018_log$NAME,fib2016_2018_rb_sn$NAME),])
fib2018vs2016_2018_rb_sn <-  fib2018vs2016_2018_rb_sn[,c(2,4)]
head(fib2018vs2016_2018_rb_sn)
fib2018vs2016_2018_rb_sn <- na.omit(fib2018vs2016_2018_rb_sn)

library(grid)
library(gridExtra)
grob_p = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(fib2018vs2016_2018_rb_sn$NES, fib2018vs2016_2018_rb_sn$NES.1,method = "pearson"), 4) ), x = 0.5, y = 0.97, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))
grob_s = grobTree(textGrob(paste("Spearman Correlation : ", round(cor(fib2018vs2016_2018_rb_sn$NES, fib2018vs2016_2018_rb_sn$NES.1,method = "spearman"), 4) ), x = 0.5, y = 0.93, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))



library(ggplot2)
ggplot(fib2018vs2016_2018_rb_sn, aes(NES, NES.1))+ 
  geom_point() + 
  
  #stat_summary(fun.data=mean_cl_normal) + 
  #add a regression line
  geom_smooth(method='lm',formula=y~x,se=TRUE, fullrange=FALSE, level=0.95) +
  #stat_cor(method = "pearson", label.x = 0, label.y = 0) +
  
  annotation_custom(grob_p) +
  
  annotation_custom(grob_s) +
  labs(title = "fib2018 vs. fib2016_2018 (sn remove) ", x="fib2016 NES", y="fib 2018 NES") +
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



# as we can see the correlation between fib 2016 or 2018 with fib 2016+2018(re-batch) will be lower, so comaparefib 2016 or 2018(re batch ) with 2016+2018(re-batch)

# load fib2016 and 2018 re batch
fib2016_fast_rb_log <- read.table("./GSEA/result_from_server2/fib2016_rb_negative_log_gsea_report_for_FAST_1546782673090.tsv", header = T,sep="\t")
fib2016_slow_rb_log <- read.table("./GSEA/result_from_server2/fib2016_rb_negative_log_gsea_report_for_SLOW_1546782673090.tsv", header = T,sep="\t")

fib2018_fast_rb_log <- read.table("./GSEA/result_from_server2/fib2018_rb_negative_log_gsea_report_for_FAST_1546780394888.tsv", header = T,sep="\t")
fib2018_slow_rb_log <- read.table("./GSEA/result_from_server2/fib2018_rb_negative_log_gsea_report_for_SLOW_1546780394888.tsv", header = T,sep="\t")

fib2016_fast_rb_log <-as.data.frame(fib2016_fast_rb_log[,c(1,6)])
fib2016_slow_rb_log <-as.data.frame(fib2016_slow_rb_log[,c(1,6)])
fib2018_fast_rb_log <-as.data.frame(fib2018_fast_rb_log[,c(1,6)])
fib2018_slow_rb_log <-as.data.frame(fib2018_slow_rb_log[,c(1,6)])


fib2016_rb_log <- rbind(fib2016_slow_rb_log,fib2016_fast_rb_log)
fib2018_rb_log <- rbind(fib2018_slow_rb_log,fib2018_fast_rb_log)

#----------plot fib 2016 vs 2108 log re batch--------------


fib2016vs2018_rb_log <- cbind(fib2016_rb_log,fib2018_rb_log[match(fib2016_rb_log$NAME,fib2018_rb_log$NAME),])
fib2016vs2018_rb_log <-  fib2016vs2018_rb_log[,c(2,4)]
nrow(fib2016vs2018_rb_log)
fib2016vs2018_rb_log <- na.omit(fib2016vs2018_rb_log)

library(grid)
library(gridExtra)
grob_p = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(fib2016vs2018_rb_log$NES, fib2016vs2018_rb_log$NES.1,method = "pearson"), 4) ), x = 0.5, y = 0.97, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))
grob_s = grobTree(textGrob(paste("Spearman Correlation : ", round(cor(fib2016vs2018_rb_log$NES, fib2016vs2018_rb_log$NES.1,method = "spearman"), 4) ), x = 0.5, y = 0.93, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))



library(ggplot2)
ggplot(fib2016vs2018_rb_log, aes(NES, NES.1))+ 
  geom_point() + 
  
  #stat_summary(fun.data=mean_cl_normal) + 
  #add a regression line
  geom_smooth(method='lm',formula=y~x,se=TRUE, fullrange=FALSE, level=0.95) +
  #stat_cor(method = "pearson", label.x = 0, label.y = 0) +
  
  annotation_custom(grob_p) +
  
  annotation_custom(grob_s) +
  
  
  
  labs(title = "fib2016 vs. fib2018 (log2 removed) ", x="fib2016 NES", y="fib 2018 NES") +
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







#-----------plot fib 2016 re  vs fib 2016+2018 re log batch---------

fib2016vs2016_2018_rb_log <- cbind(fib2016_rb_log,fib2016_2018_rb_log[match(fib2016_rb_log$NAME,fib2016_2018_rb_log$NAME),])
fib2016vs2016_2018_rb_log <-  fib2016vs2016_2018_rb_log[,c(2,4)]
nrow(fib2016vs2016_2018_rb_log)
fib2016vs2016_2018_rb_log <- na.omit(fib2016vs2016_2018_rb_log)

library(grid)
library(gridExtra)
grob_p = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(fib2016vs2016_2018_rb_log$NES, fib2016vs2016_2018_rb_log$NES.1,method = "pearson"), 4) ), x = 0.5, y = 0.97, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))
grob_s = grobTree(textGrob(paste("Spearman Correlation : ", round(cor(fib2016vs2016_2018_rb_log$NES, fib2016vs2016_2018_rb_log$NES.1,method = "spearman"), 4) ), x = 0.5, y = 0.93, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))



library(ggplot2)
ggplot(fib2016vs2016_2018_rb_log, aes(NES, NES.1))+ 
  geom_point() + 
  
  #stat_summary(fun.data=mean_cl_normal) + 
  #add a regression line
  geom_smooth(method='lm',formula=y~x,se=TRUE, fullrange=FALSE, level=0.95) +
  #stat_cor(method = "pearson", label.x = 0, label.y = 0) +
  
  annotation_custom(grob_p) +
  
  annotation_custom(grob_s) +
  
  
  
  labs(title = "fib2016 vs. fib2016+2018 (log2 removed) ", x="fib2016 NES", y="fib 2018 NES") +
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

#-----------plot fib 2018 re  vs fib 2016+2018 re log batch-----------

fib2018vs2016_2018_rb_log <- cbind(fib2018_rb_log,fib2016_2018_rb_log[match(fib2018_rb_log$NAME,fib2016_2018_rb_log$NAME),])
fib2018vs2016_2018_rb_log <-  fib2018vs2016_2018_rb_log[,c(2,4)]
nrow(fib2018vs2016_2018_rb_log)
fib2018vs2016_2018_rb_log <- na.omit(fib2018vs2016_2018_rb_log)

library(grid)
library(gridExtra)
grob_p = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(fib2018vs2016_2018_rb_log$NES, fib2018vs2016_2018_rb_log$NES.1,method = "pearson"), 4) ), x = 0.5, y = 0.97, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))
grob_s = grobTree(textGrob(paste("Spearman Correlation : ", round(cor(fib2018vs2016_2018_rb_log$NES, fib2018vs2016_2018_rb_log$NES.1,method = "spearman"), 4) ), x = 0.5, y = 0.93, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))



library(ggplot2)
ggplot(fib2018vs2016_2018_rb_log, aes(NES, NES.1))+ 
  geom_point() + 
  
  #stat_summary(fun.data=mean_cl_normal) + 
  #add a regression line
  geom_smooth(method='lm',formula=y~x,se=TRUE, fullrange=FALSE, level=0.95) +
  #stat_cor(method = "pearson", label.x = 0, label.y = 0) +
  
  annotation_custom(grob_p) +
  
  annotation_custom(grob_s) +
  
  
  
  labs(title = "fib2018 vs. fib2016+2018 (log2 removed) ", x="fib2016 NES", y="fib 2018 NES") +
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
fib2016_2018_fast_rb_sn <- read.table("./GSEA/result_from_server2/fib2016+2018_rb_sn_gsea_report_for_FAST_1546773653298.tsv", header = T,sep="\t")
fib2016_2018_slow_rb_sn <- read.table("./GSEA/result_from_server2/fib2016+2018_rb_sn_gsea_report_for_SLOW_1546773653298.tsv", header = T,sep="\t")

esc2016_2018_fast_rb_sn <- read.table("./GSEA/result_from_server2/esc2016+2018_rb_sn_gsea_report_for_FAST_1546749349230.tsv", header = T,sep="\t")
esc2016_2018_slow_rb_sn <- read.table("./GSEA/result_from_server2/esc2016+2018_rb_sn_gsea_report_for_SLOW_1546749349230.tsv", header = T,sep="\t")


fib2016_2018_fast_rb_sn <-as.data.frame(fib2016_2018_fast_rb_sn[,c(1,6)])
fib2016_2018_slow_rb_sn <-as.data.frame(fib2016_2018_slow_rb_sn[,c(1,6)])
esc2016_2018_fast_rb_sn <-as.data.frame(esc2016_2018_fast_rb_sn[,c(1,6)])
esc2016_2018_slow_rb_sn <-as.data.frame(esc2016_2018_slow_rb_sn[,c(1,6)])


fib2016_2018_rb_sn <- rbind(fib2016_2018_slow_rb_sn,fib2016_2018_fast_rb_sn)
esc2016_2018_rb_sn <- rbind(esc2016_2018_slow_rb_sn,esc2016_2018_fast_rb_sn)

#------------plot fib 2016 re  vs fib 2016+2018 re sn batch-----
fib2016vs2016_2018_rb_sn <- cbind(fib2016_rb_log,fib2016_2018_rb_sn[match(fib2016_rb_log$NAME,fib2016_2018_rb_sn$NAME),])
fib2016vs2016_2018_rb_sn <-  fib2016vs2016_2018_rb_sn[,c(2,4)]
nrow(fib2016vs2016_2018_rb_sn)
fib2016vs2016_2018_rb_sn <- na.omit(fib2016vs2016_2018_rb_sn)

library(grid)
library(gridExtra)
grob_p = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(fib2016vs2016_2018_rb_sn$NES, fib2016vs2016_2018_rb_sn$NES.1,method = "pearson"), 4) ), x = 0.5, y = 0.97, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))
grob_s = grobTree(textGrob(paste("Spearman Correlation : ", round(cor(fib2016vs2016_2018_rb_sn$NES, fib2016vs2016_2018_rb_sn$NES.1,method = "spearman"), 4) ), x = 0.5, y = 0.93, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))



library(ggplot2)
ggplot(fib2016vs2016_2018_rb_sn, aes(NES, NES.1))+ 
  geom_point() + 
  
  #stat_summary(fun.data=mean_cl_normal) + 
  #add a regression line
  geom_smooth(method='lm',formula=y~x,se=TRUE, fullrange=FALSE, level=0.95) +
  #stat_cor(method = "pearson", label.x = 0, label.y = 0) +
  
  annotation_custom(grob_p) +
  
  annotation_custom(grob_s) +
  
  
  
  labs(title = "fib2016 vs. fib2016+2018 (sn removed) ", x="fib2016 NES", y="fib 2018 NES") +
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

fib2018vs2016_2018_rb_sn <- cbind(fib2018_rb_log,fib2016_2018_rb_sn[match(fib2018_rb_log$NAME,fib2016_2018_rb_sn$NAME),])
fib2018vs2016_2018_rb_sn <-  fib2018vs2016_2018_rb_sn[,c(2,4)]
nrow(fib2018vs2016_2018_rb_sn)
fib2018vs2016_2018_rb_sn <- na.omit(fib2018vs2016_2018_rb_sn)

library(grid)
library(gridExtra)
grob_p = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(fib2018vs2016_2018_rb_sn$NES, fib2018vs2016_2018_rb_sn$NES.1,method = "pearson"), 4) ), x = 0.5, y = 0.97, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))
grob_s = grobTree(textGrob(paste("Spearman Correlation : ", round(cor(fib2018vs2016_2018_rb_sn$NES, fib2018vs2016_2018_rb_sn$NES.1,method = "spearman"), 4) ), x = 0.5, y = 0.93, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))



library(ggplot2)
ggplot(fib2018vs2016_2018_rb_sn, aes(NES, NES.1))+ 
  geom_point() + 
  
  #stat_summary(fun.data=mean_cl_normal) + 
  #add a regression line
  geom_smooth(method='lm',formula=y~x,se=TRUE, fullrange=FALSE, level=0.95) +
  #stat_cor(method = "pearson", label.x = 0, label.y = 0) +
  
  annotation_custom(grob_p) +
  
  annotation_custom(grob_s) +
  
  
  
  labs(title = "fib2018 vs. fib2016+2018 (log2 removed) ", x="fib2016 NES", y="fib 2018 NES") +
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


#--load fib rb sn ng------

fib2016_2018_fast_rb_sn_ng <- read.table("./GSEA/result_from_server2/fib2016+2018_rb_ng_sn_gsea_report_for_FAST_1546779207404.tsv", header = T,sep="\t")
fib2016_2018_slow_rb_sn_ng <- read.table("./GSEA/result_from_server2/fib2016+2018_rb_ng_sn_gsea_report_for_SLOW_1546779207404.tsv", header = T,sep="\t")



fib2016_2018_fast_rb_sn_ng <-as.data.frame(fib2016_2018_fast_rb_sn_ng[,c(1,6)])
fib2016_2018_slow_rb_sn_ng <-as.data.frame(fib2016_2018_slow_rb_sn_ng[,c(1,6)])


fib2016_2018_rb_sn_ng <- rbind(fib2016_2018_slow_rb_sn_ng,fib2016_2018_fast_rb_sn_ng)



# plot re sn vs re sn ng-------------------
fib_re_snvsfib_re_sn_ng <- cbind(fib2016_2018_rb_sn,fib2016_2018_rb_sn_ng[match(fib2016_2018_rb_sn$NAME,fib2016_2018_rb_sn_ng$NAME),])
fib_re_snvsfib_re_sn_ng <-  fib_re_snvsfib_re_sn_ng[,c(2,4)]
nrow(fib_re_snvsfib_re_sn_ng)
fib_re_snvsfib_re_sn_ng <- na.omit(fib_re_snvsfib_re_sn_ng)

library(grid)
library(gridExtra)
grob_p = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(fib_re_snvsfib_re_sn_ng$NES, fib_re_snvsfib_re_sn_ng$NES.1,method = "pearson"), 4) ), x = 0.5, y = 0.97, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))
grob_s = grobTree(textGrob(paste("Spearman Correlation : ", round(cor(fib_re_snvsfib_re_sn_ng$NES, fib_re_snvsfib_re_sn_ng$NES.1,method = "spearman"), 4) ), x = 0.5, y = 0.93, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))



library(ggplot2)
ggplot(fib_re_snvsfib_re_sn_ng, aes(NES, NES.1))+ 
  geom_point() + 
  
  #stat_summary(fun.data=mean_cl_normal) + 
  #add a regression line
  geom_smooth(method='lm',formula=y~x,se=TRUE, fullrange=FALSE, level=0.95) +
  #stat_cor(method = "pearson", label.x = 0, label.y = 0) +
  
  annotation_custom(grob_p) +
  
  annotation_custom(grob_s) +
  
  
  
  labs(title = "fib2016+2018 vs. fib2016+2018_ng (sn removed) ", x="fib2016 NES", y="fib 2018 NES") +
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



