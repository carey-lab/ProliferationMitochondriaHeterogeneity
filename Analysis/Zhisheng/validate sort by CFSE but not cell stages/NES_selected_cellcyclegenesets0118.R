# Jan 18 2019
# plot sox2 vs. ESC, GSE vs. ESC, sox2 vs. GSE
# only plot cell cycle gene sets

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

# get cell cycle gene sets

sox2_ESC <- sox2_ESC[c(grep("FISCHER_G.*",sox2_ESC$NAME),grep("WHITFIELD_.*",sox2_ESC$NAME)),]

sox2_ESC <- sox2_ESC[c(2,1,7,6,3,4,8),]
  
  


library(ggplot2)
ggplot(sox2_ESC, aes(NES, NES.1))+ 
  geom_point(aes(colour = sox2_ESC$NAME),size=4) + 
  
  scale_fill_discrete(name="gene set name") +
  
  guides(fill=guide_legend(title="New Legend Title")) +
                      # breaks=c("ctrl", "trt1", "trt2"),
                      # labels=c("Control", "Treatment 1", "Treatment 2")) +
  
  # geom_point(data=sox2_ESC[grep("FISCHER_G.*",sox2_ESC$NAME),], aes(NES, NES.1), colour="blue", size=3) +
  # 
  # geom_point(data=sox2_ESC[grep("WHITFIELD_.*",sox2_ESC$NAME),], aes(NES, NES.1), colour="green", size=3) +
  # 
  # scale_color_gradient(low = "blue", high = "yellow") +
  
  
  labs(title = "sox2 vs. ESCs",colour="gene set name") +
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
  axis.title.y = element_text(color="blue", size=14, face="bold"),
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
  legend.title = element_text(colour="black", size=14, face="bold"),
  
  legend.text = element_text(colour="black", size=10, face="bold")
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
GSE_ESC <-  GSE_ESC[,c(1,2,4)]
head(GSE_ESC)

GSE_ESC <- na.omit(GSE_ESC)
nrow(GSE_ESC)
GSE_ESC$NES.1 <- GSE_ESC$NES.1*(-1)




# get cell cycle gene sets

GSE_ESC <- GSE_ESC[c(grep("FISCHER_G.*",GSE_ESC$NAME),grep("WHITFIELD_.*",GSE_ESC$NAME)),]

GSE_ESC <- GSE_ESC[-3,]




library(ggplot2)
ggplot(GSE_ESC, aes(NES, NES.1))+ 
  geom_point(aes(colour = GSE_ESC$NAME),size=4) + 
  
  scale_fill_discrete(name="gene set name") +
  
  guides(fill=guide_legend(title="New Legend Title")) +
  # breaks=c("ctrl", "trt1", "trt2"),
  # labels=c("Control", "Treatment 1", "Treatment 2")) +
  
  # geom_point(data=GSE_ESC[grep("FISCHER_G.*",GSE_ESC$NAME),], aes(NES, NES.1), colour="blue", size=3) +
  # 
  # geom_point(data=GSE_ESC[grep("WHITFIELD_.*",GSE_ESC$NAME),], aes(NES, NES.1), colour="green", size=3) +
  # 
  # scale_color_gradient(low = "blue", high = "yellow") +
  
  
  labs(title = "GSE vs. ESCs",colour="gene set name") +
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
  axis.title.y = element_text(color="blue", size=14, face="bold"),
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
  legend.title = element_text(colour="black", size=14, face="bold"),
  
  legend.text = element_text(colour="black", size=10, face="bold")
)








#-----plot GSE vs esc2016+2018 sn ng re-------
GSE_sox2 <- cbind(GSE,sox2[match(GSE$NAME,sox2$NAME),])
GSE_sox2 <-  GSE_sox2[,c(1,2,4)]
head(GSE_sox2)

GSE_sox2 <- na.omit(GSE_sox2)
nrow(GSE_sox2)

# get cell cycle gene sets

GSE_sox2 <- GSE_sox2[c(grep("FISCHER_G.*",GSE_sox2$NAME),grep("WHITFIELD_.*",GSE_sox2$NAME)),]

GSE_sox2 <- GSE_sox2[-3,]




library(ggplot2)
ggplot(GSE_sox2, aes(NES, NES.1))+ 
  geom_point(aes(colour = GSE_sox2$NAME),size=4) + 
  
  scale_fill_discrete(name="gene set name") +
  
  guides(fill=guide_legend(title="New Legend Title")) +
  # breaks=c("ctrl", "trt1", "trt2"),
  # labels=c("Control", "Treatment 1", "Treatment 2")) +
  
  # geom_point(data=GSE_sox2[grep("FISCHER_G.*",GSE_sox2$NAME),], aes(NES, NES.1), colour="blue", size=3) +
  # 
  # geom_point(data=GSE_sox2[grep("WHITFIELD_.*",GSE_sox2$NAME),], aes(NES, NES.1), colour="green", size=3) +
  # 
  # scale_color_gradient(low = "blue", high = "yellow") +
  
  
  labs(title = "GSE vs. sox2",colour="gene set name") +
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
  axis.title.y = element_text(color="blue", size=14, face="bold"),
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
  legend.title = element_text(colour="black", size=14, face="bold"),
  
  legend.text = element_text(colour="black", size=10, face="bold")
)


