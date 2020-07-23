# Jan 18

# cell cycle gene sets

# Jan 17

# get the yeast GSEA result-----

# heatmap of GSEA cell cycle

NES_FIB_ESCs_yeast <- read.table("./NES_FIB_ESCs_yeast.tsv",header = T,sep = "\t")
head(NES_FIB_ESCs_yeast)

all_cell_cycle <- NES_FIB_ESCs_yeast[c(grep("FISCHER_G.*",NES_FIB_ESCs_yeast$NAME.FIB), grep("WHITFIELD_.*",NES_FIB_ESCs_yeast$NAME.FIB)),]

all_cell_cycle <- all_cell_cycle[-3,]

all_cell_cycle <- all_cell_cycle[,c(1,4,9,17)]

all_cell_cycle <- all_cell_cycle[c(2,1,6,4,3,5,7),]

library(ggplot2)
library(reshape2)

all_cell_cycle.m <- melt(all_cell_cycle)

all_cell_cycle.m$value <- all_cell_cycle.m$value*-1


# ggplot(all_cell_cycle.m,aes(variable,NAME.FIB)) +
#   geom_tile(aes(fill = all_cell_cycle.m$value, colour = "white") +
#               scale_fill_gradient(low = "white",high = "steelblue")) +
  
            
ggplot(all_cell_cycle.m,aes(variable,NAME.FIB,fill=value)) +
  geom_tile() +
  geom_text(aes(label = round(value, 1))) +
  scale_fill_gradient(low = "white",high = "steelblue") +
  labs(title = "NES of cell cycle gene sets",x="",y="") +
  theme(
    plot.title = element_text(lineheight=.8, size=20,face="bold",hjust = 0.5),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    # axis.title.y = element_blank(),
    # axis.title.x = element_blank(),
    #delete background
    panel.grid.major =element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
   
    #x轴标签
    axis.text.x = element_text(size=10,color='black',face="bold"),
    #y轴标签
    axis.text.y = element_text(size=10,color='black',face="bold"),
    #图例
    legend.title = element_text(colour="black", size=14, face="bold")
  )



#---------------------------------------




library(grid)
library(gridExtra)
grob_p = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(sox2_ESC$NES, sox2_ESC$NES.1,method = "pearson"), 4) ), x = 0.5, y = 0.97, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))
grob_s = grobTree(textGrob(paste("Spearman Correlation : ", round(cor(sox2_ESC$NES, sox2_ESC$NES.1,method = "spearman"), 4) ), x = 0.5, y = 0.93, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))



library(ggplot2)
ggplot(sox2_ESC, aes(NES, NES.1))+ 
  geom_point() + 
  
  geom_point(data=sox2_ESC[sox2_ESC$NAME=="HALLMARK_MYC_TARGETS_V1",], aes(NES, NES.1), colour="red", size=3) +
  
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
