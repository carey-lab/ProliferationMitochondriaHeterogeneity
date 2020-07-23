# plot yeast vs. fibroblasts ,fibroblasts vs. escs, escs vs. yeast
# with FDR less than 0.05



# Jan 20

# all yeast GSEA vs. FIB  FDR less than 0.05 ------

NES_FIB_ESCs_yeast <- read.table("../../NES_FIB_ESCs_yeast.tsv",header = T,sep = "\t")
head(NES_FIB_ESCs_yeast)
colnames(NES_FIB_ESCs_yeast)
yeast_fib <- NES_FIB_ESCs_yeast[,c(1,2,4,6,17,19)]
head(yeast_fib)

nrow(yeast_fib)
yeast_fib <- na.omit(yeast_fib)

yeast_fib <- yeast_fib[! (yeast_fib$NES_FIB == 0 | yeast_fib$NES_yeast ==0),]

yeast_fib <- yeast_fib[yeast_fib$FDR.q.val <= 0.05 & yeast_fib$FDR.q.val.2 <= 0.05, ]

yeast_fib$NES_FIB <- yeast_fib$NES_FIB*-1

yeast_fib$NES_yeast <- yeast_fib$NES_yeast*-1

nrow(yeast_fib) #182

# get the same in yeast and fib

yeast_fib <- yeast_fib[yeast_fib$NES_FIB*yeast_fib$NES_yeast>0,]

nrow(yeast_fib) #167

write.table(yeast_fib,file="yeast_fib_same.tsv",col.names = T,row.names = F,sep = "\t",quote = F)
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
ESCs_yeast <- NES_FIB_ESCs_yeast[,c(1,2,9,11,17,19)]
head(ESCs_yeast)

nrow(ESCs_yeast)
ESCs_yeast <- na.omit(ESCs_yeast)

ESCs_yeast <- ESCs_yeast[! (ESCs_yeast$NES_ESC == 0 | ESCs_yeast$NES_yeast ==0),]


ESCs_yeast <- ESCs_yeast[ESCs_yeast$FDR.q.val.1 <= 0.05 & ESCs_yeast$FDR.q.val.2 <= 0.05, ]

ESCs_yeast$NES_ESC <- ESCs_yeast$NES_ESC*-1

ESCs_yeast$NES_yeast <- ESCs_yeast$NES_yeast*-1

nrow(ESCs_yeast) #40

# get the same in ESCs and yeast

ESCs_yeast <- ESCs_yeast[ESCs_yeast$NES_ESC*ESCs_yeast$NES_yeast>0,]

nrow(ESCs_yeast) #32



write.table(ESCs_yeast,file="ESCs_yeast_same.tsv",col.names = T,row.names = F,sep = "\t",quote = F)


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




############



# all FIB GSEA vs. ESCs------
colnames(NES_FIB_ESCs_yeast)
fib_ESCs <- NES_FIB_ESCs_yeast[,c(1,2,4,6,9,11)]
head(fib_ESCs)

nrow(fib_ESCs)
fib_ESCs <- na.omit(fib_ESCs)

# row_sub = apply(fib_ESCs, 1, function(row) all(row !=0 ))
# ##Subset as usual
# fib_ESCs <- fib_ESCs[row_sub,]

fib_ESCs <- fib_ESCs[! (fib_ESCs$NES_FIB == 0 | fib_ESCs$NES_ESC ==0),]


fib_ESCs <- fib_ESCs[fib_ESCs$FDR.q.val <= 0.05 & fib_ESCs$FDR.q.val.1 <= 0.05, ]

fib_ESCs$NES_FIB <- fib_ESCs$NES_FIB*-1

fib_ESCs$NES_ESC <- fib_ESCs$NES_ESC*-1

nrow(fib_ESCs) #127

fib_ESCs <- fib_ESCs[fib_ESCs$NES_ESC*fib_ESCs$NES_FIB>0,]

nrow(fib_ESCs) #106



write.table(fib_ESCs,file="fib_ESCs_same.tsv",col.names = T,row.names = F,sep = "\t",quote = F)



#plot-------------------

library(grid)
library(gridExtra)
grob_p = grobTree(textGrob(paste("Pearson Correlation : ", round(cor(fib_ESCs$NES_FIB, fib_ESCs$NES_ESC,method = "pearson"), 4) ), x = 0.5, y = 0.97, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))
grob_s = grobTree(textGrob(paste("Spearman Correlation : ", round(cor(fib_ESCs$NES_FIB, fib_ESCs$NES_ESC,method = "spearman"), 4) ), x = 0.5, y = 0.93, hjust = 0, gp = gpar(col = "red", fontsize = 11, fontface = "bold")))



library(ggplot2)
ggplot(fib_ESCs, aes(NES_FIB, NES_ESC))+ 
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



# overlap of three 

colnames(NES_FIB_ESCs_yeast)
yeast_fib_ESCs <- NES_FIB_ESCs_yeast[,c(1,2,4,6,9,11,17,19)]

head(yeast_fib_ESCs)

nrow(yeast_fib_ESCs)
yeast_fib_ESCs <- na.omit(yeast_fib_ESCs)

# row_sub = apply(yeast_fib_ESCs, 1, function(row) all(row !=0 ))
# ##Subset as usual
# yeast_fib_ESCs <- yeast_fib_ESCs[row_sub,]

yeast_fib_ESCs <- yeast_fib_ESCs[! (yeast_fib_ESCs$NES_FIB == 0 | yeast_fib_ESCs$NES_ESC ==0 | yeast_fib_ESCs$NES_yeast == 0),]


yeast_fib_ESCs <- yeast_fib_ESCs[yeast_fib_ESCs$FDR.q.val <= 0.05 & yeast_fib_ESCs$FDR.q.val.1 <= 0.05 & yeast_fib_ESCs$FDR.q.val.2<=0.05, ]

yeast_fib_ESCs[,c(3,5,7)] <- yeast_fib_ESCs[,c(3,5,7)]*-1

yeast_fib_ESCs <- yeast_fib_ESCs[yeast_fib_ESCs$NES_ESC >0 & yeast_fib_ESCs$NES_FIB >0 & yeast_fib_ESCs$NES_yeast>0,]

nrow(yeast_fib_ESCs) 


write.table(yeast_fib_ESCs,file="yeast_fib_ESCs_same.tsv",col.names = T,row.names = F,sep = "\t",quote = F)


library(scatterplot3d)
scatterplot3d(yeast_fib_ESCs[,c(3,5,7)],
              main="yeast_fib_ESCs",
              xlab = "NES_FIB",
              ylab = "NES_ESC",
              zlab = "NES_yeast",
              pch = 16,
              color="steelblue")



# different in mammalian cells and yeast


fib_ESCs_s_yeast_d <- NES_FIB_ESCs_yeast[,c(1,2,4,6,9,11,17,19)]
head(fib_ESCs_s_yeast_d)

nrow(fib_ESCs_s_yeast_d)
# fib_ESCs_s_yeast_d <- na.omit(fib_ESCs_s_yeast_d)

# row_sub = apply(fib_ESCs, 1, function(row) all(row !=0 ))
# ##Subset as usual
# fib_ESCs <- fib_ESCs[row_sub,]

fib_ESCs_s_yeast_d <- fib_ESCs_s_yeast_d[! (fib_ESCs_s_yeast_d$NES_FIB == 0 | fib_ESCs_s_yeast_d$NES_ESC ==0),]


fib_ESCs_s_yeast_d <- fib_ESCs_s_yeast_d[fib_ESCs_s_yeast_d$FDR.q.val <= 0.05 & fib_ESCs_s_yeast_d$FDR.q.val.1 <= 0.05, ]

fib_ESCs_s_yeast_d$NES_FIB <- fib_ESCs_s_yeast_d$NES_FIB*-1

fib_ESCs_s_yeast_d$NES_ESC <- fib_ESCs_s_yeast_d$NES_ESC*-1

fib_ESCs_s_yeast_d$NES_yeast <- fib_ESCs_s_yeast_d$NES_yeast*-1

nrow(fib_ESCs_s_yeast_d) #127

fib_ESCs_s_yeast_d <- fib_ESCs_s_yeast_d[fib_ESCs_s_yeast_d$NES_ESC*fib_ESCs_s_yeast_d$NES_FIB>0,]

nrow(fib_ESCs_s_yeast_d) #106



# differ in yeast

temp <- fib_ESCs_s_yeast_d$NES_FIB*fib_ESCs_s_yeast_d$NES_yeast<0

temp[is.na(temp)] <- TRUE

fib_ESCs_s_yeast_d <- fib_ESCs_s_yeast_d[temp,]

nrow(fib_ESCs_s_yeast_d)

#fib_ESCs_s_yeast_d <- fib_ESCs_s_yeast_d[fib_ESCs_s_yeast_d$FDR.q.val.2<0.05,]


write.table(fib_ESCs_s_yeast_d,file="fib_ESCs_same_yeast_diff.tsv",col.names = T,row.names = F,sep = "\t",quote = F)


fib_ESCs_s_yeast_d <- na.omit(fib_ESCs_s_yeast_d)


write.table(fib_ESCs_s_yeast_d,file="fib_ESCs_same_yeast_diff_naremove.tsv",col.names = T,row.names = F,sep = "\t",quote = F)


# diff in fib and escs-----

fib_ESCs_diff <- NES_FIB_ESCs_yeast[,c(1,2,4,6,9,11)]
head(fib_ESCs_diff)

nrow(fib_ESCs_diff)
fib_ESCs_diff <- na.omit(fib_ESCs_diff)

# row_sub = apply(fib_ESCs_diff, 1, function(row) all(row !=0 ))
# ##Subset as usual
# fib_ESCs_diff <- fib_ESCs_diff[row_sub,]

fib_ESCs_diff <- fib_ESCs_diff[! (fib_ESCs_diff$NES_FIB == 0 | fib_ESCs_diff$NES_ESC ==0),]


fib_ESCs_diff <- fib_ESCs_diff[fib_ESCs_diff$FDR.q.val <= 0.05 & fib_ESCs_diff$FDR.q.val.1 <= 0.05, ]

fib_ESCs_diff$NES_FIB <- fib_ESCs_diff$NES_FIB*-1

fib_ESCs_diff$NES_ESC <- fib_ESCs_diff$NES_ESC*-1

nrow(fib_ESCs_diff) #127

fib_ESCs_diff <- fib_ESCs_diff[fib_ESCs_diff$NES_ESC*fib_ESCs_diff$NES_FIB<0,]

nrow(fib_ESCs_diff) #21



write.table(fib_ESCs_diff,file="fib_ESCs_diff.tsv",col.names = T,row.names = F,sep = "\t",quote = F)





