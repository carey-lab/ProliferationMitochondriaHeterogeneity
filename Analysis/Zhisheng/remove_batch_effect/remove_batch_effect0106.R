# only remove batch effect with sv and combat

library(sva)
# library(bladderbatch)
# data(bladderdata)
library(pamr)
library(limma)
require(Biobase)
#--------------------fib all--------------------------------
#set up the data
#expression data
fiball <- read.table("./GSEA/FIB_GSEA_1231.gct",header = T)
head(fiball)


rownames(fiball) <- fiball$NAME
fiball <- fiball[,c(3,4,5,6,7,8,9,10)]
fiball <- as.matrix(fiball)
rownames(fiball)


#phenotype
fib_pheno <- data.frame("batch"=rep(c(1,1,2,2),2),"celltype"=c(rep("slow",4),c(rep("fast",4))))

rownames(fib_pheno) <- colnames(fiball)

#the interested variable
mod = model.matrix(~as.factor(celltype), data=fib_pheno)
head(mod)
#adjusted variables
mod0 = model.matrix(~1,data=fib_pheno)
head(mod0)

#the number of latent factors that need to be estimated
n.sv = num.sv(fiball,mod,method="leek")
n.sv

# estimate the surrogate variables
# svobj = sva(fiball,mod,mod0)

#remove all 0
fiball2 <- fiball[rowMax(fiball)>0,]
nrow(fiball2) # 16849
# # add [0.99,1.01]
# 
# runif(8, min=0.99, max=1.01)
# apply(fiball, 1, var)!=0

# add random [0.99,1.01] to rowmin=0
fiball2[rowMin(fiball2)==0,] <- fiball2[rowMin(fiball2)==0,]+runif(8, min=0.99, max=1.01)


svobj2 = sva(fiball2,mod,mod0)

#Adjusting for surrogate variables using the f.pvalue function
pValues = f.pvalue(fiball2,mod,mod0)
qValues = p.adjust(pValues,method="BH")

#Applying the ComBat function to adjust for known batches
batch = fib_pheno$batch
modcombat = model.matrix(~1, data=fib_pheno)

#get NAN with 
combat_fiball2 = ComBat(dat=fiball2, batch=batch, mod=modcombat, par.prior=T, prior.plots= F)

head(fiball2)
head(combat_fiball2) #16849
nrow(fiball2)

pValuesComBat = f.pvalue(combat_fiball2,mod,mod0)
qValuesComBat = p.adjust(pValuesComBat,method="BH")

length(qValues)
length(qValues[qValues<0.05])
length(qValuesComBat[qValuesComBat<0.05])

#---plot PCA compare result-----------------------------

library(ggplot2)
library(ggfortify)

tfiball2 <- as.data.frame(t(fiball2)) 
tfiball2$name <- rownames(tfiball2)

tcombat_fiball2 <- as.data.frame(t(combat_fiball2))
tcombat_fiball2$name <- rownames(tcombat_fiball2)

autoplot(prcomp(tfiball2[,-16850]),data = tfiball2, colour = 'name',label = TRUE, label.size = 3)

autoplot(prcomp(tcombat_fiball2[,-16850]),data = tcombat_fiball2, colour = 'name',label = TRUE, label.size = 3)

#---plot rlog transformed PCA----

library(DESeq2)
library(gplots)
library(RColorBrewer)
library(ggplot2)

coldata <- read.table(file="./PCA/fib gene id counts coldata.csv", header=T,sep = ',')
head(coldata)

fiball2 <- round(fiball2)
nrow(fiball2) #48177
head(fiball2)
dds <- DESeqDataSetFromMatrix(fiball2,colData=coldata,design=~condition)

dds<-DESeq(dds)
head(counts(dds))

dds_rLogTransformed <- rlog( dds ) # regularized-log transform the data , for PCA & clustering

# apply PCA
plotPCA( dds_rLogTransformed , intgroup = c( "condition","year") ) 


# ESC_slow1 ESC_slow2 ESC_slow3 ESC_slow4 ESC_fast1 ESC_fast2 ESC_fast3 ESC_fast4
# RPS27         28       430        14        19       437       407        18        27
# TMEM72       156      1155        13        13       111       246         7         7


# save remove batch effect fib with negative values-------

fib_rb <- as.data.frame(combat_fiball2) 
fib_rb$NAME <- rownames(fib_rb)
fib_rb$Description <- rep("na",nrow(fib_rb))
head(fib_rb)



## all
fib_rb <- fib_rb[,c(9,10,1:8)]
colSums(fib_rb[,c(3,4,5,6,7,8,9,10)])
nrow(fib_rb)


write.table(fib_rb,file="./GSEA/fib2016+2018_rb_negative.gct",sep = "\t",row.names = F,col.names = T,quote = F)

# 0106 fib_rb_2016
fib_rb_2016 <-fib_rb[,c(9,10,1,2,5,6)]
fib_rb_2018 <-fib_rb[,c(9,10,3,4,7,8)]
write.table(fib_rb_2016,file="./GSEA/fib2016_rb_negative.gct",sep = "\t",row.names = F,col.names = T,quote = F)
write.table(fib_rb_2018,file="./GSEA/fib2018_rb_negative.gct",sep = "\t",row.names = F,col.names = T,quote = F)



# remove abnormal genes
temp <- combat_fiball2[rowMin(combat_fiball2)<0,]
#         FIB_slow1 FIB_slow2 FIB_slow3 FIB_slow4 FIB_fast1 FIB_fast2 FIB_fast3 FIB_fast4
# RPS27       -29       135        86        89       171       176        95        89
# PTPRN        13        13        -1        -1        18        19        41        43
# ELN         286       237       156       158         1        -9        65        66

combat_fiball2 <- combat_fiball2[ !(c(rownames(combat_fiball2) %in% c("RPS27","PTPRN","ELN"))),]

# save remove batch effect fib-----
fib_rb <- as.data.frame(combat_fiball2) 
fib_rb$NAME <- rownames(fib_rb)
fib_rb$Description <- rep("na",nrow(fib_rb))
head(fib_rb)

fib_rb <- fib_rb[,c(9,10,1:8)]
colSums(fib_rb[,c(3,4,5,6,7,8,9,10)])
nrow(fib_rb)


write.table(fib_rb,file="./GSEA/fib2016+2018_rb.gct",sep = "\t",row.names = F,col.names = T,quote = F)





# rlog PCA after removing batch effect

combat_fiball2 <- round(combat_fiball2)


dds <- DESeqDataSetFromMatrix(combat_fiball2,colData=coldata,design=~condition)

dds<-DESeq(dds)
head(counts(dds))

dds_rLogTransformed <- rlog( dds ) # regularized-log transform the data , for PCA & clustering

# apply PCA
plotPCA( dds_rLogTransformed , intgroup = c( "condition","year") ) 




#-----------------looks good---to be continue---------







#--------------------esc all--------------------------------
#set up the data
#expression data
escall <- read.table("./GSEA/ESC_GSEA_0104.gct",header = T)
head(escall)


rownames(escall) <- escall$NAME
escall <- escall[,c(3,4,5,6,7,8,9,10)]
escall <- as.matrix(escall)
rownames(escall)


#phenotype
esc_pheno <- data.frame("batch"=rep(c(1,1,2,2),2),"celltype"=c(rep("slow",4),c(rep("fast",4))))

rownames(esc_pheno) <- colnames(escall)

#the interested variable
mod = model.matrix(~as.factor(celltype), data=esc_pheno)
head(mod)
#adjusted variables
mod0 = model.matrix(~1,data=esc_pheno)
head(mod0)

#the number of latent factors that need to be estimated
n.sv = num.sv(escall,mod,method="leek")
n.sv

# estimate the surrogate variables
# svobj = sva(escall,mod,mod0)

#try if remove 0 values , it works , I do not know why
#also remove n.sv
# need to add 0 later




#remove all 0
escall2 <- escall[rowMax(escall)>0,]
nrow(escall2) # 17765
# # add [0.99,1.01]
# 
# runif(8, min=0.99, max=1.01)
# apply(escall, 1, var)!=0

# add random [0.99,1.01] to rowmin=0
escall2[rowMin(escall2)==0,] <- escall2[rowMin(escall2)==0,]+runif(8, min=0.99, max=1.01)


svobj2 = sva(escall2,mod,mod0)

#Adjusting for surrogate variables using the f.pvalue function
pValues = f.pvalue(escall2,mod,mod0)
qValues = p.adjust(pValues,method="BH")

#Applying the ComBat function to adjust for known batches
batch = esc_pheno$batch
modcombat = model.matrix(~1, data=esc_pheno)

#get NAN with 
combat_escall2 = ComBat(dat=escall2, batch=batch, mod=modcombat, par.prior=T, prior.plots= F)

head(escall2)
head(combat_escall2) #17765
nrow(escall2)

pValuesComBat = f.pvalue(combat_escall2,mod,mod0)
qValuesComBat = p.adjust(pValuesComBat,method="BH")

length(qValues)
length(qValues[qValues<0.05])
length(qValuesComBat[qValuesComBat<0.05])

#---plot PCA compare result-----------------------------

library(ggplot2)
library(ggfortify)

tescall2 <- as.data.frame(t(escall2)) 
tescall2$name <- rownames(tescall2)

tcombat_escall2 <- as.data.frame(t(combat_escall2))
tcombat_escall2$name <- rownames(tcombat_escall2)

autoplot(prcomp(tescall2[,-17766]),data = tescall2, colour = 'name',label = TRUE, label.size = 3)

autoplot(prcomp(tcombat_escall2[,-17766]),data = tcombat_escall2, colour = 'name',label = TRUE, label.size = 3)

#---plot rlog transformed PCA----

library(DESeq2)
library(gplots)
library(RColorBrewer)
library(ggplot2)

coldata <- read.table(file="./PCA/esc gene id counts coldata.csv", header=T,sep = ',')
head(coldata)


escall2 <- round(escall2)
nrow(escall2) #48177
head(escall2)
dds <- DESeqDataSetFromMatrix(escall2,colData=coldata,design=~condition)

dds<-DESeq(dds)
head(counts(dds))

dds_rLogTransformed <- rlog( dds ) # regularized-log transform the data , for PCA & clustering

# apply PCA
plotPCA( dds_rLogTransformed , intgroup = c( "condition","year") ) 


# ESC_slow1 ESC_slow2 ESC_slow3 ESC_slow4 ESC_fast1 ESC_fast2 ESC_fast3 ESC_fast4
# RPS27         28       430        14        19       437       407        18        27
# TMEM72       156      1155        13        13       111       246         7         7

# save remove batch effect esc with negative values genes
esc_rb <- as.data.frame(combat_escall2) 
esc_rb$NAME <- rownames(esc_rb)
esc_rb$Description <- rep("na",nrow(esc_rb))
head(esc_rb)



esc_rb <- esc_rb[,c(9,10,1:8)]
colSums(esc_rb[,c(3,4,5,6,7,8,9,10)])
nrow(esc_rb) #17765


write.table(esc_rb,file="./GSEA/esc2016+2018_rb_negative.gct",sep = "\t",row.names = F,col.names = T,quote = F)

# esc 2016 with negative values genes
esc_rb2016 <- esc_rb[,c(9,10,1,2,5,6)]
esc_rb2018 <- esc_rb[,c(9,10,3,4,7,8)]

write.table(esc_rb2016,file="./GSEA/esc2016_rb_negative.gct",sep = "\t",row.names = F,col.names = T,quote = F)

write.table(esc_rb2018,file="./GSEA/esc2018_rb_negative.gct",sep = "\t",row.names = F,col.names = T,quote = F)



# remove abnormal genes
temp <- combat_escall2[rowMin(combat_escall2)<0,]

# RPS27   -42.13915702 272.2710356 153.5637155 159.2479186 277.79072931 254.2166071 157.0885473 167.4066880
# TLX2     -0.04157660   0.5291253   0.2642053   0.3294612  -0.04324054   0.9773697   0.3895477   0.2632913
# ASB11     0.36456946   0.5591347   0.5778927   0.2553075  -0.01442152   0.3932309   0.1962083   0.1588706
# CACNA1S  -0.02913937   3.9632941   1.2076605   1.0997671  -0.07006414   0.2848215   0.8010982   0.5900857
# TMEM72   32.79124938 815.5192207 199.4185382 198.8441270  -2.29449504 103.2662812 192.0283019 192.5665615
# CD4      -0.18802488  14.9771348   8.7123860   5.9881610  -0.39285449  15.6987297   7.1008790   5.0335930
# TMEM92    0.37572733   4.7040205   2.2417659   1.2267891  -0.10986403   0.5670041   0.8916038   0.6018619
# CCSER1    0.20001383   3.4916845   0.7195113   1.1897499  -0.07875113   0.6036809   1.0124338   0.7603013
combat_escall2 <- combat_escall2[ !(c(rownames(combat_escall2) %in% c("RPS27","TLX2","ASB11","CACNA1S","TMEM72","CD4","TMEM92","CCSER1"))),]

# save remove batch effect esc
esc_rb <- as.data.frame(combat_escall2) 
esc_rb$NAME <- rownames(esc_rb)
esc_rb$Description <- rep("na",nrow(esc_rb))
head(esc_rb)

esc_rb <- esc_rb[,c(9,10,1:8)]
colSums(esc_rb[,c(3,4,5,6,7,8,9,10)])
nrow(esc_rb) #17757


write.table(esc_rb,file="./GSEA/esc2016+2018_rb.gct",sep = "\t",row.names = F,col.names = T,quote = F)




# rlog PCA after removing batch effect

combat_escall2 <- round(combat_escall2)


dds <- DESeqDataSetFromMatrix(combat_escall2,colData=coldata,design=~condition)

dds<-DESeq(dds)
head(counts(dds))

dds_rLogTransformed <- rlog( dds ) # regularized-log transform the data , for PCA & clustering

# apply PCA
plotPCA( dds_rLogTransformed , intgroup = c( "condition","year") ) 














