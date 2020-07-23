# only remove batch effect with sv and combat

library(sva)
library(bladderbatch)
data(bladderdata)
library(pamr)
library(limma)

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

#try if remove 0 values , it works , I do not know why
#also remove n.sv
# need to add 0 later
fiball2 <- fiball[rowMin(fiball)>0,]
fiball0 <- fiball[rowMin(fiball)<=0,]
nrow(fiball2)
svobj2 = sva(fiball2,mod,mod0)

#Adjusting for surrogate variables using the f.pvalue function
pValues = f.pvalue(fiball2,mod,mod0)
qValues = p.adjust(pValues,method="BH")

#Applying the ComBat function to adjust for known batches
batch = fib_pheno$batch
modcombat = model.matrix(~1, data=fib_pheno)


combat_fiball2 = ComBat(dat=fiball2, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots= F)

head(fiball2)
head(combat_fiball2)
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

autoplot(prcomp(tfiball2[,-13180]),data = tfiball2, colour = 'name',label = TRUE, label.size = 3)

autoplot(prcomp(tcombat_fiball2[,-13180]),data = tcombat_fiball2, colour = 'name',label = TRUE, label.size = 3)

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


# convert negative to 1
combat_fiball2 <- round(combat_fiball2)
combat_fiball2[combat_fiball2<0] <- 1
                                                                   

dds <- DESeqDataSetFromMatrix(combat_fiball2,colData=coldata,design=~condition)

dds<-DESeq(dds)
head(counts(dds))
head(counts(dds_rLogTransformed))

dds_rLogTransformed <- rlog( dds ) # regularized-log transform the data , for PCA & clustering

# apply PCA
plotPCA( dds_rLogTransformed , intgroup = c( "condition","year") ) 








# add removed 0 values

fib_rb <- rbind(combat_fiball2,fiball0)
head(fib_rb)


# transfer to .gct file and save
fib_rb <- as.data.frame(fib_rb)
fib_rb$NAME <- rownames(fib_rb)
fib_rb$Description <- rep("na",nrow(fib_rb))
head(fib_rb)

fib_rb <- fib_rb[,c(9,10,1,2,5,6,3,4,7,8)]
colSums(fib_rb[,c(3,4,5,6,7,8,9,10)])
nrow(fib_rb)
fib_rb[fib_rb<0] <- 0

write.table(fib_rb,file="./GSEA/fib2016+2018_rb.gct",sep = "\t",row.names = F,col.names = T,quote = F)

#plot fib rb
fib_rb <- round(fib_rb[,c(3:10)])
fib_rb[fib_rb<0] <- 1

nrow(fib_rb)
dds <- DESeqDataSetFromMatrix(fib_rb,colData=coldata,design=~condition)

dds<-DESeq(dds)
head(counts(dds))
head(counts(dds_rLogTransformed))

dds_rLogTransformed <- rlog( dds ) # regularized-log transform the data , for PCA & clustering

# apply PCA
plotPCA( dds_rLogTransformed , intgroup = c( "condition","year") ) 






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
escall2 <- escall[rowMin(escall)>0,]
escall0 <- escall[rowMin(escall)<=0,]
nrow(escall2)
svobj2 = sva(escall2,mod,mod0)

#Adjusting for surrogate variables using the f.pvalue function
pValues = f.pvalue(escall2,mod,mod0)
qValues = p.adjust(pValues,method="BH")

#Applying the ComBat function to adjust for known batches
batch = esc_pheno$batch
modcombat = model.matrix(~1, data=esc_pheno)


combat_escall2 = ComBat(dat=escall2, batch=batch, mod=modcombat, par.prior=TRUE, prior.plots= F)

head(escall2)
head(combat_escall2)
nrow(escall2)

pValuesComBat = f.pvalue(combat_escall2,mod,mod0)
qValuesComBat = p.adjust(pValuesComBat,method="BH")

length(qValues)
length(qValues[qValues<0.05])
length(qValuesComBat[qValuesComBat<0.05])

#---plot PCA compare result----------------------------

library(ggplot2)
library(ggfortify)

tescall2 <- as.data.frame(t(escall2)) 
tescall2$name <- rownames(tescall2)

tcombat_escall2 <- as.data.frame(t(combat_escall2))
tcombat_escall2$name <- rownames(tcombat_escall2)

autoplot(prcomp(tescall2[,-14613]),data = tescall2, colour = 'name',label = TRUE, label.size = 3)

autoplot(prcomp(tcombat_escall2[,-14613]),data = tcombat_escall2, colour = 'name',label = TRUE, label.size = 3)


#----plot rlog transformed PCA--------------------------------------

library(DESeq2)
library(gplots)
library(RColorBrewer)
library(ggplot2)

coldata <- read.table(file="./PCA/ESC gene id counts coldata.csv", header=T,sep = ',')
head(coldata)

escall2 <- round(escall2)
nrow(escall2) #48177
head(escall2)
dds <- DESeqDataSetFromMatrix(escall2,colData=coldata,design=~condition)

dds<-DESeq(dds)
head(counts(dds))
head(counts(dds_rLogTransformed))

dds_rLogTransformed <- rlog( dds ) # regularized-log transform the data , for PCA & clustering

# apply PCA
plotPCA( dds_rLogTransformed , intgroup = c( "condition","year") ) 

escall2[c("RPS27","TMEM72"),]
# ESC_slow1 ESC_slow2 ESC_slow3 ESC_slow4 ESC_fast1 ESC_fast2 ESC_fast3 ESC_fast4
# RPS27         28       430        14        19       437       407        18        27
# TMEM72       156      1155        13        13       111       246         7         7


# remove outlier
combat_escall2 <- combat_escall2[!(rownames(combat_escall2) %in% c("RPS27","TMEM72")),]
#combat_escall2
combat_escall2 <- round(combat_escall2)
nrow(combat_escall2) #14610
head(combat_escall2)

# temp <- combat_escall2[rowMin(combat_escall2) <0,]
# head(temp)
# ESC_slow1 ESC_slow2 ESC_slow3 ESC_slow4 ESC_fast1 ESC_fast2 ESC_fast3 ESC_fast4
# RPS27        -44       273       153       159       278       255       157       167
# TMEM72        31       820       199       198        -4       102       191       192

dds <- DESeqDataSetFromMatrix(combat_escall2,colData=coldata,design=~condition)

dds<-DESeq(dds)
head(counts(dds))
head(counts(dds_rLogTransformed))

dds_rLogTransformed <- rlog( dds ) # regularized-log transform the data , for PCA & clustering

# apply PCA
plotPCA( dds_rLogTransformed , intgroup = c( "condition","year") ) 





# add removed 0 values

esc_rb <- rbind(combat_escall2,escall0)
head(esc_rb)


# transfer to .gct file and save
esc_rb <- as.data.frame(esc_rb)
esc_rb$NAME <- rownames(esc_rb)
esc_rb$Description <- rep("na",nrow(esc_rb))
head(esc_rb)

esc_rb <- esc_rb[,c(9,10,1,2,5,6,3,4,7,8)]
colSums(esc_rb[,c(3,4,5,6,7,8,9,10)])
nrow(esc_rb)

write.table(esc_rb,file="./GSEA/esc2016+2018_rb.gct",sep = "\t",row.names = F,col.names = T,quote = F)













