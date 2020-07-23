#from lucas
---
  title: "TMRE sorting expression analysis"
---
  Code to load the TMRE sorted FIB & ESCs and do some differential expression analysis

(1) load packages
(2) set the working directory where the sample table and Kalliso .h5 files are
```{r setup}
# load libraries
library(tximportData)
library(tidyverse)
library(GenomicFeatures)
library(DESeq2)
library(gplots)
library(RColorBrewer)
library(tximport)
KallistoFilesFolder <- "~/Downloads/Proliferation_Expression__SequencingData/DataTMRE/"
ProjectGitHubFolder <- '~/Develop/GrowthPredictionTMRE/'
knitr::opts_knit$set(root.dir = ProjectGitHubFolder )
```

build Transcript2Gene table for mouse
```{r}
# load mouse table from Stephen Turner ( https://zenodo.org/record/1324497#.XCnCDM8zbOQ )
txdb <- loadDb('~/Data/Mouse/gencode.vM18.annotation.sqlite')
tx2gene <- select(txdb, keys(txdb, keytype = "TXNAME") , "GENEID", "TXNAME")
```


This function is necessary to read kallisto .h5 files that are not named abundance.h5
```{r}
read_kallisto_h5 <- function(fpath, ...) {
  if (!requireNamespace("rhdf5", quietly=TRUE)) {
    stop("reading kallisto results from hdf5 files requires Bioconductor package `rhdf5`")
  }
  counts <- rhdf5::h5read(fpath, "est_counts")
  ids <- rhdf5::h5read(fpath, "aux/ids")
  efflens <- rhdf5::h5read(fpath, "aux/eff_lengths")
  # as suggested by https://support.bioconductor.org/p/96958/#101090
  ids <- as.character(ids)
  
  stopifnot(length(counts) == length(ids)) 
  stopifnot(length(efflens) == length(ids))
  result <- data.frame(target_id = ids,
                       eff_length = efflens,
                       est_counts = counts,
                       stringsAsFactors = FALSE)
  normfac <- with(result, (1e6)/sum(est_counts/eff_length))
  result$tpm <- with(result, normfac*(est_counts/eff_length))
  return(result)
}
```

Load data sequencing data and the meta data (experiment info) table
```{r}
samples <- read.table( paste(ProjectGitHubFolder , "Data/SamplesTMRE2018.txt" , sep=''), header = TRUE)
kallisto_files <- paste( KallistoFilesFolder , as.character(samples$FileName) , sep='/') 
data <- tximport( kallisto_files , type = "kallisto", txOut = FALSE , tx2gene=tx2gene , ignoreAfterBar=TRUE , importer = read_kallisto_h5 )
```


```{r}
# use DESeq2 to find genes differentially expressed by TMRE level across both cell types
# http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#differential-expression-analysis
# 
rownames(samples) <- colnames(data$counts)
dds <- DESeqDataSetFromTximport(data, samples, design = ~TMRE+CellType+TMRE:CellType)
dds <- DESeq(dds)
dds_rLogTransformed <- rlog( dds ) # regularized-log transform the data , for PCA & clustering
resTMRE <- results(dds, contrast=c("TMRE","high","low"))
resCellType <- results(dds, contrast=c("CellType","ES","FIB"))
```

Do cell types and treatment cluster together? 
  ```{r}
# using a distance matrix plot
sampleDists <- dist( t( assay(dds_rLogTransformed) ) )
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( dds_rLogTransformed$TMRE,dds_rLogTransformed$CellType,dds_rLogTransformed$Rep, sep="-" )
colnames(sampleDistMatrix) <- NULL
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatmap.2( sampleDistMatrix, trace="none", col=colours)
# using PCA -- results don't look very promising
plotPCA( dds_rLogTransformed , intgroup = c( "CellType", "TMRE") ) 
```
# only check a single cell type for PCA
```{r} 
samples_CT <- samples[samples$CellType=='FIB',]
kallisto_files_CT <- paste( KallistoFilesFolder , as.character(samples_CT$FileName) , sep='/') 
data_CT <- tximport( kallisto_files_CT , type = "kallisto", txOut = FALSE , tx2gene=tx2gene , ignoreAfterBar=TRUE , importer = read_kallisto_h5 )
dds_CT <- DESeqDataSetFromTximport(data_CT, samples_CT, design = ~TMRE)
dds_CT <- DESeq(dds_CT)
dds_rLogTransformed_CT <- rlog( dds_CT ) # regularized-log transform the data , for PCA & clustering
rownames(samples_CT) <- colnames(data_CT$counts)
# here they look better
plotPCA( dds_rLogTransformed_CT , intgroup = c( "CellType", "TMRE") ) 
# and for ESCs
samples_CT <- samples[samples$CellType=='ES' & samples$Rep>1,]
kallisto_files_CT <- paste( KallistoFilesFolder , as.character(samples_CT$FileName) , sep='/') 
data_CT <- tximport( kallisto_files_CT , type = "kallisto", txOut = FALSE , tx2gene=tx2gene , ignoreAfterBar=TRUE , importer = read_kallisto_h5 )
dds_CT <- DESeqDataSetFromTximport(data_CT, samples_CT, design = ~TMRE)
dds_CT <- DESeq(dds_CT)
dds_rLogTransformed_CT <- rlog( dds_CT ) # regularized-log transform the data , for PCA & clustering
rownames(samples_CT) <- colnames(data_CT$counts)
# here they look better
plotPCA( dds_rLogTransformed_CT , intgroup = c( "CellType", "TMRE","Rep") ) 
```