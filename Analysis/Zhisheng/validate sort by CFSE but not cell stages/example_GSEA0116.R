
# check ribosome
 nes <- read.table("NES_all.tsv",header = T,sep = "\t") 

head(nes)

temp <- nes[grep("ribosom",nes$NAME.FIB,ignore.case = T),]

