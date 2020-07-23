# filter gene sets 

# yeast vs. fib

yeast_fib <- read.table("./yeast_fib.tsv",header = T,sep = "\t")
head(yeast_fib)

#ribosome
yeast_fib_ribosome <- yeast_fib[grep("ribosom",yeast_fib$NAME.FIB,ignore.case = T),]

#mitoch
yeast_fib_mitoch <- yeast_fib[grep("mitoch",yeast_fib$NAME.FIB,ignore.case = T),]

#dna damage
yeast_fib_dna <- yeast_fib[grep("dna",yeast_fib$NAME.FIB,ignore.case = T),]
yeast_fib_dnadamage <- yeast_fib_dna[grep("damage",yeast_fib_dna$NAME.FIB,ignore.case = T),]
# damage
yeast_fib_damage <- yeast_fib[grep("damage",yeast_fib$NAME.FIB,ignore.case = T),]


#repair
yeast_fib_repair <- yeast_fib[grep("repair",yeast_fib$NAME.FIB,ignore.case = T),]

#stress
yeast_fib_stress <- yeast_fib[grep("stress",yeast_fib$NAME.FIB,ignore.case = T),]

#CCN
yeast_fib_CCN <- yeast_fib[grep("ccn",yeast_fib$NAME.FIB,ignore.case = T),]
#E2F
yeast_fib_E2F <- yeast_fib[grep("E2F",yeast_fib$NAME.FIB,ignore.case = T),]
#p53
yeast_fib_p53 <- yeast_fib[grep("p53",yeast_fib$NAME.FIB,ignore.case = T),]
#CDK
yeast_fib_CDK<- yeast_fib[grep("CDK",yeast_fib$NAME.FIB,ignore.case = T),]


yeast_fib_sets <- rbind(yeast_fib_ribosome,yeast_fib_mitoch,yeast_fib_dnadamage,yeast_fib_damage,yeast_fib_repair,yeast_fib_stress,yeast_fib_CCN,yeast_fib_E2F,yeast_fib_p53,yeast_fib_CDK)

write.table(yeast_fib_sets,file="./yeast_fib_sets.tsv",col.names = T,row.names = F,sep = "\t",quote=F)


# fib  VS. ESCS

fib_ESCs <- read.table("./fib_ESCs.tsv",header = T,sep = "\t")
head(fib_ESCs)

#ribosome
fib_ESCs_ribosome <- fib_ESCs[grep("ribo",fib_ESCs$NAME.FIB,ignore.case = T),]

#mitoch
fib_ESCs_mitoch <- fib_ESCs[grep("mitoch",fib_ESCs$NAME.FIB,ignore.case = T),]

#dna damage
fib_ESCs_dna <- fib_ESCs[grep("dna",fib_ESCs$NAME.FIB,ignore.case = T),]
fib_ESCs_dnadamage <- fib_ESCs_dna[grep("damage",fib_ESCs_dna$NAME.FIB,ignore.case = T),]
# damage
fib_ESCs_damage <- fib_ESCs[grep("damage",fib_ESCs$NAME.FIB,ignore.case = T),]


#repair
fib_ESCs_repair <- fib_ESCs[grep("repair",fib_ESCs$NAME.FIB,ignore.case = T),]

#stress
fib_ESCs_stress <- fib_ESCs[grep("stress",fib_ESCs$NAME.FIB,ignore.case = T),]

#CCN
fib_ESCs_CCN <- fib_ESCs[grep("ccn",fib_ESCs$NAME.FIB,ignore.case = T),]
#E2F
fib_ESCs_E2F <- fib_ESCs[grep("E2F",fib_ESCs$NAME.FIB,ignore.case = T),]
#p53
fib_ESCs_p53 <- fib_ESCs[grep("p53",fib_ESCs$NAME.FIB,ignore.case = T),]
#CDK
fib_ESCs_CDK<- fib_ESCs[grep("CDK",fib_ESCs$NAME.FIB,ignore.case = T),]


fib_ESCs_sets <- rbind(fib_ESCs_ribosome,fib_ESCs_mitoch,fib_ESCs_dnadamage,fib_ESCs_damage,fib_ESCs_repair,fib_ESCs_stress,fib_ESCs_CCN,fib_ESCs_E2F,fib_ESCs_p53,fib_ESCs_CDK)

write.table(fib_ESCs_sets,file="./fib_ESCs_sets.tsv",col.names = T,row.names = F,sep = "\t",quote=F)


# all
# yeast vs. fib

yeast_fib_ESCs <- read.table("./yeast_fib_ESCs.tsv",header = T,sep = "\t")
head(yeast_fib_ESCs)

#ribosome
yeast_fib_ESCs_ribosome <- yeast_fib_ESCs[grep("ribosom",yeast_fib_ESCs$NAME.FIB,ignore.case = T),]

#mitoch
yeast_fib_ESCs_mitoch <- yeast_fib_ESCs[grep("mitoch",yeast_fib_ESCs$NAME.FIB,ignore.case = T),]

#dna damage
yeast_fib_ESCs_dna <- yeast_fib_ESCs[grep("dna",yeast_fib_ESCs$NAME.FIB,ignore.case = T),]
yeast_fib_ESCs_dnadamage <- yeast_fib_ESCs_dna[grep("damage",yeast_fib_ESCs_dna$NAME.FIB,ignore.case = T),]
# damage
yeast_fib_ESCs_damage <- yeast_fib_ESCs[grep("damage",yeast_fib_ESCs$NAME.FIB,ignore.case = T),]


#repair
yeast_fib_ESCs_repair <- yeast_fib_ESCs[grep("repair",yeast_fib_ESCs$NAME.FIB,ignore.case = T),]

#stress
yeast_fib_ESCs_stress <- yeast_fib_ESCs[grep("stress",yeast_fib_ESCs$NAME.FIB,ignore.case = T),]

#CCN
yeast_fib_ESCs_CCN <- yeast_fib_ESCs[grep("ccn",yeast_fib_ESCs$NAME.FIB,ignore.case = T),]
#E2F
yeast_fib_ESCs_E2F <- yeast_fib_ESCs[grep("E2F",yeast_fib_ESCs$NAME.FIB,ignore.case = T),]
#p53
yeast_fib_ESCs_p53 <- yeast_fib_ESCs[grep("p53",yeast_fib_ESCs$NAME.FIB,ignore.case = T),]
#CDK
yeast_fib_ESCs_CDK<- yeast_fib_ESCs[grep("CDK",yeast_fib_ESCs$NAME.FIB,ignore.case = T),]


yeast_fib_ESCs_sets <- rbind(yeast_fib_ESCs_ribosome,yeast_fib_ESCs_mitoch,yeast_fib_ESCs_dnadamage,yeast_fib_ESCs_damage,yeast_fib_ESCs_repair,yeast_fib_ESCs_stress,yeast_fib_ESCs_CCN,yeast_fib_ESCs_E2F,yeast_fib_ESCs_p53,yeast_fib_ESCs_CDK)

write.table(yeast_fib_ESCs_sets,file="./yeast_fib_ESCs_sets.tsv",col.names = T,row.names = F,sep = "\t",quote=F)
