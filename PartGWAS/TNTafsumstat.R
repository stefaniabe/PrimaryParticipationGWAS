args = commandArgs(trailingOnly = TRUE)
date=Sys.Date()
date=gsub("-","_",date)

#allele-frequency path
afpath=paste(args[1])
#variance path
varpath=paste(args[2])
#mafpath
mafpath=paste(args[4])
#bimpath
bimpath=paste(args[3])
#outputpath
outpath=paste(args[5])

S1 <- as.data.frame(matrix(ncol=11,nrow=0))
names(S1) <- c("SNP","AF","CHR","POS","REF","ALT","ibd1n","ib1nmean","ibd1s","ib1smean",'ib1diffvar')

for(i in 1:22)
{
S1s=read.csv(paste(afpath,"/TNT_af1s_",i,".csv.gz",sep=""),header=TRUE)
S1n=read.csv(paste(afpath,"/TNT_af1n_",i,".csv.gz",sep=""),header=TRUE)
S1var=read.csv(paste(varpath,"/TNT_var_",i,".csv.gz",sep=""),header=TRUE)
maf=read.table(paste(mafpath,"/chr",i,".frq",sep=""),header=TRUE)
maf <- merge(maf,S1s[,c(1,2)],by="SNP")
maf$A1 <- as.character(maf$A1)
maf$REF <- as.character(maf$REF)
maf$AF <- ifelse(maf$A1==maf$REF,1-maf$MAF,maf$MAF)
pos <- read.table(paste(bimpath,"/chr",i,".bim",sep=""),header=FALSE)
names(pos) <- c("CHR","SNP","b","POS","A1","A2")

S2 <- merge(S1n,S1s,by=c("SNP","REF","ALT"))
Sns <- merge(S2,S1var[,c(1,2,3,5)],by=c("SNP","REF","ALT"))
PM <- merge(pos[which(colnames(pos)=="SNP"|colnames(pos)=="CHR"|colnames(pos)=="POS")],maf[which(colnames(maf)=="SNP"|colnames(maf)=="AF")],by="SNP")
S <- merge(PM,Sns,by="SNP")  

D <- rbind(S1,S)
S1 <- D
}

D <- D[order(D$CHR,D$POS),]

write.csv(D,file=gzfile(paste(outpath,"/TNT_AF_sumstat_",date,".csv.gz",sep="")),row.names=FALSE,quote=FALSE)
