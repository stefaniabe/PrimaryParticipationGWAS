args = commandArgs(trailingOnly = TRUE)
date=Sys.Date()
date=gsub("-","_",date)

#dherror path
dhpath=paste(args[1])
#mafpath
mafpath=paste(args[3])
#bimpath
bimpath=paste(args[2])
#outputpath
outpath=paste(args[4])

S1 <- as.data.frame(matrix(ncol=10,nrow=0))
names(S1) <- c("SNP","AF","CHR","POS","REF","ALT","n_DH_other1","n_DH_other0","error_1","error_0")

for(i in 1:22)
{
S1s=read.csv(paste(dhpath,"/phase_error_trios_",i,".csv.gz",sep=""),header=TRUE)
maf=read.table(paste(mafpath,"/chr",i,".frq",sep=""),header=TRUE)
maf <- merge(maf,S1s[,c(1,2)],by="SNP")
maf$A1 <- as.character(maf$A1)
maf$REF <- as.character(maf$REF)
maf$AF <- ifelse(maf$A1==maf$REF,1-maf$MAF,maf$MAF)
pos <- read.table(paste(bimpath,"/chr",i,".bim",sep=""),header=FALSE)
names(pos) <- c("CHR","SNP","b","POS","A1","A2")

PM <- merge(maf[which(colnames(maf)=="SNP"|colnames(maf)=="AF")],pos[which(colnames(pos)=="SNP"|colnames(pos)=="CHR"|colnames(pos)=="POS")],by="SNP")
S <- merge(PM,S1s,by="SNP")  

D <- rbind(S1,S)
S1 <- D
}

D <- D[order(D$CHR,D$POS),]

write.csv(D,file=gzfile(paste(outpath,"/TNT_DH_sumstat_",date,".csv.gz",sep="")),row.names=FALSE,quote=FALSE)
