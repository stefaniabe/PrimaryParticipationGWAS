args = commandArgs(trailingOnly = TRUE)
trim <- paste(args[1])
pathmaf <- paste(args[5])
pathbim <- paste(args[4])
outpath <- paste(args[6])
afpath <- paste(paste(args[2]),"/trim",trim,"_",sep="")
varpath <- paste(paste(args[3]),"/trim",trim,"_",sep="")

date=Sys.Date()
date=gsub("-","_",date)

S1 <- as.data.frame(matrix(ncol=17,nrow=0))
names(S1) <- c("SNP","CHR","POS","AF","REF","ALT","ibd0","ib0mean","ib0meanvar","ibd1n","ib1nmean","ibd1s","ib1smean","ibd1diffvar","ibd2","ib2mean","ib2meanvar")

for(i in 1:22)
{
S0=read.csv(paste(afpath,"af0_",i,".csv.gz",sep=""),header=TRUE)
S1s=read.csv(paste(afpath,"af1s_",i,".csv.gz",sep=""),header=TRUE)
S1n=read.csv(paste(afpath,"af1n_",i,".csv.gz",sep=""),header=TRUE)
S2=read.csv(paste(afpath,"af2_",i,".csv.gz",sep=""),header=TRUE)
V0=read.csv(paste(varpath,"var0_",i,".csv.gz",sep=""),header=TRUE)
V1=read.csv(paste(varpath,"var1_",i,".csv.gz",sep=""),header=TRUE)
V2=read.csv(paste(varpath,"var2_",i,".csv.gz",sep=""),header=TRUE)
maf=read.table(paste(pathmaf,"/chr",i,".frq",sep=""),header=TRUE)
maf <- merge(maf,S0[,c(1,2)],by="SNP")
maf$A1 <- as.character(maf$A1)
maf$REF <- as.character(maf$REF)
maf$AF <- ifelse(maf$A1==maf$REF,1-maf$MAF,maf$MAF)
pos <- read.table(paste(pathbim,"/chr",i,".bim",sep=""),header=FALSE)
names(pos) <- c("CHR","SNP","b","POS","A1","A2")

S0 <- merge(S0,V0[which(colnames(V0)=="SNP"|colnames(V0)=="ib0meanvar")],by="SNP")
S1s <- merge(S1s,V1[which(colnames(V1)=="SNP"|colnames(V1)=="ib1diffvar")],by="SNP")
S2 <- merge(S2,V2[which(colnames(V2)=="SNP"|colnames(V2)=="ib2meanvar")],by="SNP")  
  
Sn <- merge(S0,S1n,by=c("SNP","REF","ALT"))
Ss <- merge(S1s,S2,by=c("SNP","REF","ALT"))
Sns <- merge(Sn,Ss,by=c("SNP","REF","ALT"))
PM <- merge(pos[which(colnames(pos)=="SNP"|colnames(pos)=="CHR"|colnames(pos)=="POS")],maf[which(colnames(maf)=="SNP"|colnames(maf)=="AF")],by="SNP")
S <- merge(PM,Sns,by="SNP")

D <- rbind(S1,S)
S1 <- D
}

D <- D[order(D$CHR,D$POS),]
             
write.csv(D,file=gzfile(paste(outpath,"/trim",trim,"_AF_sumstat_",date,".csv.gz",sep="")),row.names=FALSE,quote=FALSE)
