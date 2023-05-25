args = commandArgs(trailingOnly = TRUE)
#path to the folder with all the phenotype lists, the lists should be called "phenotype.txt" and have no header and two columns, IID and phenotype
path=paste(args[1])
#output path
outpath=paste(args[2])
#full path to the pgs
#Needs to have columns IID and SCORESUM (which is the PGS)
pgspath=paste(args[3])
#Name of pgs for output files
pgs=paste(args[4])
#full path to individuals to include, tab seperated files, needs to have column IID or eid
indpath=paste(args[5])
#full path to file with a list of phenotypes to include (no header), has two two columns IID and phenotype value
phenopath=paste(args[6])
#full path to a file with information about chip-type and sex
#Needs to have columns IID, sex and chip
pcpath=paste(args[7])
#Path the folder with files with LD score regression intercept
#Need to be named L2_all.txt, L2_women.txt and L2_men.txt and include columns with names pheno and L2 (the corresponding LD score regression intercept).
L2path=paste(args[8])
#path to educational attainment phenotype list, should have two columns (no header) with IID and phenotype
edupath=paste(args[9])

#Individuals to include
d <- read.table(indpath,header=TRUE)
names(d) <- gsub("eid","IID",names(d))

#List with phenotypes you want to perform linear regression for
phenos <- read.table(phenopath, header=FALSE)
names(phenos) <- "PP"

#Make a matrix with all the transformed quantitative phenotypes
i=1
P <- read.table(paste(path,"/",phenos$PP[i],".txt",sep=""),header=FALSE)
names(P) <- c("IID",paste(phenos$PP[i]))

for(i in 2:nrow(phenos))
{
p1 <- read.table(paste(path,"/",phenos$PP[i],".txt",sep=""),header=FALSE)
names(p1) <- c("IID",paste(phenos$PP[i]))
P <- merge(P,p1,by="IID",all=TRUE)
}

#Add matrix with information about chip-type (needs to have columns called IID and chip)
PC <- read.table(pcpath,header=TRUE,nrows=1)
classes <- sapply(PC,class)
PC <-  read.table(pcpath,header=TRUE,colClasses=classes)
P <- merge(P,PC,by="IID",all.x=TRUE)

#Add the pgs
D1 <- read.table(pgspath,header=TRUE)
P <- merge(D1,P,by="IID")
P$chip <- as.factor(P$chip)

P <- P[which(P$IID %in% d$IID),]

#Output - not adjusted for pheno
PRTc1 <- as.data.frame(matrix(ncol=5,nrow=(nrow(phenos))))
names(PRTc1) <- c("pheno","T","beta","P","N")
PRTc1$pheno <- phenos$PP
PRTc1$pgs <- paste(pgs)

#Output - adjsuted for pheno
PRTc2 <- as.data.frame(matrix(ncol=5,nrow=(nrow(phenos))))
names(PRTc2) <- c("pheno","T","beta","P","N")
PRTc2$pheno <- phenos$PP
PRTc2$pgs <- paste(pgs,"_adjEDU",sep="")

EDU <- read.table(edupath,header=FALSE)
names(EDU) <- c("IID","EDU")
P <- merge(P,EDU,by="IID",all.x=TRUE)

P$PGS <- (P$SCORESUM-mean(P$SCORESUM))/sd(P$SCORESUM)

##Want to do this for all, women and men
sex <- c("all","women","men")

for(j in 1:length(sex))
{
group=sex[j]
#Read in table with LD score regression intercept
#header is pheno and L2
E <- read.table(paste(L2path,"/L2_",group,".txt",sep=""),header=TRUE)
if(group=="women"){ind=PC[PC$sex==2,]$IID}else{if(group=="men"){ind=PC[PC$sex==1,]$IID}else{ind=PC$IID}}
PP <- P[which(P$IID %in% ind),]
for(i in 1:length(phenos$PP))
{
pheno <- phenos$PP[i]
l2=E[E$pheno==paste(pheno),]$L2
DP1 <- PP[which(colnames(PP)==paste(pheno)|colnames(PP)=="PGS"|colnames(PP)=="chip")]
DP1 <- na.omit(DP1)
if(nrow(DP1)>100){
names(DP1) <- gsub(paste(pheno),"yyy",names(DP1))
sum1 <- summary(lm(DP1$yyy~DP1$PGS+DP1$chip))

#degrees of freedom
f2=nrow(DP1)-3
PRTc1$beta[i] <- sum1$coefficients[2,1]
PRTc1$T[i] <- sqrt(sum1$coefficients[2,3]^2/l2)
PRTc1$P[i] <- 2*pt(PRTc1$T[i],df=f2,lower.tail=FALSE)
PRTc1$N[i] <- nrow(DP1)

DP1 <- PP[which(colnames(PP)==paste(pheno)|colnames(PP)=="PGS"|colnames(PP)=="chip"|colnames(PP)=="EDU")]
DP1$edu <- (DP1$EDU-mean(DP1$EDU,na.rm=TRUE))/sd(DP1$EDU,na.rm=TRUE)
DP1 <- na.omit(DP1)
names(DP1) <- gsub(paste(pheno),"yyy",names(DP1))
sum1 <- summary(lm(DP1$yyy~DP1$PGS+DP1$edu+DP1$chip))

#degrees of freedom
f2=nrow(DP1)-4
PRTc2$beta[i] <- sum1$coefficients[2,1]
PRTc2$T[i] <- sqrt(sum1$coefficients[2,3]^2/l2)
PRTc2$P[i] <- 2*pt(PRTc2$T[i],df=f2,lower.tail=FALSE)
PRTc2$N[i] <- nrow(DP1)
}else{
PRTc1$beta[i] <- NA
PRTc1$T[i] <- NA
PRTc1$P[i] <- NA
PRTc1$N[i] <- NA
PRTc2$beta[i] <- NA
PRTc2$T[i] <- NA
PRTc2$P[i] <- NA
PRTc2$N[i] <- NA
}
}


PRc=rbind(PRTc1,PRTc2)
print(PRc)

date <- Sys.Date()
date <- gsub("-","_",date)

write.table(PRc,paste(outpath,"/",pgs,"_PGS_pheno_",group,"_",date,".csv",sep=""),sep=",",quote=FALSE,row.names=FALSE)
}
