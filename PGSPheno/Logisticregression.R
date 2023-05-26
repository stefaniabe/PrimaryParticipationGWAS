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
#full path to file with a list of phenotypes to include (no header)
phenopath=paste(args[6])
#full path to a file with information about covariates
pcpath=paste(args[7])
#Path to the folder with files with LD score regression intercepts
#Need to be named L2_all.txt, L2_women.txt and L2_men.txt and include columns with names pheno and L2 (the corresponding LD score regression intercept).
L2path=paste(args[8])
#path to educational attainment phenotype list, should have two columns (no header) with IID and phenotype
edupath=paste(args[9])

#Individuals to include
d <- read.table(indpath,header=TRUE)
names(d) <- gsub("eid","IID",names(d))

#List with phenotypes you want to perform logistic regression for
phenos <- read.table(phenopath, header=FALSE)
names(phenos) <- "PP"

#Make a matrix with all the binary phenotypes
i=1
P <- read.table(paste(path,"/",phenos$PP[i],".txt",sep=""),header=FALSE)
names(P) <- c("IID",paste(phenos$PP[i]))

for(i in 2:nrow(phenos))
{
p1 <- read.table(paste(path,"/",phenos$PP[i],".txt",sep=""),header=FALSE)
names(p1) <- c("IID",paste(phenos$PP[i]))
P <- merge(P,p1,by="IID",all=TRUE)
}

#Add matrix with information about chip-type (chip), sex, Age, YOB and PCs (Identifier column is called "IID").
PC <- read.table(pcpath,header=TRUE,nrows=1)
classes <- sapply(PC,class)
PC <-  read.table(pcpath,header=TRUE,colClasses=classes)
P <- merge(P,PC,by="IID",all.x=TRUE)
P$Age2=P$Age^2
P$Age3=P$Age^3

#Add the pgs
D1 <- read.table(pgspath,header=TRUE)
P <- merge(D1,P,by="IID")
P$chip <- as.factor(P$chip)
P$sex <- as.factor(P$sex)

P <- P[which(P$IID %in% d$IID),]

#Output - not adjusted for EDU
PRTc1 <- as.data.frame(matrix(ncol=5,nrow=(nrow(phenos))))
names(PRTc1) <- c("pheno","T","beta","P","N")
PRTc1$pheno <- phenos$PP
PRTc1$pgs <- paste(pgs)

#Output - adjusted for EDU
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
#header has to have columns named pheno and L2
E <- read.table(paste(L2path,"/L2_",group,".txt",sep=""),header=TRUE)
if(group=="women"){ind=PC[PC$sex==2,]$IID}else{if(group=="men"){ind=PC[PC$sex==1,]$IID}else{ind=PC$IID}}
PP <- P[which(P$IID %in% ind),]
for(i in 1:length(phenos$PP))
{
pheno <- phenos$PP[i]
l2=E[E$pheno==paste(pheno),]$L2
DP1 <- na.omit(PP[which(colnames(PP)==paste(pheno)|colnames(PP)=="PGS"|grepl("PC",colnames(PP))|colnames(PP)=="sex"|colnames(PP)=="YOB"|colnames(PP)=="Age"|colnames(PP)=="Age2"|colnames(PP)=="Age3"|colnames(PP)=="chip")])

if(nrow(DP1)>100){

names(DP1) <- gsub(paste(pheno),"yyy",names(DP1))
DP1$sex <- as.factor(DP1$sex)

if(group=="all"){
sum1 <- summary(glm(data=DP1,yyy~PGS+sex+YOB+Age+Age2+Age3+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+PC21+PC22+PC23+PC24+PC25+PC26+PC27+PC28+PC29+PC30+PC31+PC32+PC33+PC34+PC35+PC36+PC37+PC38+PC39+PC40+chip,family="binomial"))
}else{
sum1 <- summary(glm(data=DP1,yyy~PGS+YOB+Age+Age2+Age3+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+PC21+PC22+PC23+PC24+PC25+PC26+PC27+PC28+PC29+PC30+PC31+PC32+PC33+PC34+PC35+PC36+PC37+PC38+PC39+PC40+chip,family="binomial"))
}

PRTc1$beta[i] <- sum1$coefficients[2,1]
PRTc1$T[i] <- sqrt(sum1$coefficients[2,3]^2/l2)
PRTc1$P[i] <- 2*pnorm(PRTc1$T[i],lower.tail=FALSE)
PRTc1$N[i] <- paste(table(DP1$yyy)[2],"/",table(DP1$yyy)[1],sep="")

##
DP1 <- na.omit(PP[which(colnames(PP)==paste(pheno)|colnames(PP)=="PGS" |colnames(PP)=="EDU"|grepl("PC",colnames(PP))|colnames(PP)=="sex"|colnames(PP)=="YOB"|colnames(PP)=="Age"|colnames(PP)=="Age2"|colnames(PP)=="Age3"|colnames(PP)=="chip")])
names(DP1) <- gsub(paste(pheno),"yyy",names(DP1))
DP1 <- na.omit(DP1)
if(group=="all"){
sum1 <- summary(glm(data=DP1,yyy~PGS+EDU+sex+YOB+Age+Age2+Age3+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+PC21+PC22+PC23+PC24+PC25+PC26+PC27+PC28+PC29+PC30+PC31+PC32+PC33+PC34+PC35+PC36+PC37+PC38+PC39+PC40+chip,family="binomial"))
}else{
sum1 <- summary(glm(data=DP1,yyy~PGS+EDU+YOB+Age+Age2+Age3+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+PC21+PC22+PC23+PC24+PC25+PC26+PC27+PC28+PC29+PC30+PC31+PC32+PC33+PC34+PC35+PC36+PC37+PC38+PC39+PC40+chip,family="binomial"))
}

PRTc2$beta[i] <- sum1$coefficients[2,1]
PRTc2$T[i] <- sqrt(sum1$coefficients[2,3]^2/l2)
PRTc2$P[i] <- 2*pnorm(PRTc2$T[i],lower.tail=FALSE)
PRTc2$N[i] <- paste(table(DP1$yyy)[2],"/",table(DP1$yyy)[1],sep="")
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

write.table(PRc,paste(outpath,"/",pgs,"_PGS_pheno_cc_",group,"_",date,".csv",sep=""),sep=",",quote=FALSE,row.names=FALSE)
}
