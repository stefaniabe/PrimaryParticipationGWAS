args = commandArgs(trailingOnly = TRUE)
#path to the folder with a matrix with raw phenotype values, should have columns "IID" and then the other column names should be the phenoytpes
datapath=paste(args[1])
#output path, path to a folder
path=paste(args[2],"/",sep="")
#full path to individuals to include, tab seperated files, needs to have column IID or eid
indpath=paste(args[3])
#full path to a file with information about chip-type and sex
#Needs to have columns IID, sex, Age, YOB, PC1-PC40 and chip
pcpath=paste(args[4])
#phenotypes we want to transform, need to match to the columnheaders in the phenotype matrix
phenopath=paste(args[5])
#Phenotypes that we want to RINT and adjust for Age up to the order of 3 and YOB up to the order of 1, column name ="pheno"
phenopathage=paste(args[6])
#Phenotypes we do not want to RINT (just standardise) and use YOB up to the order of 3 and Age up to the order of 1 (this is just EA in the paper)
phenonotrank=paste(args[7])
#Phenotypes that we want to RINT and adjust for YOB up to the order of 3 and Age up to the order of 1, column name ="pheno"
phenopathyob=paste(args[8])

library(plyr)
date <- Sys.Date()
date <- gsub("-","_",date)

#Functions for RINT
#Ties are equal in rank
std <- function(x)
{
W <- !is.na(x)
z <- qnorm(rank(x[W])/(length(na.omit(x))+1))
y <- x
y[W] <- z
return(y)
}

#Ties are randomized in the rank
rstd <- function(x)
{
W <- !is.na(x)
z <- qnorm(rank(x[W],ties.method="random")/(length(na.omit(x))+1))
y <- x
y[W] <- z
return(y)
}

#list of individuals to include in the transformed phenotype list
fam <- read.table(indpath,header=TRUE)
names(fam) <- gsub("eid","IID",names(fam))

####################################
#Read in matrix with information about phenotypic measurements (one column for each phenotype), 
#Note, for blood biomarkers, in this matrix, missing values due to reportability range should be coded as -7777777, -8888888 or -9999999.
M <- read.csv(datapath,header=TRUE,nrows=10)
classes <- sapply(M,class)
M <- read.csv(datapath,header=TRUE,colClasses=classes)
names(M) <- gsub("eid","IID",names(M))

#Read in matrix with information about year of birth (YOB), age at recruitment (Age),
#sex (sex=1/2) and 40 PCs (PC1-PC40).
PC <- read.table(pcpath,header=TRUE,nrows=1)
classes <- sapply(PC,class)
PC <-  read.table(pcpath,header=TRUE,colClasses=classes)
M <- merge(M,PC,by="IID")

#Read in vector with phenotype names that I want to transform (could be a subset of the ones in table M)
#Needs to have a header with column name "pheno"
phenos=read.table(phenopath,header=TRUE)
#Want to count how often per blood biomarker we have missing values due to reportability range
dm <- as.data.frame(matrix(nrow=length(phenos$pheno),ncol=2))
names(dm) <- c("pheno","low")

for (i in 1:nrow(dm))
{dm$pheno[i]=paste(phenos$pheno[i])
k=which(colnames(M)==paste(phenos$pheno[i]))
dm$low[i]=nrow(M[which(M[,k]< -7777776),])
}

#Phenotypes defined here with low missing (fract<2500/272409)
Pheno1a <- dm[dm$low<=2500,]
#Phenotypes defined here as high missing (fract>2500/272409)
Pheno2 <- dm[dm$low>2500,]$pheno

M$Age2 <- M$Age^2
M$Age3 <- M$Age^3
M$YOB2 <- M$YOB^2
M$YOB3 <- M$YOB^3

M <- M[which(M$IID %in% fam$IID),]

#Phenotypes that we want to transform and adjust for Age up to the order of 3 and YOB up to the order of 1, column name ="pheno"
Pheno1b=read.table(phenopathage,header=TRUE)
Pheno1c <- merge(Pheno1a,Pheno1b,by="pheno")
Pheno1=Pheno1c$pheno
if(length(Pheno1)>0){
for (i in 1:length(Pheno1))
{
phs <-  Pheno1[i]
kk <- which(colnames(M)==phs)
M[,kk] <- as.numeric(as.character(M[,kk]))
DD1 <- na.omit(M[which(colnames(M)=="IID"|colnames(M)=="YOB"|colnames(M)=="sex"|colnames(M)=="Age"|colnames(M)=="Age2"|colnames(M)=="Age3"|colnames(M)==phs|grepl("PC",colnames(M)))])
#Make sure that we have at least 1000 individuals
if(nrow(DD1)>1000)
{
#Exclude values for a sex if there are fewer than 50 non-na values for that sex
dl <- table(DD1$sex)
dl2 <- dl[dl<50]
if(length(dl2>0)){D <- DD1[DD1$sex!=names(dl2),]}else{D<-DD1}

k <- which(colnames(D)==phs)
D$y <- D[,k]
D$y <- as.numeric(as.character(D$y))
###Impute "below reportable range" (-9999999) and "above reportable range" (-8888888) values for the biomarkers
miny <- min(D[D$y> -7777777,]$y,na.rm=TRUE)-0.001
maxy <- max(D[D$y> -7777777,]$y,na.rm=TRUE)+0.001
D$y <- ifelse(D$y=="-9999999",miny,ifelse(D$y=="-8888888",maxy,ifelse(D$y=="-7777777",NA,D$y)))
D$sex <- factor(D$sex)

if(length(table(D$sex))==2){
d1 <- na.omit(D[D$sex==1,])
d2 <- na.omit(D[D$sex==2,])

d1$yy <- std(d1$y)
attach(d1)
d1$resid <- resid(lm(data=d1,yy~Age+Age2+Age3+YOB+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+PC21+PC22+PC23+PC24+PC25+PC26+PC27+PC28+PC29+PC30+PC31+PC32+PC33+PC34+PC35+PC36+PC37+PC38+PC39+PC40))
d1$res=std(d1$resid)

d2$yy <- std(d2$y)
attach(d2)
d2$resid <- resid(lm(data=d2,yy~Age+Age2+Age3+YOB+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+PC21+PC22+PC23+PC24+PC25+PC26+PC27+PC28+PC29+PC30+PC31+PC32+PC33+PC34+PC35+PC36+PC37+PC38+PC39+PC40))
d2$res=std(d2$resid)

d <- rbind(d1[which(colnames(d1)=="IID"|colnames(d1)=="res")],d2[which(colnames(d2)=="IID"|colnames(d2)=="res")])

write.table(d,paste(path,phs,"_sex_age3_PC_rank_",date,".txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")

}else{ d1 <- na.omit(D)
d1$yy <- std(d1$y)
attach(d1)
d1$resid <- resid(lm(data=d1,yy~Age+Age2+Age3+YOB+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+PC21+PC22+PC23+PC24+PC25+PC26+PC27+PC28+PC29+PC30+PC31+PC32+PC33+PC34+PC35+PC36+PC37+PC38+PC39+PC40))
d1$res=std(d1$resid)

d <- d1[which(colnames(d1)=="IID"|colnames(d1)=="res")]

write.table(d,paste(path,phs,"_onesex_age3_PC_rank_",date,".txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
}
}
}
}

####Phenotypes we want to use ties=random. This group only includes blood biomarkers so we know that we want to adjust for Age up to the order of 3 and YOB up to the order of 1
if(length(Pheno2)>0){
for (i in 1:length(Pheno2))
{
phs <-  Pheno2[i]
kk <- which(colnames(M)==phs)
M[,kk] <- as.numeric(as.character(M[,kk]))
DD1 <- na.omit(M[which(colnames(M)=="IID"|colnames(M)=="YOB"|colnames(M)=="sex"|colnames(M)=="Age"|colnames(M)=="Age2"|colnames(M)=="Age3"|colnames(M)==phs|grepl("PC",colnames(M)))])
#Make sure that we have at least 1000 individuals
if(nrow(DD1)>1000)
{
dl <- table(DD1$sex)
dl2 <- dl[dl<50]
if(length(dl2>0)){D <- DD1[DD1$sex!=names(dl2),]}else{D<-DD1}

k <- which(colnames(D)==phs)
D$y <- D[,k]
D$y <- as.numeric(as.character(D$y))
###Impute "too low" and "too high" values for the biomarkers
miny <- min(D[D$y> -7777777,]$y,na.rm=TRUE)-0.001
maxy <- max(D[D$y> -7777777,]$y,na.rm=TRUE)+0.001
D$y <- ifelse(D$y=="-9999999",miny,ifelse(D$y=="-8888888",maxy,ifelse(D$y=="-7777777",NA,D$y)))

D$sex <- factor(D$sex)

if(length(table(D$sex))==2){
d1 <- na.omit(D[D$sex==1,])
d2 <- na.omit(D[D$sex==2,])

d1$yy <- std(d1$y)
attach(d1)
d1$resid <- resid(lm(data=d1,yy~Age+Age2+Age3+YOB+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+PC21+PC22+PC23+PC24+PC25+PC26+PC27+PC28+PC29+PC30+PC31+PC32+PC33+PC34+PC35+PC36+PC37+PC38+PC39+PC40))
d1$res=rstd(d1$resid)

d2$yy <- std(d2$y)
attach(d2)
d2$resid <- resid(lm(data=d2,yy~Age+Age2+Age3+YOB+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+PC21+PC22+PC23+PC24+PC25+PC26+PC27+PC28+PC29+PC30+PC31+PC32+PC33+PC34+PC35+PC36+PC37+PC38+PC39+PC40))
d2$res=rstd(d2$resid)

d <- rbind(d1[which(colnames(d1)=="IID"|colnames(d1)=="res")],d2[which(colnames(d2)=="IID"|colnames(d2)=="res")])

write.table(d,paste(path,phs,"_sex_age3_PC_randrank_",date,".txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")

}else{ d1 <- na.omit(D)
d1$yy <- std(d1$y)
attach(d1)
d1$resid <- resid(lm(data=d1,yy~Age+Age2+Age3+YOB+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+PC21+PC22+PC23+PC24+PC25+PC26+PC27+PC28+PC29+PC30+PC31+PC32+PC33+PC34+PC35+PC36+PC37+PC38+PC39+PC40))
d1$res=rstd(d1$resid)

d <- d1[which(colnames(d1)=="IID"|colnames(d1)=="res")]

write.table(d,paste(path,phs,"_onesex_age3_PC_randrank_",date,".txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
}
}
}
}

###Now make lists for phenotypes that we do not want to RINT and use YOB up to the order of 3 and Age up to the order of 1
Pheno3b=read.table(phenonotrank,header=TRUE)
Pheno3=Pheno3b$pheno

if(length(Pheno3)>0){
for (i in 1:length(Pheno3))
{
phs <-  Pheno3[i]
kk <- which(colnames(M)==phs)
M[,kk] <- as.numeric(as.character(M[,kk]))
DD1 <- na.omit(M[which(colnames(M)=="IID"|colnames(M)=="YOB"|colnames(M)=="sex"|colnames(M)=="Age"|colnames(M)=="YOB2"|colnames(M)=="YOB3"|colnames(M)==phs|grepl("PC",colnames(M)))])
#Make sure that we have at least 1000 individuals
if(nrow(DD1)>1000)
{
dl <- table(DD1$sex)
dl2 <- dl[dl<50]
if(length(dl2>0)){D <- DD1[DD1$sex!=names(dl2),]}else{D<-DD1}

k <- which(colnames(D)==phs)
D$y <- D[,k]
D$y <- as.numeric(as.character(D$y))

D$sex <- factor(D$sex)

if(length(table(D$sex))==2){
d1 <- na.omit(D[D$sex==1,])
d2 <- na.omit(D[D$sex==2,])

d1$yy <- d1$y
attach(d1)
d1$resid <- resid(lm(data=d1,yy~YOB+YOB2+YOB3+Age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+PC21+PC22+PC23+PC24+PC25+PC26+PC27+PC28+PC29+PC30+PC31+PC32+PC33+PC34+PC35+PC36+PC37+PC38+PC39+PC40))
d1$res=(d1$resid-mean(d1$resid))/sd(d1$resid)

d2$yy <- d2$y
attach(d2)
d2$resid <- resid(lm(data=d2,yy~YOB+YOB2+YOB3+Age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+PC21+PC22+PC23+PC24+PC25+PC26+PC27+PC28+PC29+PC30+PC31+PC32+PC33+PC34+PC35+PC36+PC37+PC38+PC39+PC40))
d2$res=(d2$resid-mean(d2$resid))/sd(d2$resid)

d <- rbind(d1[which(colnames(d1)=="IID"|colnames(d1)=="res")],d2[which(colnames(d2)=="IID"|colnames(d2)=="res")])

write.table(d,paste(path,phs,"_sex_yob3_PC_notrank_",date,".txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
}
else{ d1 <- na.omit(D)

d1$yy <- d1$y
attach(d1)
d1$resid <- resid(lm(data=d1,yy~YOB+YOB2+YOB3+Age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+PC21+PC22+PC23+PC24+PC25+PC26+PC27+PC28+PC29+PC30+PC31+PC32+PC33+PC34+PC35+PC36+PC37+PC38+PC39+PC40))
d1$res=(d1$resid-mean(d1$resid))/sd(d1$resid)

d <- d1[which(colnames(d1)=="IID"|colnames(d1)=="res")]

write.table(d,paste(path,phs,"_onesex_yob3_PC_notrank_",date,".txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
}
}
}
}

#Phenotypes that we want to RINT and adjust for YOB up to the order of 3 and Age up to the order of 1, column name ="pheno"
Pheno4b=read.table(phenopathyob,header=TRUE)
Pheno4c <- merge(Pheno1a,Pheno4b,by="pheno")
Pheno4=Pheno4c$pheno

if(length(Pheno4)>0){
for (i in 1:length(Pheno4))
{
phs <-  Pheno4[i]
kk <- which(colnames(M)==phs)
M[,kk] <- as.numeric(as.character(M[,kk]))
DD1 <- na.omit(M[which(colnames(M)=="IID"|colnames(M)=="YOB"|colnames(M)=="sex"|colnames(M)=="Age"|colnames(M)=="YOB2"|colnames(M)=="YOB3"|colnames(M)==phs|grepl("PC",colnames(M)))])
#Make sure that we have at least 1000 individuals
if(nrow(DD1)>1000)
{

dl <- table(DD1$sex)
dl2 <- dl[dl<50]
if(length(dl2>0)){D <- DD1[DD1$sex!=names(dl2),]}else{D<-DD1}

k <- which(colnames(D)==phs)
D$y <- D[,k]
D$y <- as.numeric(as.character(D$y))
D$sex <- factor(D$sex)

if(length(table(D$sex))==2){
d1 <- na.omit(D[D$sex=="1",])
d2 <- na.omit(D[D$sex=="2",])

d1$yy <- std(d1$y)
attach(d1)
d1$resid <- resid(lm(data=d1,yy~YOB+YOB2+YOB3+Age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+PC21+PC22+PC23+PC24+PC25+PC26+PC27+PC28+PC29+PC30+PC31+PC32+PC33+PC34+PC35+PC36+PC37+PC38+PC39+PC40))
d1$res=std(d1$resid)

d2$yy <- std(d2$y)
attach(d2)
d2$resid <- resid(lm(data=d2,yy~YOB+YOB2+YOB3+Age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+PC21+PC22+PC23+PC24+PC25+PC26+PC27+PC28+PC29+PC30+PC31+PC32+PC33+PC34+PC35+PC36+PC37+PC38+PC39+PC40))
d2$res=std(d2$resid)

d <- rbind(d1[which(colnames(d1)=="IID"|colnames(d1)=="res")],d2[which(colnames(d2)=="IID"|colnames(d2)=="res")])

write.table(d,paste(path,phs,"_sex_yob3_PC_rank_",date,".txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")

}else{ d1 <- na.omit(D)
d1$yy <- std(d1$y)
attach(d1)
d1$resid <- resid(lm(data=d1,yy~YOB+YOB2+YOB3+Age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10+PC11+PC12+PC13+PC14+PC15+PC16+PC17+PC18+PC19+PC20+PC21+PC22+PC23+PC24+PC25+PC26+PC27+PC28+PC29+PC30+PC31+PC32+PC33+PC34+PC35+PC36+PC37+PC38+PC39+PC40))
d1$res=std(d1$resid)

d <- d1[which(colnames(d1)=="IID"|colnames(d1)=="res")]

write.table(d,paste(path,phs,"_onesex_yob3_PC_rank_",date,".txt",sep=""),row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")
}
}
}
}
