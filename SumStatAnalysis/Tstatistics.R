rm(list=ls())
args = commandArgs(trailingOnly = TRUE)
date=Sys.Date()
date=gsub("-","_",date)
library(ggplot2)

#Path to the primary participatin summary statistic file based on siblings
pathsib=paste(args[1])
#Path to the primary participation summary statistic file based on PO pairs
pathPO=paste(args[2])
#Path to the error summary statistics based on the PO trios
pathDH=paste(args[3])
#Dataframe showing the rsnames of the 500,632 high quality SNPs in a column called SNP
pathSNPs=paste(args[4])
#Path to folder where to store results
outpath=paste(args[5])

#Which SNPs to use
W <- read.table(pathSNPs, header=TRUE)

#################################################################################
#Sibling T-statistics
D <- read.csv(paste(pathsib),header=TRUE)
###calculating f difference and its variance, then the t-stat, for ibd1 and then ibd20
#Frequency difference
D$df1<- D$ib1smean-D$ib1nmean

#Frequency difference
D$df20<-D$ib2mean-D$ib0mean
#t-statistic computed with the right variance
D$tdf1 <- D$df1/sqrt(D$ib1diffvar/(D$ibd1n))
D$tdf20 <- D$df20/sqrt((D$ib0meanvar/D$ibd0)+(D$ib2meanvar/D$ibd2))

#################################################################################
#TNT T-statistics
DPO<-read.csv(paste(pathPO,sep=""), header=TRUE)
	
#Frequency difference
DPO$df1<-DPO$ib1smean-DPO$ib1nmean
DPO$tdf1<-DPO$df1/sqrt(DPO$ib1diffvar/(DPO$ibd1n))
####################################################################################
#SE^2
DPO$vpo <- DPO$ib1diffvar/DPO$ibd1n
D$vsib1 <- D$ib1diffvar/D$ibd1n
D$vsib2 <- (D$ib0meanvar/D$ibd0)+(D$ib2meanvar/D$ibd2)
####################################################################################
############Adjusting for allele frequency for IBD=1 sib and PO#####################
############################Centred allele frequency################################
###Compute the overall sibling IBD1 allele frequency
#ibd1 frequency
names(DPO) <- gsub("tdf1","tdf1_po",names(DPO))
D$fo1 <- (D$ib1nmean+D$ib1smean)/2
#Compute the overall allele frequency in the parents
DPO$fo1 <- (DPO$ib1nmean+DPO$ib1smean)/2
#centered allele frequency
D$cfo1 <- D$fo1-0.5
DPO$cfo1 <- DPO$fo1-0.5

#Only use the high-quality sequence variants in the regression
modelAtdf1=summary(lm(D[D$SNP %in% W$SNP,]$tdf1 ~ 0 + D[D$SNP %in% W$SNP,]$cfo1))
modelAtdf1po=summary(lm(DPO[DPO$SNP %in% W$SNP,]$tdf1_po ~ 0 + DPO[DPO$SNP %in% W$SNP,]$cfo1))

D$t1=D$tdf1-D$cfo1*modelAtdf1$coefficients[1,1]
DPO$t1_po=DPO$tdf1_po-DPO$cfo1*modelAtdf1po$coefficients[1,1]
###########################Minor allele frequency##################################
D$maf <- ifelse(D$fo1>0.5,1-D$fo1,D$fo1)
DPO$maf <- ifelse(DPO$fo1>0.5,1-DPO$fo1,DPO$fo1)
D$maf2=D$maf^2
DPO$maf2=DPO$maf^2
D$maf3=D$maf^3
DPO$maf3=DPO$maf^3

D$chisq=D$t1^2
DPO$chisq=DPO$t1_po^2

chisqamodel=summary(lm(D[D$SNP %in% W$SNP,]$chisq~D[D$SNP %in% W$SNP,]$maf+D[D$SNP %in% W$SNP,]$maf2+D[D$SNP %in% W$SNP,]$maf3))

D$fit=chisqamodel$coefficients[1,1]+chisqamodel$coefficients[2,1]*D$maf+chisqamodel$coefficients[3,1]*D$maf2+chisqamodel$coefficients[4,1]*D$maf3
D$chisqx=D$chisq/D$fit
D$t1_b=sqrt(D$chisqx)*sign(D$t1)

chisqbmodel=summary(lm(DPO[DPO$SNP %in% W$SNP,]$chisq~DPO[DPO$SNP %in% W$SNP,]$maf+DPO[DPO$SNP %in% W$SNP,]$maf2+DPO[DPO$SNP %in% W$SNP,]$maf3))
DPO$fit=chisqbmodel$coefficients[1,1]+chisqbmodel$coefficients[2,1]*DPO$maf+chisqbmodel$coefficients[3,1]*DPO$maf2+chisqbmodel$coefficients[4,1]*DPO$maf3
DPO$chisqx=DPO$chisq/DPO$fit
DPO$t1_pob=sqrt(DPO$chisqx)*sign(DPO$t1_po)

#Combining the T-statistics using their inverse SEs as weights
SC2=merge(D[which(colnames(D)=="SNP"|colnames(D)=="CHR"|colnames(D)=="POS"|colnames(D)=="REF"|colnames(D)=="ALT"|colnames(D)=="t1_b"|colnames(D)=="tdf20"|colnames(D)=="vsib1"|colnames(D)=="vsib2")],DPO[which(colnames(DPO)=="SNP"|colnames(DPO)=="CHR"|colnames(DPO)=="POS"|colnames(DPO)=="REF"|colnames(DPO)=="ALT"|colnames(DPO)=="t1_pob"|colnames(DPO)=="vpo")],by=c("SNP","CHR","POS","REF","ALT"))
SC2=SC2[order(SC2$CHR,SC2$POS),]
SC2$Combined <- (SC2$t1_b*sqrt(1/SC2$vsib1)+(SC2$tdf20*sqrt(1/SC2$vsib2))+(SC2$t1_pob*sqrt(1/SC2$vpo)))/sqrt((1/SC2$vsib1)+(1/SC2$vsib2)+(1/SC2$vpo))

#Write out a file with the T/Z-statistic used for computing PGS and P-values
names(SC2) <- gsub("tdf20","BSPC",names(SC2))
names(SC2) <- gsub("t1_b","WSPC",names(SC2))
names(SC2) <- gsub("t1_pob","TNTC",names(SC2))

SC2 <- SC2[which(colnames(SC2)=="SNP"|colnames(SC2)=="CHR"|colnames(SC2)=="POS"|colnames(SC2)=="REF"|colnames(SC2)=="ALT"|colnames(SC2)=="BSPC"|colnames(SC2)=="WSPC"|colnames(SC2)=="TNTC"|colnames(SC2)=="Combined")]
write.table(SC2,paste(outpath,"/T_stats_",date,".txt",sep=""),row.names=FALSE,quote=FALSE)

###################Make a figure with the fit for the paper#######################
#Need coefficients for the t20 statistic
Ds <- D[D$SNP %in% W$SNP,]
Ds$chisq2=Ds$tdf20^2
Ds$afibd2=(Ds$ib2mean+Ds$ib0mean)/2
Ds$mafibd2=ifelse(Ds$afibd2>0.5,1-Ds$afibd2,Ds$afibd2)
chisqcmodel=summary(lm(Ds$chisq2~Ds$mafibd2))

DD1=as.data.frame(matrix(nrow=501,ncol=2))
DD2=as.data.frame(matrix(nrow=501,ncol=2))
DD3=as.data.frame(matrix(nrow=501,ncol=2))

names(DD1)=c("maf","chis")
DD1$maf=seq(0,500)/1000
DD1$chis=chisqamodel$coefficients[1,1]+chisqamodel$coefficients[2,1]*DD1$maf+chisqamodel$coefficients[3,1]*DD1$maf^2+chisqamodel$coefficients[4,1]*DD1$maf^3

names(DD2)=c("maf","chis")
DD2$maf=seq(0,500)/1000
DD2$chis=chisqbmodel$coefficients[1,1]+chisqbmodel$coefficients[2,1]*DD2$maf+chisqbmodel$coefficients[3,1]*DD2$maf^2+chisqbmodel$coefficients[4,1]*DD2$maf^3

names(DD2)=c("maf","chis")
names(DD3)=c("maf","chis")
DD3$maf=seq(0,500)/1000
DD3$chis=chisqcmodel$coefficients[1,1]+chisqcmodel$coefficients[2,1]*DD3$maf

DD1$Statistic=c("WSPC")
DD2$Statistic=c("TNTC")
DD3$Statistic=c("BSPC")
 
DD <- rbind(DD1,DD2,DD3)
DD$Statistic=factor(DD$Statistic,levels=c("WSPC","TNTC","BSPC"))

s3=expression(Chi^2~paste("(fitted average)"))

jpeg(paste(outpath,"/ExtFig7.jpg",sep=""),width=8,height=8,res=600,units='in')
ggplot(data=DD,aes(x=maf,y=chis, colour=Statistic,linetype=Statistic))+
geom_line(alpha=1)+
scale_colour_manual(values=c("blue","red","black"))+
scale_linetype_manual(values=c(1,1,2))+
theme_classic(base_size=15)+
labs(x="MAF",y=s3)+
labs(color="",linetype="")+
ylim(1,1.14)
dev.off()

ph <- read.csv(pathDH,header=TRUE)
ph <- ph[(ph$SNP %in% W$SNP),]
##Fraction of times when the shared allele is truly 1 (based on the genotype of the other parent) but is called 0
ph$e0 <- ph$error_0/ph$n_DH_other0
##Fraction of times when the shared allele is truly 0 (based on the genotype of the other parent) but is called 1
ph$e1 <- ph$error_1/ph$n_DH_other1
D1=ph[which(colnames(ph)=="AF"|colnames(ph)=="e1"|colnames(ph)=="n_DH_other1")]
D0=ph[which(colnames(ph)=="AF"|colnames(ph)=="e0"|colnames(ph)=="n_DH_other0")]
names(D1) = c("AF","n","error")
names(D0) = c("AF","n","error")
D1$group=c("True shared allele is 0")
D0$group=c("True shared allele is 1")
D1$AF2=1-D1$AF
D2 <- D1[,c(5,2,3,4)]
names(D2) <- gsub("AF2","AF",names(D2))
DD <- rbind(D2,D0)
DD <- DD[DD$n>0,]
DD$AF2=DD$AF^2
DD$AF3=DD$AF^3
model1=summary(lm(DD$error~DD$AF+DD$AF2+DD$AF3,weights=DD$n))
DD$fit2=model1$coefficients[1,1]+model1$coefficients[2,1]*DD$AF+model1$coefficients[3,1]*DD$AF2+model1$coefficients[4,1]*DD$AF3

jpeg(paste(outpath,"/extFig5.jpg",sep=""),width=6,height=6,unit='in',res=600)
ggplot(data=DD,aes(x=AF,y=fit2))+
geom_smooth(se=F)+
theme_classic(base_size=15)+
labs(x="Frequency of allele coded as 1",y="Error rate when shared allele is 1")+
labs(color=" ")
dev.off()

#Need the unadjusted WSPC statistic
DD <- merge(D,ph,by=c("SNP","CHR","POS","REF","ALT"))
DD$cfo3=DD$cfo1^3
DD$cfo5=DD$cfo1^5

#Relationship between the unadjusted WSPC statistic and centred allele frequency
model2=summary(lm(DD$tdf1~DD$cfo1+DD$cfo3-1))

DD$SE=sqrt(DD$vsib1)
DD$AF22=sqrt(DD$fo1*(1-DD$fo1)*2)
#Model the relationship between expected SE under HWE and actual SE
model3=summary(lm(DD$SE~DD$AF22))

#Plot the relationship on a grid
E <- as.data.frame(matrix(nrow=1000,ncol=3))
names(E) <- c("AF","empirical","phasing")
E$AF=seq(1:1000)/1000
E$CF=E$AF-0.5
E$CF3=E$CF^3
E$AF22=sqrt(2*E$AF*(1-E$AF))

#Use coefficients from fit in model 3 to get expected SE given the allele frequency
E$SE1=model3$coefficients[1,1]+model3$coefficients[2,1]*E$AF22

#Frequency for the allele coded as 1
E$AF2=E$AF^2
E$AF3=E$AF^3
#Frequency for the allele coded as 0
E$aAF=(1-E$AF)
E$aAF2=(1-E$AF)^2
E$aAF3=(1-E$AF)^3

E$ef=model1$coefficients[1,1]+model1$coefficients[2,1]*E$AF+model1$coefficients[3,1]*E$AF2+model1$coefficients[4,1]*E$AF3
E$ef1=model1$coefficients[1,1]+model1$coefficients[2,1]*E$aAF+model1$coefficients[3,1]*E$aAF2+model1$coefficients[4,1]*E$aAF3

E$phasing=(2*E$AF*(1-E$AF))*(E$AF*E$ef1-((1-E$AF)*E$ef))/(E$SE1)
E$empirical=model2$coefficients[1,1]*E$CF+model2$coefficients[2,1]*E$CF3

d1 <- E[which(colnames(E)=="AF"|colnames(E)=="empirical")]
d2 <- E[which(colnames(E)=="AF"|colnames(E)=="phasing")]
names(d1)=c("AF","eff")
names(d2)=c("AF","eff")
d1$group=c("empirical")
d2$group=c("phasing")
dd <- rbind(d1,d2)
dd$group=as.factor(dd$group)

#Restrict from 1 to 99%
d <- dd[dd$AF>=0.01&dd$AF<=0.99,]

jpeg(paste(outpath,"/extFig6.jpg",sep=""),width=6.5,height=6,unit='in',res=600)
ggplot(data=d,aes(x=AF,y=eff,linetype=group))+
geom_line()+
scale_colour_manual(values=c("orange","blue"))+
theme_classic(base_size=15)+
labs(x="Frequency of allele coded as 1",y="Bias in the WSPC t-statistic")+
labs(linetype=" ")
dev.off()
