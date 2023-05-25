args = commandArgs(trailingOnly = TRUE)
library(ggplot2)
mytheme=list(
theme_classic(base_size=30)+
theme(axis.text.y=element_text(face="italic")))

#Name of the pgs
pgs=paste(args[1])
#Full path to the women specific linear regression output
WL=paste(args[2])
#Full path to the women specific logistic regression output
WC=paste(args[3])
#Full path to the men specific linear regression output
ML=paste(args[4])
#Full path to the men specific logistic regression output
MC=paste(args[5])
#Path to folder to store the figure
outpath=paste(args[6])

D1 <- read.csv(WL,header=TRUE)
D1$sex=c("women")
D2 <- read.csv(WC,header=TRUE)
d <- data.frame(do.call('rbind',strsplit(as.character(D2$N),"/")))
d$N <- as.numeric(as.character(d$X1))+as.numeric(as.character(d$X2))
D2$N=d$N
D2$sex=c("women")
D3 <- read.csv(ML,header=TRUE)
D3$sex=c("men")
D4 <- read.csv(MC,header=TRUE)
d <- data.frame(do.call('rbind',strsplit(as.character(D4$N),"/")))
d$N <- as.numeric(as.character(d$X1))+as.numeric(as.character(d$X2))
D4$N=d$N
D4$sex=c("men")

DD <- rbind(D1,D2,D3,D4)
DD <- DD[which(DD$pgs==paste(pgs)),]
DD$SE = abs(DD$beta)/abs(DD$T)
DD$SE1 <- DD$beta - 1.96*DD$SE
DD$SE2 <- DD$beta + 1.96*DD$SE
DD <- na.omit(DD)
ym1=min(DD$SE1)-0.015
ym1b=ym1+0.001
ym2=max(DD$SE2)+0.003

DD$label=paste("n=",DD$N,sep="")
mytitle=expression(paste("    Estimated effect/log(OR) per 1-SD higher ",italic(pPGS)))

pdf(paste(outpath,"/",pgs,"_Fig3.pdf",sep=""),width=32,height=15)
print(ggplot(DD,aes(x=pheno,y=beta,colour=sex,label=label))+
geom_errorbar(aes(ymin=SE1,ymax=SE2),width=.2,alpha=0.8,position = position_dodge(width = 0.4))+
geom_point(alpha=1,size=3,position = position_dodge(width = 0.4))+
mytheme+
scale_color_manual(values=c("cornflowerblue","red"))+
geom_hline(yintercept=0,linetype="dashed")+
ylim(ym1,ym2)+
coord_flip()+
geom_text(aes(y=ym1b),show.legend=FALSE,size=7.5,position = position_dodge(width = 0.8))+
labs(y=mytitle,x="Phenotype",colour="sex"))
dev.off()
