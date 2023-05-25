args = commandArgs(trailingOnly = TRUE)
#N sibling pairs
N <-as.numeric(paste(args[1]))
#What is your target participation rate
a=as.numeric(paste(args[2]))
#output path
outpath <- paste(args[3])
#name of run
name <- paste(args[4])
#What is the allele frequency of the single SNP
p=as.numeric(paste(args[5]))
##How many simulations
sim=as.numeric(paste(args[6]))
r=sim-1

library(ggplot2)

#Focusing on the effect of one SNP with its standardized genotype denoted by z, the liability of sib ij is modelled as
#X_ij= w_1*z_ij+ w_A*A_i+w_B*B_ij
#where A_i and B_ij are standard normal variables.
#The variables  A_i,B_i1 and B_i2 are assumed to be independent of z_i1 and z_i2, and independent of each other.
#The term A_i  is meant to capture effects from shared environment as well as shared genetic factors other than g. 
#We assume w_1^2+w_A^2+w_B^2=1 so that var(X_ij )=1. Because cor(z_i1,z_i2 )=0.5, cor(X_i1,X_i2 )=w_A^2+(w_1^2)/2.

#Fix w_1^2=0.001 and choose 16 different values of w_A:

wAs=c(0.0005,0.025,0.050,0.075,0.100,0.125,0.150,0.175,0.193,0.225,0.250,0.275,0.300,0.325,0.350,0.375)
w_1=sqrt(0.001)

dir.create(paste(outpath,"/ibd_a_",a,"_p_",p,"_N_",N,sep=""))

########################################################
#simulating the genotypes

D0=as.data.frame(matrix(ncol=5,nrow=0.25*N))  
names(D0)=c("s1m","s1p","s2m","s2p","ibd")
D0$s1m=rbinom(0.25*N,1,p)
D0$s1p=rbinom(0.25*N,1,p)
D0$s2m=rbinom(0.25*N,1,p)
D0$s2p=rbinom(0.25*N,1,p)
D0$ibd=c("ibd0")
  
D1=as.data.frame(matrix(ncol=5,nrow=0.5*N))  
names(D1)=c("s1m","s1p","s2m","s2p","ibd")
D1$s1m=rbinom(0.5*N,1,p)
D1$s1p=rbinom(0.5*N,1,p)
D1$s2m=D1$s1m
D1$s2p=rbinom(0.5*N,1,p)
D1$ibd=c("ibd1")  

D2=as.data.frame(matrix(ncol=5,nrow=0.25*N))  
names(D2)=c("s1m","s1p","s2m","s2p","ibd")
D2$s1m=rbinom(0.25*N,1,p)
D2$s1p=rbinom(0.25*N,1,p)
D2$s2m=D2$s1m
D2$s2p=D2$s1p
D2$ibd=c("ibd2")
  
D <- rbind(D0,D1,D2)
D$pair=paste("s",1:nrow(D),sep="")

D$z1=(D$s1m+D$s1p-2*p)/sqrt(2*p*(1-p))
D$z2=(D$s2m+D$s2p-2*p)/sqrt(2*p*(1-p))
########################################################################
#simulating the shared complex component
D$A = rnorm(nrow(D))
#Simulating the non-shared complex components
D$B1 = rnorm(nrow(D))
D$B2 = rnorm(nrow(D))
           
SR <- as.data.frame(matrix(ncol=22,nrow=length(wAs)))
colnames(SR) <-  c("partfract","partfractsib","sibenrich","EAF","b0","b1","EAF_sample","EAF_notsample",
                 "EAF_sample_sib","EAF_sample_notsib","corsample",
                 "ibd0P","ibd0N","ibd0_N","ibd1sP","ibd1sN","ibd1nP","ibd1nN","ibd1_N","ibd2P","ibd2N","ibd2_N")

tt=qnorm(1-a)

for(s in 1:length(wAs)){
print(s)
wA=sqrt(wAs[s]-w_1^2/2)
wB=sqrt(1-wA^2-w_1^2)
##Defining X
D$x1=w_1*D$z1+wA*D$A+wB*D$B1
D$x2=w_1*D$z2+wA*D$A+wB*D$B2
D$part1=ifelse(D$x1>tt,1,0)
D$part2=ifelse(D$x2>tt,1,0)
  
W0=D[D$ibd=="ibd0",]
W1=D[D$ibd=="ibd1",]
W2=D[D$ibd=="ibd2",]
  
SR$partfract[s]=(mean(D$part1)+mean(D$part2))/2
SR$partfractsib[s]=(nrow(D[D$part1==1&D$part2==1,]))/nrow(D)
SR$sibenrich[s]=SR$partfractsib[s]/(SR$partfract[s]^2)
SR$b0[s]=wA
SR$b1[s]=wB
SR$EAF[s]=p
SR$EAF_sample[s]= mean(c(D[D$part1==1,]$s1m,D[D$part1==1,]$s1p,D[D$part2==1,]$s2m,D[D$part2==1,]$s2p))

SR$EAF_sample_sib[s]= mean(c(D[D$part1==1&D$part2==1,]$s1m,D[D$part1==1&D$part2==1,]$s1p,D[D$part1==1&D$part2==1,]$s2m,D[D$part1==1&D$part2==1,]$s2p))
SR$EAF_sample_notsib[s]= mean(c(D[D$part1==1&D$part2==0,]$s1m,D[D$part1==1&D$part2==0,]$s1p,D[D$part1==0&D$part2==1,]$s2m,D[D$part1==0&D$part2==1,]$s2p))
  
SR$EAF_notsample[s]= mean(c(D[D$part1==0,]$s1m,D[D$part1==0,]$s1p,D[D$part2==0,]$s2m,D[D$part2==0,]$s2p))
SR$corsample[s]=(cor(D[D$part1==1,]$s1m,D[D$part1==1,]$s1p)+cor(D[D$part2==1,]$s2m,D[D$part2==1,]$s2p))/2

SR$ibd0P[s]=(mean(c(W0[W0$part1==1&W0$part2==1,]$s1m,W0[W0$part1==1&W0$part2==1,]$s1p))+
             mean(c(W0[W0$part1==1&W0$part2==1,]$s2m,W0[W0$part1==1&W0$part2==1,]$s2p)))/2

SR$ibd0N[s]=(mean(c(W0[W0$part1==1&W0$part2==0,]$s1m,W0[W0$part1==1&W0$part2==0,]$s1p))+
             mean(c(W0[W0$part1==0&W0$part2==1,]$s2m,W0[W0$part1==0&W0$part2==1,]$s2p)))/2

SR$ibd1sP[s]=mean(c(W1[W1$part1==1&W1$part2==1,]$s1m))

SR$ibd1sN[s]=(mean(c(W1[W1$part1==1&W1$part2==0,]$s1m))+mean(c(W1[W1$part1==0&W1$part2==1,]$s1m)))/2

SR$ibd1nP[s]=(mean(c(W1[W1$part1==1&W1$part2==1,]$s1p))+mean(c(W1[W1$part1==1&W1$part2==1,]$s2p)))/2

SR$ibd1nN[s]=(mean(c(W1[W1$part1==1&W1$part2==0,]$s1p))+
             mean(c(W1[W1$part1==0&W1$part2==1,]$s2p)))/2
  
  
SR$ibd2P[s]=mean(c(W2[W2$part1==1&W2$part2==1,]$s1m,W2[W2$part1==1&W2$part2==1,]$s1p))

SR$ibd2N[s]=(mean(c(W2[W2$part1==1&W2$part2==0,]$s1m,W2[W2$part1==1&W2$part2==0,]$s1p))
             +mean(c(W2[W2$part1==0&W2$part2==1,]$s1m,W2[W2$part1==0&W2$part2==1,]$s1p)))/2

  SR$ibd2_N[s]=nrow(W2[W2$part1==1&W2$part2==1,])
  SR$ibd1_N[s]=nrow(W1[W1$part1==1&W1$part2==1,])
   SR$ibd0_N[s]=nrow(W0[W0$part1==1&W0$part2==1,])
}

SR3=SR

for(j in 1:r)
{
D0=as.data.frame(matrix(ncol=5,nrow=0.25*N))  
names(D0)=c("s1m","s1p","s2m","s2p","ibd")
D0$s1m=rbinom(0.25*N,1,p)
D0$s1p=rbinom(0.25*N,1,p)
D0$s2m=rbinom(0.25*N,1,p)
D0$s2p=rbinom(0.25*N,1,p)
D0$ibd=c("ibd0")
  
D1=as.data.frame(matrix(ncol=5,nrow=0.5*N))  
names(D1)=c("s1m","s1p","s2m","s2p","ibd")
D1$s1m=rbinom(0.5*N,1,p)
D1$s1p=rbinom(0.5*N,1,p)
D1$s2m=D1$s1m
D1$s2p=rbinom(0.5*N,1,p)
D1$ibd=c("ibd1")
  

D2=as.data.frame(matrix(ncol=5,nrow=0.25*N))  
names(D2)=c("s1m","s1p","s2m","s2p","ibd")
D2$s1m=rbinom(0.25*N,1,p)
D2$s1p=rbinom(0.25*N,1,p)
D2$s2m=D2$s1m
D2$s2p=D2$s1p
D2$ibd=c("ibd2")
  
D <- rbind(D0,D1,D2)
D$pair=paste("s",1:nrow(D),sep="")

D$z1=(D$s1m+D$s1p-2*p)/sqrt(2*p*(1-p))
D$z2=(D$s2m+D$s2p-2*p)/sqrt(2*p*(1-p))
########################################################################
#simulating the shared complex component
D$A = rnorm(nrow(D))
#Simulating the non-shared complex components
D$B1 = rnorm(nrow(D))
D$B2 = rnorm(nrow(D))
           
SR <- as.data.frame(matrix(ncol=22,nrow=length(wAs)))
colnames(SR) <-  c("partfract","partfractsib","sibenrich","EAF","b0","b1","EAF_sample","EAF_notsample",
                 "EAF_sample_sib","EAF_sample_notsib","corsample",
                 "ibd0P","ibd0N","ibd0_N","ibd1sP","ibd1sN","ibd1nP","ibd1nN","ibd1_N","ibd2P","ibd2N","ibd2_N")

tt=qnorm(1-a)

for(s in 1:length(wAs)){
print(s)
wA=sqrt(wAs[s]-w_1^2/2)
wB=sqrt(1-wA^2-w_1^2)
##Defining X
D$x1=w_1*D$z1+wA*D$A+wB*D$B1
D$x2=w_1*D$z2+wA*D$A+wB*D$B2
D$part1=ifelse(D$x1>tt,1,0)
D$part2=ifelse(D$x2>tt,1,0)
  
W0=D[D$ibd=="ibd0",]
W1=D[D$ibd=="ibd1",]
W2=D[D$ibd=="ibd2",]
  
SR$partfract[s]=(mean(D$part1)+mean(D$part2))/2
SR$partfractsib[s]=(nrow(D[D$part1==1&D$part2==1,]))/nrow(D)
SR$sibenrich[s]=SR$partfractsib[s]/(SR$partfract[s]^2)
SR$b0[s]=wA
SR$b1[s]=wB
SR$EAF[s]=p
SR$EAF_sample[s]= mean(c(D[D$part1==1,]$s1m,D[D$part1==1,]$s1p,D[D$part2==1,]$s2m,D[D$part2==1,]$s2p))

SR$EAF_sample_sib[s]= mean(c(D[D$part1==1&D$part2==1,]$s1m,D[D$part1==1&D$part2==1,]$s1p,D[D$part1==1&D$part2==1,]$s2m,D[D$part1==1&D$part2==1,]$s2p))
SR$EAF_sample_notsib[s]= mean(c(D[D$part1==1&D$part2==0,]$s1m,D[D$part1==1&D$part2==0,]$s1p,D[D$part1==0&D$part2==1,]$s2m,D[D$part1==0&D$part2==1,]$s2p))
  
SR$EAF_notsample[s]= mean(c(D[D$part1==0,]$s1m,D[D$part1==0,]$s1p,D[D$part2==0,]$s2m,D[D$part2==0,]$s2p))
SR$corsample[s]=(cor(D[D$part1==1,]$s1m,D[D$part1==1,]$s1p)+cor(D[D$part2==1,]$s2m,D[D$part2==1,]$s2p))/2

SR$ibd0P[s]=(mean(c(W0[W0$part1==1&W0$part2==1,]$s1m,W0[W0$part1==1&W0$part2==1,]$s1p))+
             mean(c(W0[W0$part1==1&W0$part2==1,]$s2m,W0[W0$part1==1&W0$part2==1,]$s2p)))/2

SR$ibd0N[s]=(mean(c(W0[W0$part1==1&W0$part2==0,]$s1m,W0[W0$part1==1&W0$part2==0,]$s1p))+
             mean(c(W0[W0$part1==0&W0$part2==1,]$s2m,W0[W0$part1==0&W0$part2==1,]$s2p)))/2

SR$ibd1sP[s]=mean(c(W1[W1$part1==1&W1$part2==1,]$s1m))

SR$ibd1sN[s]=(mean(c(W1[W1$part1==1&W1$part2==0,]$s1m))+mean(c(W1[W1$part1==0&W1$part2==1,]$s1m)))/2

SR$ibd1nP[s]=(mean(c(W1[W1$part1==1&W1$part2==1,]$s1p))+mean(c(W1[W1$part1==1&W1$part2==1,]$s2p)))/2

SR$ibd1nN[s]=(mean(c(W1[W1$part1==1&W1$part2==0,]$s1p))+
             mean(c(W1[W1$part1==0&W1$part2==1,]$s2p)))/2
  
  
SR$ibd2P[s]=mean(c(W2[W2$part1==1&W2$part2==1,]$s1m,W2[W2$part1==1&W2$part2==1,]$s1p))

SR$ibd2N[s]=(mean(c(W2[W2$part1==1&W2$part2==0,]$s1m,W2[W2$part1==1&W2$part2==0,]$s1p))
             +mean(c(W2[W2$part1==0&W2$part2==1,]$s1m,W2[W2$part1==0&W2$part2==1,]$s1p)))/2

  SR$ibd2_N[s]=nrow(W2[W2$part1==1&W2$part2==1,])
  SR$ibd1_N[s]=nrow(W1[W1$part1==1&W1$part2==1,])
   SR$ibd0_N[s]=nrow(W0[W0$part1==1&W0$part2==1,])
}
BR2=SR+SR3
SR3=BR2
}
SR=SR3/(r+1)

write.table(SR,paste(outpath,"/ibd_a_",a,"_p_",p,"_N_",N,"/ibd_p_",p,"_N_",N,"_",name,".txt",sep=""),
            row.names=FALSE,col.names=TRUE,quote=FALSE)


########################Making figure 4 from the output##########################
#Note in the paper the simulations were repeated 500 times before the figures were made
#If computational resources are limited one might want to simulate first many SR matrices parallele, then join them in the end and make the figures
D <- SR
D$lambda=D$sibenrich
D$y1=D$EAF_sample-p
D$sib=(D$EAF_sample_sib-D$EAF_sample_notsib)/D$y1
D$ibd2=(D$ibd2P-D$ibd0P)/D$y1
D$ibd1=(D$ibd1sP-D$ibd1nP)/D$y1

D1=D[which(colnames(D)=="lambda"|colnames(D)=="sib")]
names(D1)=c("xx","yy")
D1$group="siblings - singletons"

D2=D[which(colnames(D)=="lambda"|colnames(D)=="ibd2")]
names(D2)=c("xx","yy")
D2$group="IBD2-IBD0"

D3=D[which(colnames(D)=="lambda"|colnames(D)=="ibd1")]
names(D3)=c("xx","yy")
D3$group="IBD1s-IBD1n"

my <- "Relative frequency difference"
mx <- expression(lambda[s])

s1=expression((F[IBD1S]-F[IBD1NS])/(F[samp]-f[pop]))
s2=expression((F[IBD2]-F[IBD0])/(F[samp]-f[pop])    )
s3=expression((F[SIBS]-F[SING])/(F[samp]-f[pop]))

DD <- rbind(D1,D2,D3)

library(ggplot2)

jpeg(paste(outpath,"/ibd_a_",a,"_p_",p,"_N_",N,"/ibd_p_",p,"_N_",N,"_",name,"_figure_4.jpg",sep=""),width=14,height=8,units="in",res=400)
ggplot(DD[DD$xx<3.5,],aes(x=xx,y=yy,colour=group,shape=group))+
geom_point(size=3,alpha=1)+
geom_smooth(se=FALSE,linetype="dashed",size=0.2)+
xlab(mx)+
ylab(my)+
ylim(0,1.1)+
geom_vline(xintercept=2,linetype="dashed",size=0.2)+
scale_color_manual(labels=c("siblings - singletons"=s3,"IBD2-IBD0"=s2,"IBD1s-IBD1n"=s1),values=c("red","cornflowerblue","#999999"))+
scale_shape_manual(values=c(0,17,16),labels=c("siblings - singletons"=s3,"IBD2-IBD0"=s2,"IBD1s-IBD1n"=s1))+
labs(colour="",shape="")+
theme_classic(base_size=20)
dev.off()

########################Making Supplementary figure 8 from the output##########################
##Figure with all the ratios
D$sib=(D$EAF_sample_sib-0.5)/D$y1
D$sing=(D$EAF_sample_notsib-0.5)/D$y1
D$ibd2=(D$ibd2P-0.5)/D$y1
D$ibd0=(D$ibd0P-0.5)/D$y1
D$ibd1s=(D$ibd1sP-0.5)/D$y1
D$ibd1n=(D$ibd1nP-0.5)/D$y1

D1a=D[which(colnames(D)=="lambda"|colnames(D)=="sib")]
names(D1a)=c("xx","yy")
D1a$group="siblings"

D1b=D[which(colnames(D)=="lambda"|colnames(D)=="sing")]
names(D1b)=c("xx","yy")
D1b$group="singletons"

D2a=D[which(colnames(D)=="lambda"|colnames(D)=="ibd2")]
names(D2a)=c("xx","yy")
D2a$group="IBD2"

D2b=D[which(colnames(D)=="lambda"|colnames(D)=="ibd0")]
names(D2b)=c("xx","yy")
D2b$group="IBD0"

D3a=D[which(colnames(D)=="lambda"|colnames(D)=="ibd1s")]
names(D3a)=c("xx","yy")
D3a$group="IBD1s"

D3b=D[which(colnames(D)=="lambda"|colnames(D)=="ibd1n")]
names(D3b)=c("xx","yy")
D3b$group="IBD1n"

my <- "Relative frequency difference"
mx <- expression(lambda[s])

s1a=expression((F[IBD1S]-f[pop])/(F[samp]-f[pop]) )
s1b=expression((F[IBD1NS]-f[pop])/(F[samp]-f[pop]))
s2a=expression((F[IBD2]-f[pop])/(F[samp]-f[pop]) )
s2b=expression((F[IBD0]-f[pop])/(F[samp]-f[pop])  )
s3a=expression((F[SIBS]-f[pop])/(F[samp]-f[pop])  )
s3b=expression((F[SING]-f[pop])/(F[samp]-f[pop])  )

DD <- rbind(D1a,D1b,D2a,D2b,D3a,D3b)
DD$group <- factor(DD$group,levels=c("IBD2","IBD1s","siblings","singletons","IBD0","IBD1n"))

jpeg(paste(outpath,"/ibd_a_",a,"_p_",p,"_N_",N,"/ibd_p_",p,"_N_",N,"_",name,"_figure_S8.jpg",sep=""),width=14,height=8,units="in",res=400)
ggplot(DD[DD$xx<3.5,],aes(x=xx,y=yy,colour=group,shape=group))+
geom_point(size=3,alpha=1)+
geom_smooth(se=FALSE,linetype="dashed",size=0.2)+
xlab(mx)+
ylab(my)+
#ylim(0,1.1)+
geom_vline(xintercept=2,size=0.2,colour="grey")+
geom_hline(yintercept=1,size=0.2,colour="grey")+
scale_color_manual(labels=c("siblings"=s3a,"singletons"=s3b,"IBD2"=s2a,"IBD0"=s2b,"IBD1s"=s1a,"IBD1n"=s1b),values=c("red","cornflowerblue","black","red","cornflowerblue","black"))+
scale_shape_manual(values=c(0,17,8,4,16,5),labels=c("siblings"=s3a,"singletons"=s3b,"IBD2"=s2a,"IBD0"=s2b,"IBD1s"=s1a,"IBD1n"=s1b))+
labs(colour="",shape="")+
theme_classic(base_size=20)
dev.off()
