## Spectrum paper
## Function plot prev; per prev
analysis<-function(dat,pheno){
  
  modComp  <- summary( glm(CC~NRA,family="binomial",data=dat) )$coefficients
  dat$grp  <- round(dat$NRA) #cut(dat$NRA,breaks = seq(min(dat$NRA),max(dat$NRA),by=2),include.lowest=TRUE)
  prev_p   <- aggregate(dat$CC,by=list(dat$grp),FUN=mean,na.rm=T)
  prev_n   <- aggregate(dat$CC,by=list(dat$grp),FUN=sum,na.rm=T)
  prev_N   <- aggregate(dat$CC,by=list(dat$grp),FUN=function(x) length(!is.na(x)))
  allele_c <- aggregate(dat$NRA,by=list(dat$grp),FUN=mean)
  d        <- cbind.data.frame(N=prev_N$x,n=prev_n$x,p=prev_p$x,numAllele=allele_c$x)
  d$p_se   <- sqrt(d$p*(1-d$p)/d$N)
  d        <- d[which(d$N>=200),]
  d$px     <- d$p / d$numAllele
  d$px_se  <- d$p_se / d$numAllele
  
  modelLogit  <- glm(cbind(n,N)~numAllele,data=d,family = binomial(link="logit"))
  modelProbit <- glm(cbind(n,N)~numAllele,data=d,family = binomial(link="probit"))
  modelLog    <- glm(cbind(n,N)~numAllele,data=d,family = binomial(link="log"))
  modelLinear <- lm(p~numAllele,data=d,weights = N)
  
  predLogit   <- predict(modelLogit,type = "response")
  predProbit  <- predict(modelProbit,type = "response")
  predLog     <- predict(modelLog,type = "response")
  predLin     <- predict(modelLinear,type="response")
  
  Rlogit      <- round(unlist(cor.test(d$p,predLogit)[c("estimate","conf.int")]),2)
  Rprobit     <- round(unlist(cor.test(d$p,predProbit)[c("estimate","conf.int")]),2)
  Rlog        <- round(unlist(cor.test(d$p,predLog)[c("estimate","conf.int")]),2)
  Rlin        <- round(unlist(cor.test(d$p,predLin)[c("estimate","conf.int")]),2)
  
  Cols <- c("seagreen","coral1","goldenrod","blue")
  l    <- which(d[,"p_se"]>0)
  
  #  png(paste0("plots/",ds,".prevalence.fit.png"),width=1500,height=1250,res=250)
  #pdf("test.pdf")
  # plot(d[l,"numAllele"],d[l,"p"],pch=19,axes=FALSE,ylim=c(0,max(d$p+d$p_se)),main=ds,xlab="Mean number of risk alleles",ylab="Prevalence");axis(1);axis(2)
  plot(d[l,"numAllele"],d[l,"p"],pch=19,axes=FALSE,ylim=c(0,max(d$p+d$p_se)),xlim=c(min(d$numAllele-sd(d$numAllele)),max(d$numAllele+sd(d$numAllele))),main=pheno,
       xlab="Mean number of risk alleles",ylab="Prevalence");axis(1);axis(2)
  for(i in 1:nrow(d)){
    if(d[i,"p_se"]>0){
      arrows(d[i,"numAllele"],d[i,"p"]-d[i,"p_se"],d[i,"numAllele"],d[i,"p"]+d[i,"p_se"],col="grey",angle=90,len=0.05)
      arrows(d[i,"numAllele"],d[i,"p"]+d[i,"p_se"],d[i,"numAllele"],d[i,"p"]-d[i,"p_se"],col="grey",angle=90,len=0.05)
    }
  }
  #  matlines(d$numAllele,cbind(predLogit,predProbit,predLog),col=Cols,lwd=2)
  #  abline(a=modelLinear$coefficients[1],b=modelLinear$coefficients[2],col=4,lty=3)
  matlines(d$numAllele,cbind(predLogit,predLog),col=Cols,lwd=2)
  abline(a=modelLog$coefficients[1],b=modelLog$coefficients[2],col=4,lty=3)
  quant<-quantile(dat$NRA,probs=c(0.1,0.9))
  polygon(x=c(quant[1],quant[1],quant[2],quant[2]),y=c(0,max(d$p),max(d$p),0),col=rgb(0.2,0.2,0.2,0.2),border = NA)
  legend("topleft",
         legend=c(paste0("Logit: R=",Rlogit[1]," (95%CI: [",Rlogit[2],"-",Rlogit[3],"])"),
                  #paste0("Probit: R=",Rprobit[1]," (95%CI: [",Rprobit[2],"-",Rprobit[3],"])"),
                  #paste0("Linear: R=",Rlin[1]," (95%CI: [",Rlin[2],"-",Rlin[3],"])"),
                  paste0("Log: R=",Rlog[1]," (95%CI: [",Rlog[2],"-",Rlog[3],"])")),
         #       col=c(Cols,"blue"),cex=0.9,box.lty=0,lty=c(1,1,1,3),lwd=2,bg="transparent")
         col=Cols,cex=0.9,box.lty=0,lwd=2,bg="transparent")
  #dev.off()
  
  ## Second  plot
  #png(paste0("plots/",ds,".prevalence_per_allele.png"),width=1500,height=1250,res=250)
  #pdf("test2.pdf")
  #plot(d[l,"numAllele"],d[l,"px"],pch=19,axes=FALSE,ylim=c(0,max(d$px+d$px_se)),main=ds,xlab="Mean number of risk alleles",ylab="Prevalence/Allele");axis(1);axis(2)
  plot(d[l,"numAllele"],d[l,"px"],pch=19,axes=FALSE,ylim=c(0,max(d$px+d$px_se)),xlim=c(min(d$numAllele-sd(d$numAllele)),max(d$numAllele+sd(d$numAllele))),main=pheno,
       xlab="Mean number of risk alleles",ylab="Prevalence/Allele");axis(1);axis(2)
  for(i in 1:nrow(d)){
    if(d[i,"px_se"]>0){
      arrows(d[i,"numAllele"],d[i,"px"]-d[i,"px_se"],d[i,"numAllele"],d[i,"px"]+d[i,"px_se"],col="grey",angle=90,len=0.05)
      arrows(d[i,"numAllele"],d[i,"px"]+d[i,"px_se"],d[i,"numAllele"],d[i,"px"]-d[i,"px_se"],col="grey",angle=90,len=0.05)
    }
  }
  #md <- lm(d[l,"px"]~d[l,"numAllele"],weights = d[l,"N"])
  #mx <- coef(md)
  #abline(a=mx[1],b=mx[2],col="grey",lty=2)
  
  #modelLinear <- lm(px~numAllele,data=d,weights = N)
  #predLin     <- predict(modelLinear,type="response")
  #Rlin        <- round(unlist(cor.test(d$p,predLin)[c("estimate","conf.int")]),2)
  modelLogit  <- glm(px~numAllele,data=d,family = binomial(link="logit"))
  modelLog    <- glm(px~numAllele,data=d,family = binomial(link="log"))
  
  predLogit   <- predict(modelLogit,type = "response")
  predLog     <- predict(modelLog,type = "response")
  
  Rlogit      <- round(unlist(cor.test(d$p,predLogit)[c("estimate","conf.int")]),2)
  Rlog        <- round(unlist(cor.test(d$p,predLog)[c("estimate","conf.int")]),2)
  matlines(d$numAllele,cbind(predLogit,predLog),col=Cols,lwd=2)
  quant<-quantile(dat$NRA,probs=c(0.1,0.9))
  polygon(x=c(quant[1],quant[1],quant[2],quant[2]),y=c(0,max(d$px),max(d$px),0),col=rgb(0.2,0.2,0.2,0.2),border = NA)
  #legend("topleft",legend=paste0("Linear: R=",Rlin[1]," (95%CI: [",Rlin[2],"-",Rlin[3],"])"),box.lty=0,lty=2,col="grey")
  legend("topleft",
         legend=c(paste0("Logit: R=",Rlogit[1]," (95%CI: [",Rlogit[2],"-",Rlogit[3],"])"),
                  #paste0("Probit: R=",Rprobit[1]," (95%CI: [",Rprobit[2],"-",Rprobit[3],"])"),
                  #paste0("Linear: R=",Rlin[1]," (95%CI: [",Rlin[2],"-",Rlin[3],"])"),
                  paste0("Log: R=",Rlog[1]," (95%CI: [",Rlog[2],"-",Rlog[3],"])")),
         col=Cols,cex=0.9,box.lty=0,lwd=2,bg="transparent")
}

## Function plot prev; per prev -wGRS
wanalysis<-function(dat,pheno){
  
  modComp  <- summary( glm(CC~wgrs,family="binomial",data=dat) )$coefficients
  #dat$grp  <- round(dat$NRA) #
  dat$grp<-cut(dat$wgrs,breaks = seq(min(dat$wgrs),max(dat$wgrs),by=0.05),include.lowest=TRUE)
  prev_p   <- aggregate(dat$CC,by=list(dat$grp),FUN=mean,na.rm=T)
  prev_n   <- aggregate(dat$CC,by=list(dat$grp),FUN=sum,na.rm=T)
  prev_N   <- aggregate(dat$CC,by=list(dat$grp),FUN=function(x) length(!is.na(x)))
  allele_c <- aggregate(dat$NRA,by=list(dat$grp),FUN=mean)
  wgrs<-aggregate(dat$wgrs,by=list(dat$grp),FUN=mean)
  d        <- cbind.data.frame(N=prev_N$x,n=prev_n$x,p=prev_p$x,numAllele=allele_c$x,wgrs=wgrs$x)
  d$p_se   <- sqrt(d$p*(1-d$p)/d$N)
  d        <- d[which(d$N>=200),]
  d$px     <- d$p / d$numAllele
  d$px_se  <- d$p_se / d$numAllele
  
  modelLogit  <- glm(cbind(n,N)~wgrs,data=d,family = binomial(link="logit"))
  modelProbit <- glm(cbind(n,N)~wgrs,data=d,family = binomial(link="probit"))
  modelLog    <- glm(cbind(n,N)~wgrs,data=d,family = binomial(link="log"))
  modelLinear <- lm(p~wgrs,data=d,weights = N)
  
  predLogit   <- predict(modelLogit,type = "response")
  predProbit  <- predict(modelProbit,type = "response")
  predLog     <- predict(modelLog,type = "response")
  predLin     <- predict(modelLinear,type="response")
  
  Rlogit      <- round(unlist(cor.test(d$p,predLogit)[c("estimate","conf.int")]),2)
  Rprobit     <- round(unlist(cor.test(d$p,predProbit)[c("estimate","conf.int")]),2)
  Rlog        <- round(unlist(cor.test(d$p,predLog)[c("estimate","conf.int")]),2)
  Rlin        <- round(unlist(cor.test(d$p,predLin)[c("estimate","conf.int")]),2)
  
  Cols <- c("seagreen","coral1","goldenrod")
  l    <- which(d[,"p_se"]>0)
  
  #  png(paste0("plots/",ds,".prevalence.fit.png"),width=1500,height=1250,res=250)
  #pdf("test.pdf")
  # plot(d[l,"numAllele"],d[l,"p"],pch=19,axes=FALSE,ylim=c(0,max(d$p+d$p_se)),main=ds,xlab="Mean number of risk alleles",ylab="Prevalence");axis(1);axis(2)
  plot(d[l,"wgrs"],d[l,"p"],pch=19,axes=FALSE,ylim=c(0,max(d$p+d$p_se)),xlim=c(min(d$wgrs-sd(d$wgrs)),max(d$wgrs+sd(d$wgrs))),main=pheno,
       xlab="Mean number of weighted GRS",ylab="Prevalence");axis(1,at=(0.5*(ceiling(min(d$wgrs-sd(d$wgrs)/2)/0.5):ceiling(max(d$wgrs+sd(d$wgrs)/2)/0.5))));axis(2)
  for(i in 1:nrow(d)){
    if(d[i,"p_se"]>0){
      arrows(d[i,"wgrs"],d[i,"p"]-d[i,"p_se"],d[i,"wgrs"],d[i,"p"]+d[i,"p_se"],col="grey",angle=90,len=0.05)
      arrows(d[i,"wgrs"],d[i,"p"]+d[i,"p_se"],d[i,"wgrs"],d[i,"p"]-d[i,"p_se"],col="grey",angle=90,len=0.05)
    }
  }
  #  matlines(d$numAllele,cbind(predLogit,predProbit,predLog),col=Cols,lwd=2)
  #  abline(a=modelLinear$coefficients[1],b=modelLinear$coefficients[2],col=4,lty=3)
  matlines(d$wgrs,cbind(predLogit,predLog),col=Cols,lwd=2)
  abline(a=modelLog$coefficients[1],b=modelLog$coefficients[2],col=4,lty=3)
  quant<-quantile(dat$wgrs,probs=c(0.1,0.9))
  polygon(x=c(quant[1],quant[1],quant[2],quant[2]),y=c(0,max(d$p),max(d$p),0),col=rgb(0.2,0.2,0.2,0.2),border = NA)
  legend("topleft",
         legend=c(paste0("Logit: R=",Rlogit[1]," (95%CI: [",Rlogit[2],"-",Rlogit[3],"])"),
                  #                paste0("Probit: R=",Rprobit[1]," (95%CI: [",Rprobit[2],"-",Rprobit[3],"])"),
                  paste0("Log: R=",Rlog[1]," (95%CI: [",Rlog[2],"-",Rlog[3],"])")),
         #                paste0("Linear: R=",Rlin[1]," (95%CI: [",Rlin[2],"-",Rlin[3],"])")),
         #       col=c(Cols,"blue"),cex=0.9,box.lty=0,lty=c(1,1,1,3),lwd=2,bg="transparent")
         col=Cols,cex=0.9,box.lty=0,lwd=2,bg="transparent")
  #dev.off()
  
  ## Second  plot
  #png(paste0("plots/",ds,".prevalence_per_allele.png"),width=1500,height=1250,res=250)
  #pdf("test2.pdf")
  #plot(d[l,"numAllele"],d[l,"px"],pch=19,axes=FALSE,ylim=c(0,max(d$px+d$px_se)),main=ds,xlab="Mean number of risk alleles",ylab="Prevalence/Allele");axis(1);axis(2)
  plot(d[l,"wgrs"],d[l,"px"],pch=19,axes=FALSE,ylim=c(0,max(d$px+d$px_se)),xlim=c(min(d$wgrs-sd(d$wgrs)),max(d$wgrs+sd(d$wgrs))),main=pheno,
       xlab="Mean number of weighted GRS",ylab="Prevalence/Allele");axis(1,at=(0.5*((ceiling(min(d$wgrs-sd(d$wgrs)/2)/0.5):ceiling(max(d$wgrs+sd(d$wgrs)/2)/0.5)))));axis(2)
  for(i in 1:nrow(d)){
    if(d[i,"px_se"]>0){
      arrows(d[i,"wgrs"],d[i,"px"]-d[i,"px_se"],d[i,"wgrs"],d[i,"px"]+d[i,"px_se"],col="grey",angle=90,len=0.05)
      arrows(d[i,"wgrs"],d[i,"px"]+d[i,"px_se"],d[i,"wgrs"],d[i,"px"]-d[i,"px_se"],col="grey",angle=90,len=0.05)
    }
  }
  #md <- lm(d[l,"px"]~d[l,"numAllele"],weights = d[l,"N"])
  #mx <- coef(md)
  #abline(a=mx[1],b=mx[2],col="grey",lty=2)
  
  #modelLinear <- lm(px~numAllele,data=d,weights = N)
  #predLin     <- predict(modelLinear,type="response")
  #Rlin        <- round(unlist(cor.test(d$p,predLin)[c("estimate","conf.int")]),2)
  modelLogit  <- glm(px~wgrs,data=d,family = binomial(link="logit"))
  modelLog    <- glm(px~wgrs,data=d,family = binomial(link="log"))
  
  predLogit   <- predict(modelLogit,type = "response")
  predLog     <- predict(modelLog,type = "response")
  
  Rlogit      <- round(unlist(cor.test(d$p,predLogit)[c("estimate","conf.int")]),2)
  Rlog        <- round(unlist(cor.test(d$p,predLog)[c("estimate","conf.int")]),2)
  matlines(d$wgrs,cbind(predLogit,predLog),col=Cols,lwd=2)
  quant<-quantile(dat$wgrs,probs=c(0.1,0.9))
  polygon(x=c(quant[1],quant[1],quant[2],quant[2]),y=c(0,max(d$px),max(d$px),0),col=rgb(0.2,0.2,0.2,0.2),border = NA)
  #legend("topleft",legend=paste0("Linear: R=",Rlin[1]," (95%CI: [",Rlin[2],"-",Rlin[3],"])"),box.lty=0,lty=2,col="grey")
  legend("topleft",
         legend=c(paste0("Logit: R=",Rlogit[1]," (95%CI: [",Rlogit[2],"-",Rlogit[3],"])"),
                  paste0("Log: R=",Rlog[1]," (95%CI: [",Rlog[2],"-",Rlog[3],"])")),
         col=Cols,cex=0.9,box.lty=0,lwd=2,bg="transparent")
}

setwd("/home/pang/Desktop/Liability_zeng_HS_FB_20190923/UKBB/202_OR_1.3/flip/198")
library("tidyverse")
# read uGRS
PRS<-read.table("UKBB_198.profile",header=TRUE,as.is=TRUE,sep=" ")
dat<-PRS %>% select(1,6,3)
dat$CC<-dat$PHENO-1
dat<-dat[,-3]
colnames(dat)<-c("sample","NRA","CC")
analysis(dat,"UKBB_CAD")
dat2<-arrange(dat,NRA)
# match wGRS to uGRS
t<-fread("UKBB_198.all.score")
t2<-t[match(dat2$sample,t$FID),]
dat2$wgrs<-t2$`5e-08`
#wanalysis(dat2,"UKBB_CAD")
# Figure size 463*564
#random
PRS<-read.table("../random/UKBB_91A.profile",header=TRUE,as.is=TRUE,sep=" ")
t<-fread("../random/UKBB_91A.all.score")
PRS<-read.table("../random/UKBB_91B.profile",header=TRUE,as.is=TRUE,sep=" ")
t<-fread("../random/UKBB_91B.all.score")
#remove Lipids SNPs
PRS<-read.table("Lipids/UKBB_182_ex_LP.profile",header=TRUE,as.is=TRUE,sep=" ")
t<-fread("Lipids/UKBB_182_ex_LP.all.score")
PRS<-read.table("BloodPressure/UKBB_182_ex_BP.profile",header=TRUE,as.is=TRUE,sep=" ")
t<-fread("BloodPressure/UKBB_182_ex_BP.all.score")
#Breast Cancer
PRS<-read.table("../../../BreastCancer/UKBB_BC_female2.profile",header=TRUE,as.is=TRUE,sep=" ")
t<-fread("../../../BreastCancer/UKBB_BC_female2.all.score")
#Prostate Cancer
PRS<-read.table("../../../ProstateCancer/UKBB_PC_male.profile",header=TRUE,as.is=TRUE,sep=" ")
t<-fread("../../../ProstateCancer/UKBB_PC_male.all.score")
# merge figure into one plot
#cairo_ps(filename='test.eps', width=11.06, height=12.39,pointsize = 12,family="Univers")
cairo_ps(filename='BC_wGRS_per_01.tif', width=11.06, height=12.39,pointsize = 12,family="Univers")
nf <- layout(matrix(c(1,3,2,4),ncol=2), width=c(463,463,463,463),height=c(564,564,564,564),TRUE)
analysis(dat,"UKBB BC")
wanalysis(dat2,"UKBB BC")
dev.off()
#normal distribution test
ks.test(bcdat2$NRA, "pnorm", mean=mean(bcdat2$NRA), sd=sd(bcdat2$NRA))


# uGRS
quant<-quantile(dat2$NRA,probs = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
dat2$uorder<-10
dat2[dat2$NRA<=quant[1],5]<-1
dat2[quant[1]<dat2$NRA & dat2$NRA<=quant[2],5]<-2
dat2[quant[2]<dat2$NRA & dat2$NRA<=quant[3],5]<-3
dat2[quant[3]<dat2$NRA & dat2$NRA<=quant[4],5]<-4
dat2[quant[4]<dat2$NRA & dat2$NRA<=quant[5],5]<-5
dat2[quant[5]<dat2$NRA & dat2$NRA<=quant[6],5]<-6
dat2[quant[6]<dat2$NRA & dat2$NRA<=quant[7],5]<-7
dat2[quant[7]<dat2$NRA & dat2$NRA<=quant[8],5]<-8
dat2[quant[8]<dat2$NRA & dat2$NRA<=quant[9],5]<-9
dat2l<-split(dat2,dat2$uorder)
ures<-data.frame("nra"=0,"sd"=0,"nrac"=0,"sdc"=0,"nraco"=0,"sdco"=0,"p"=0,"ppa"=0,"N"=0,"nc"=0,"nco"=0,"wgrs"=0,"sdw"=0,"wgrsc"=0,"sdwc"=0,"wgrsco"=0,"sdwco"=0)
for (i in 1:10){
  t<-dat2l[[i]]
  nra<-mean(t$NRA) # number of risk alleles
  sd<-sd(t$NRA) # sd of nRA
  c<-t[t$CC==1,]
  nrac<-mean(c$NRA) # number of risk alleles case
  sdc<-sd(c$NRA) # sd of nRAc
  co<-t[t$CC==0,]
  nraco<-mean(co$NRA) # number of risk alleles control
  sdco<-sd(co$NRA) # sd of nRAco
  p<-nrow(c)/nrow(t) # prevalence
  ppa<-p/nra # prevalence per allele
  N<-nrow(t) # Sample N
  nc<-nrow(c) # Case
  nco<-nrow(co) # Control
  wgrs<-mean(t$wgrs) # weighted GRS
  sdw<-sd(t$wgrs) # sd of wGRS
  wgrsc<-mean(c$wgrs) # weighted GRS case
  sdwc<-sd(c$wgrs) # sd of wGRS case
  wgrsco<-mean(co$wgrs) # weighted GRS control
  sdwco<-sd(co$wgrs) # sd of wGRS control
  ures[i,]<-c(nra,sd,nrac,sdc,nraco,sdco,p,ppa,N,nc,nco,wgrs,sdw,wgrsc,sdwc,wgrsco,sdwco)
}
# Figure size 463*564
# wGRS
quant<-quantile(dat2$wgrs,probs = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
dat2$worder<-10
dat2[dat2$wgrs<=quant[1],6]<-1
dat2[quant[1]<dat2$wgrs & dat2$wgrs<=quant[2],6]<-2
dat2[quant[2]<dat2$wgrs & dat2$wgrs<=quant[3],6]<-3
dat2[quant[3]<dat2$wgrs & dat2$wgrs<=quant[4],6]<-4
dat2[quant[4]<dat2$wgrs & dat2$wgrs<=quant[5],6]<-5
dat2[quant[5]<dat2$wgrs & dat2$wgrs<=quant[6],6]<-6
dat2[quant[6]<dat2$wgrs & dat2$wgrs<=quant[7],6]<-7
dat2[quant[7]<dat2$wgrs & dat2$wgrs<=quant[8],6]<-8
dat2[quant[8]<dat2$wgrs & dat2$wgrs<=quant[9],6]<-9
dat2l<-split(dat2,dat2$worder)
wres<-data.frame("nra"=0,"sd"=0,"nrac"=0,"sdc"=0,"nraco"=0,"sdco"=0,"p"=0,"ppa"=0,"N"=0,"nc"=0,"nco"=0,"wgrs"=0,"sdw"=0,"wgrsc"=0,"sdwc"=0,"wgrsco"=0,"sdwco"=0)
for (i in 1:10){
  t<-dat2l[[i]]
  nra<-mean(t$NRA) # number of risk alleles
  sd<-sd(t$NRA) # sd of nRA
  c<-t[t$CC==1,]
  nrac<-mean(c$NRA) # number of risk alleles case
  sdc<-sd(c$NRA) # sd of nRAc
  co<-t[t$CC==0,]
  nraco<-mean(co$NRA) # number of risk alleles control
  sdco<-sd(co$NRA) # sd of nRAco
  p<-nrow(c)/nrow(t) # prevalence
  ppa<-p/nra # prevalence per allele
  N<-nrow(t) # Sample N
  nc<-nrow(c) # Case
  nco<-nrow(co) # Control
  wgrs<-mean(t$wgrs) # weighted GRS
  sdw<-sd(t$wgrs) # sd of wGRS
  wgrsc<-mean(c$wgrs) # weighted GRS case
  sdwc<-sd(c$wgrs) # sd of wGRS case
  wgrsco<-mean(co$wgrs) # weighted GRS control
  sdwco<-sd(co$wgrs) # sd of wGRS control
  wres[i,]<-c(nra,sd,nrac,sdc,nraco,sdco,p,ppa,N,nc,nco,wgrs,sdw,wgrsc,sdwc,wgrsco,sdwco)
}
# Log oddsratio
# uGRS
lor<-dat2 %>%mutate(GRE=cut_number(NRA,10))
#regression like you did, calculate lor for all quintiles
fit <- glm(CC ~ 0+GRE,data=lor,family="binomial")
# results like you have
results = coefficients(summary(fit))
# rename second column, SE for plotting
colnames(results)[2] = "SE"
results<-as.data.frame(results)
results$Estimate.ad<-results$Estimate+abs(results$Estimate[1])
results$OR<-exp(results$Estimate.ad)
results$CI.L<-exp(results$Estimate.ad-1.96*results$SE)       
results$CI.U<-exp(results$Estimate.ad+1.96*results$SE)
# wGRS
lor<-dat2 %>%mutate(GRE=cut_number(wgrs,10))
#regression like you did, calculate lor for all quintiles
fit <- glm(CC ~ 0+GRE,data=lor,family="binomial")
# results like you have
results = coefficients(summary(fit))
# rename second column, SE for plotting
colnames(results)[2] = "SE"
results<-as.data.frame(results)
results$Estimate.ad<-results$Estimate+abs(results$Estimate[1])
results$OR<-exp(results$Estimate.ad)
results$CI.L<-exp(results$Estimate.ad-1.96*results$SE)       
results$CI.U<-exp(results$Estimate.ad+1.96*results$SE)


# Figure 1 histogram
caddat2$pheno<-"CAD"
t2ddat2$pheno<-"T2D"
dist<-rbind(caddat2,t2ddat2)
ggplot(dist)+geom_histogram(aes(x=NRA,group=pheno,fill=pheno),binwidth = 1,color="white",size=0.2,position="identity")+xlab("Number of risk alleles")+ylab("Number of samples")+theme_classic()+labs(fill="Phenotype")+ annotate(geom="label", x = mean(caddat2$NRA)+15, y = 20000, label = round(mean(caddat2$NRA),digits = 2),fill="#42AD60")+ annotate(geom="label", x = mean(t2ddat2$NRA)-15, y = 5000, label = round(mean(t2ddat2$NRA),digits = 2),fill="#6DAFD8")+  scale_fill_manual(values=c("#42AD60","#6DAFD8"))+ geom_vline(xintercept = mean(caddat2$NRA), linetype="longdash",color = "black", alpha=0.8)+geom_vline(xintercept = mean(t2ddat2$NRA), linetype="longdash",color = "black",alpha=0.8)+scale_y_continuous(expand = c(0, 0))
#+ geom_rect(aes(xmin=quant1[1], xmax=quant1[2], ymin=0, ymax=21000))
#558x306
ggplot(dist)+geom_histogram(aes(x=wgrs,group=pheno,fill=pheno),binwidth = 0.1,color="white",size=0.2,position="identity")+xlab("Number of weighted GRS")+ylab("Number of samples")+theme_classic()+labs(fill="Phenotype")+ annotate(geom="label", x = mean(caddat2$wgrs)+1, y = 500, label = round(mean(caddat2$wgrs),digits = 2),fill="#42AD60")+ annotate(geom="label", x = mean(t2ddat2$wgrs)-1, y = 15000, label = round(mean(t2ddat2$wgrs),digits = 2),fill="#6DAFD8")+  scale_fill_manual(values=c("#42AD60","#6DAFD8"))+ geom_vline(xintercept = mean(caddat2$wgrs), linetype="longdash",color = "black", alpha=0.8)+geom_vline(xintercept = mean(t2ddat2$wgrs), linetype="longdash",color = "black",alpha=0.8)+scale_y_continuous(expand = c(0, 0))

#Figure 1 histogram case controls
#cairo_ps(filename='Figure1_case_control_BC.tif', width=11.06, height=12.39,pointsize = 12,family="Univers")
#nf <- layout(matrix(c(1,3,2,4),ncol=2), width=c(558,558,558,558),height=c(306,306,306,306),TRUE)
cairo_ps(filename='F1B_BC.tif', width=5.58, height=3.06,pointsize = 12,family="Univers")
cairo_pdf(filename='F1B_BC.pdf', width=5.58, height=3.06,pointsize = 12,family="Univers")
dev.off()
# same distribution
ggplot(data=caddat2,aes(x=NRA,group=CC,fill=CC))+geom_histogram(aes(y=..density..*406996),binwidth = 1,color="white",size=0.2,alpha=0.7,position="identity")+xlab("Number of risk alleles")+ylab("Number of samples")+theme_classic()+labs(fill="Phenotype")+ annotate(geom="label", x = mean(caddat2[caddat2$CC==0,]$NRA)-10, y = 18000, label = round(mean(caddat2[caddat2$CC==0,]$NRA),digits = 1),fill="#42AD60")+ annotate(geom="label", x = mean(caddat2[caddat2$CC==1,]$NRA)+20, y = 4000, label = round(mean(caddat2[caddat2$CC==1,]$NRA),digits = 1),fill="#6DAFD8")+scale_fill_manual(values=c("#42AD60","#6DAFD8"),labels = c("control","case"))+ geom_vline(xintercept = mean(caddat2[caddat2$CC==1,]$NRA), linetype="longdash",color = "black", alpha=0.8)+geom_vline(xintercept = mean(caddat2[caddat2$CC==0,]$NRA), linetype="longdash",color = "black",alpha=0.8)+ ggtitle("UKBB CAD")+scale_y_continuous(name="Number of samples (control)",sec.axis=sec_axis(~./20.03919,name="Number of samples (case)"),expand = c(0, 0))
ggplot(data=bcdat2,aes(x=NRA,group=CC,fill=CC))+geom_histogram(aes(y=..density..*205560),binwidth = 1,color="white",size=0.2,alpha=0.7,position="identity")+xlab("Number of risk alleles")+theme_classic()+labs(fill="Phenotype")+ annotate(geom="label", x = mean(bcdat2[bcdat2$CC==0,]$NRA)-10, y = 10000, label = round(mean(bcdat2[bcdat2$CC==0,]$NRA),digits = 1),fill="#42AD60")+ annotate(geom="label", x = mean(bcdat2[bcdat2$CC==1,]$NRA)+20, y = 2000, label = round(mean(bcdat2[bcdat2$CC==1,]$NRA),digits = 1),fill="#6DAFD8")+scale_fill_manual(values=c("#42AD60","#6DAFD8"),labels = c("control","case"))+ geom_vline(xintercept = mean(bcdat2[bcdat2$CC==1,]$NRA), linetype="longdash",color = "black", alpha=0.8)+geom_vline(xintercept = mean(bcdat2[bcdat2$CC==0,]$NRA), linetype="longdash",color = "black",alpha=0.8)+ ggtitle("UKBB BC")+scale_y_continuous(name="Number of samples (control)",sec.axis=sec_axis(~./15.54799,name="Number of samples (case)"),expand = c(0, 0))
ggplot(data=pcdat2,aes(x=NRA,group=CC,fill=CC))+geom_histogram(aes(y=..density..*197792),binwidth = 1,color="white",size=0.2,alpha=0.7,position="identity")+xlab("Number of risk alleles")+ylab("Number of samples")+theme_classic()+labs(fill="Phenotype")+ annotate(geom="label", x = mean(pcdat2[pcdat2$CC==0,]$NRA)-10, y = 10000, label = round(mean(pcdat2[pcdat2$CC==0,]$NRA),digits = 1),fill="#42AD60")+ annotate(geom="label", x = mean(pcdat2[pcdat2$CC==1,]$NRA)+20, y = 2000, label = round(mean(pcdat2[pcdat2$CC==1,]$NRA),digits = 1),fill="#6DAFD8")+scale_fill_manual(values=c("#42AD60","#6DAFD8"),labels = c("control","case"))+ geom_vline(xintercept = mean(pcdat2[pcdat2$CC==1,]$NRA), linetype="longdash",color = "black", alpha=0.8)+geom_vline(xintercept = mean(pcdat2[pcdat2$CC==0,]$NRA), linetype="longdash",color = "black",alpha=0.8)+ ggtitle("UKBB PC")+scale_y_continuous(name="Number of samples (control)",sec.axis=sec_axis(~./25.25434,name="Number of samples (case)"),expand = c(0, 0))

# uGRS
ggplot(data=caddat2)+geom_histogram(aes(x=NRA,group=CC,fill=CC),binwidth = 1,color="white",size=0.2,position="identity")+xlab("Number of risk alleles")+ylab("Number of samples")+theme_classic()+labs(fill="Phenotype")+ annotate(geom="label", x = mean(caddat2[caddat2$CC==0,]$NRA)-10, y = 18000, label = round(mean(caddat2[caddat2$CC==0,]$NRA),digits = 2),fill="#42AD60")+ annotate(geom="label", x = mean(caddat2[caddat2$CC==1,]$NRA)+10, y = 3000, label = round(mean(caddat2[caddat2$CC==1,]$NRA),digits = 2),fill="#6DAFD8")+scale_fill_manual(values=c("#42AD60","#6DAFD8"),labels = c("control","case"))+ geom_vline(xintercept = mean(caddat2[caddat2$CC==1,]$NRA), linetype="longdash",color = "black", alpha=0.8)+geom_vline(xintercept = mean(caddat2[caddat2$CC==0,]$NRA), linetype="longdash",color = "black",alpha=0.8)+scale_y_continuous(expand = c(0, 0))+ ggtitle("UKBB_CAD")
ggplot(data=bcdat2)+geom_histogram(aes(x=NRA,group=CC,fill=CC),binwidth = 1,color="white",size=0.2,position="identity")+xlab("Number of risk alleles")+ylab("Number of samples")+theme_classic()+labs(fill="Phenotype")+ annotate(geom="label", x = mean(bcdat2[bcdat2$CC==0,]$NRA)-10, y = 10000, label = round(mean(bcdat2[bcdat2$CC==0,]$NRA),digits = 2),fill="#42AD60")+ annotate(geom="label", x = mean(bcdat2[bcdat2$CC==1,]$NRA)+10, y = 2000, label = round(mean(bcdat2[bcdat2$CC==1,]$NRA),digits = 2),fill="#6DAFD8")+scale_fill_manual(values=c("#42AD60","#6DAFD8"),labels = c("control","case"))+ geom_vline(xintercept = mean(bcdat2[bcdat2$CC==1,]$NRA), linetype="longdash",color = "black", alpha=0.8)+geom_vline(xintercept = mean(bcdat2[bcdat2$CC==0,]$NRA), linetype="longdash",color = "black",alpha=0.8)+scale_y_continuous(expand = c(0, 0))+ ggtitle("UKBB_BC")
ggplot(data=pcdat2)+geom_histogram(aes(x=NRA,group=CC,fill=CC),binwidth = 1,color="white",size=0.2,position="identity")+xlab("Number of risk alleles")+ylab("Number of samples")+theme_classic()+labs(fill="Phenotype")+ annotate(geom="label", x = mean(pcdat2[pcdat2$CC==0,]$NRA)-10, y = 10000, label = round(mean(pcdat2[pcdat2$CC==0,]$NRA),digits = 2),fill="#42AD60")+ annotate(geom="label", x = mean(pcdat2[pcdat2$CC==1,]$NRA)+10, y = 1000, label = round(mean(pcdat2[pcdat2$CC==1,]$NRA),digits = 2),fill="#6DAFD8")+scale_fill_manual(values=c("#42AD60","#6DAFD8"),labels = c("control","case"))+ geom_vline(xintercept = mean(pcdat2[pcdat2$CC==1,]$NRA), linetype="longdash",color = "black", alpha=0.8)+geom_vline(xintercept = mean(pcdat2[pcdat2$CC==0,]$NRA), linetype="longdash",color = "black",alpha=0.8)+scale_y_continuous(expand = c(0, 0))+ ggtitle("UKBB_PC")
#wGRS
ggplot(data=caddat2)+geom_histogram(aes(x=wgrs,group=CC,fill=CC),binwidth = 0.1,color="white",size=0.2,position="identity")+xlab("Number of weighted GRS")+ylab("Number of samples")+theme_classic()+labs(fill="Phenotype")+ annotate(geom="label", x = mean(caddat2[caddat2$CC==0,]$wgrs)-1, y = 18000, label = round(mean(caddat2[caddat2$CC==0,]$wgrs),digits = 2),fill="#42AD60")+ annotate(geom="label", x = mean(caddat2[caddat2$CC==1,]$wgrs)+1, y = 3000, label = round(mean(caddat2[caddat2$CC==1,]$wgrs),digits = 2),fill="#6DAFD8")+scale_fill_manual(values=c("#42AD60","#6DAFD8"),labels = c("control","case"))+ geom_vline(xintercept = mean(caddat2[caddat2$CC==1,]$wgrs), linetype="longdash",color = "black", alpha=0.8)+geom_vline(xintercept = mean(caddat2[caddat2$CC==0,]$wgrs), linetype="longdash",color = "black",alpha=0.8)+scale_y_continuous(expand = c(0, 0))+ ggtitle("UKBB_CAD")
ggplot(data=bcdat2)+geom_histogram(aes(x=wgrs,group=CC,fill=CC),binwidth = 0.1,color="white",size=0.2,position="identity")+xlab("Number of weighted GRS")+ylab("Number of samples")+theme_classic()+labs(fill="Phenotype")+ annotate(geom="label", x = mean(bcdat2[bcdat2$CC==0,]$wgrs)-1, y = 10000, label = round(mean(bcdat2[bcdat2$CC==0,]$wgrs),digits = 2),fill="#42AD60")+ annotate(geom="label", x = mean(bcdat2[bcdat2$CC==1,]$wgrs)+1, y = 2000, label = round(mean(bcdat2[bcdat2$CC==1,]$wgrs),digits = 2),fill="#6DAFD8")+scale_fill_manual(values=c("#42AD60","#6DAFD8"),labels = c("control","case"))+ geom_vline(xintercept = mean(bcdat2[bcdat2$CC==1,]$wgrs), linetype="longdash",color = "black", alpha=0.8)+geom_vline(xintercept = mean(bcdat2[bcdat2$CC==0,]$wgrs), linetype="longdash",color = "black",alpha=0.8)+scale_y_continuous(expand = c(0, 0))+ ggtitle("UKBB_BC")
ggplot(data=pcdat2)+geom_histogram(aes(x=wgrs,group=CC,fill=CC),binwidth = 0.1,color="white",size=0.2,position="identity")+xlab("Number of weighted GRS")+ylab("Number of samples")+theme_classic()+labs(fill="Phenotype")+ annotate(geom="label", x = mean(pcdat2[pcdat2$CC==0,]$wgrs)-1, y = 10000, label = round(mean(pcdat2[pcdat2$CC==0,]$wgrs),digits = 2),fill="#42AD60")+ annotate(geom="label", x = mean(pcdat2[pcdat2$CC==1,]$wgrs)+1, y = 1000, label = round(mean(pcdat2[pcdat2$CC==1,]$wgrs),digits = 2),fill="#6DAFD8")+scale_fill_manual(values=c("#42AD60","#6DAFD8"),labels = c("control","case"))+ geom_vline(xintercept = mean(pcdat2[pcdat2$CC==1,]$wgrs), linetype="longdash",color = "black", alpha=0.8)+geom_vline(xintercept = mean(pcdat2[pcdat2$CC==0,]$wgrs), linetype="longdash",color = "black",alpha=0.8)+scale_y_continuous(expand = c(0, 0))+ ggtitle("UKBB_PC")
#library("ggpubr")
#ggarrange(cadp, bcp, pcp, labels = c("A", "B", "C"),ncol = 2, nrow = 2)

ggplot(dist)+geom_histogram(aes(x=NRA,group=pheno,fill=pheno),binwidth = 2,color="white",size=0.2,position="identity")+xlab("Number of risk alleles")+ylab("Number of samples")+theme_classic()+labs(fill="Phenotype")+ annotate(geom="label", x = mean(caddat2$NRA)+15, y = 18000, label = round(mean(caddat2$NRA),digits = 2),fill="#00BA38")+ annotate(geom="label", x = mean(pcdat2$NRA)-15, y = 10000, label = round(mean(pcdat2$NRA),digits = 2),fill="#619cff")+ geom_vline(xintercept = mean(caddat2$NRA), linetype="longdash",color = "black", alpha=0.8)+geom_vline(xintercept = mean(pcdat2$NRA), linetype="longdash",color = "black",alpha=0.8)+annotate(geom="label", x = mean(bcdat2$NRA)+15, y = 10000, label = round(mean(bcdat2$NRA),digits = 2),fill="#F8766D")+ geom_vline(xintercept = mean(bcdat2$NRA), linetype="longdash",color = "black", alpha=0.8)+scale_y_continuous(expand = c(0, 0))
#+ geom_rect(aes(xmin=quant1[1], xmax=quant1[2], ymin=0, ymax=21000))
#558x306
ggplot(dist)+geom_histogram(aes(x=wgrs,group=pheno,fill=pheno),binwidth = 0.15,color="white",size=0.2,position="identity")+xlab("Number of weighted GRS")+ylab("Number of samples")+theme_classic()+labs(fill="Phenotype")+ annotate(geom="label", x = mean(caddat2$wgrs)-1, y = 18000, label = round(mean(caddat2$wgrs),digits = 2),fill="#00BA38")+ annotate(geom="label", x = mean(bcdat2$wgrs)-1, y = 10000, label = round(mean(bcdat2$wgrs),digits = 2),fill="#F8766D")+annotate(geom="label", x = mean(pcdat2$wgrs)+1, y = 10000, label = round(mean(pcdat2$wgrs),digits = 2),fill="#619cff")+geom_vline(xintercept = mean(pcdat2$wgrs), linetype="longdash",color = "black", alpha=0.8)+geom_vline(xintercept = mean(caddat2$wgrs), linetype="longdash",color = "black", alpha=0.8)+geom_vline(xintercept = mean(bcdat2$wgrs), linetype="longdash",color = "black",alpha=0.8)+scale_y_continuous(expand = c(0, 0))

#density
ggplot(dist,aes(x=NRA,group=pheno,fill=pheno))+geom_histogram(aes(y=..density..),binwidth = 2,color="white",size=0.2,position="identity")+geom_density(aes(y=..density..),alpha=0.2)+xlab("Number of risk alleles")+ylab("Number of samples")+theme_classic()+labs(fill="Phenotype")+ annotate(geom="label", x = mean(caddat2$NRA)+15, y = 0.04, label = round(mean(caddat2$NRA),digits = 2),fill="#00BA38")+ annotate(geom="label", x = mean(pcdat2$NRA)-15, y = 0.04, label = round(mean(pcdat2$NRA),digits = 2),fill="#619cff")+ geom_vline(xintercept = mean(caddat2$NRA), linetype="longdash",color = "black", alpha=0.8)+geom_vline(xintercept = mean(pcdat2$NRA), linetype="longdash",color = "black",alpha=0.8)+annotate(geom="label", x = mean(bcdat2$NRA)+15, y = 0.04, label = round(mean(bcdat2$NRA),digits = 2),fill="#F8766D")+ geom_vline(xintercept = mean(bcdat2$NRA), linetype="longdash",color = "black", alpha=0.8)+scale_y_continuous(expand = c(0, 0))
#+ geom_rect(aes(xmin=quant1[1], xmax=quant1[2], ymin=0, ymax=21000))
#558x306
ggplot(dist,aes(x=wgrs,group=pheno,fill=pheno))+geom_histogram(aes(y=..density..),binwidth = 0.15,color="white",size=0.2,position="identity")+geom_density(aes(y=..density..),alpha=0.2)+xlab("Number of weighted GRS")+ylab("Number of samples")+theme_classic()+labs(fill="Phenotype")+ annotate(geom="label", x = mean(caddat2$wgrs)-1, y = 0.7, label = round(mean(caddat2$wgrs),digits = 2),fill="#00BA38")+ annotate(geom="label", x = mean(bcdat2$wgrs)-1, y = 0.8, label = round(mean(bcdat2$wgrs),digits = 2),fill="#F8766D")+annotate(geom="label", x = mean(pcdat2$wgrs)+1, y = 0.5, label = round(mean(pcdat2$wgrs),digits = 2),fill="#619cff")+geom_vline(xintercept = mean(pcdat2$wgrs), linetype="longdash",color = "black", alpha=0.8)+geom_vline(xintercept = mean(caddat2$wgrs), linetype="longdash",color = "black", alpha=0.8)+geom_vline(xintercept = mean(bcdat2$wgrs), linetype="longdash",color = "black",alpha=0.8)+scale_y_continuous(expand = c(0, 0))

# Phenotypes
#sex
sex<-fread("UKBB_182_final.fam")
sex<-sex %>% select(1,5)
colnames(sex)<-c("f.eid","sex")
t2<-sex[match(caddat2$sample,sex$f.eid),]
caddat2$sex<-t2$sex
#age
age<-as.data.frame(age)
t2<-age[match(caddat2$sample,age$f.eid),]
caddat2$age<-t2$age
#BMI
t2<-bmi[match(caddat2$sample,bmi$f.eid),]
caddat2$bmi<-t2$bmi
#blood pressure
t2<-bp[match(caddat2$sample,bp$f.eid),]
caddat2$BloodPressure<-t2$BloodPressure
#LDL
t2<-ldl[match(caddat2$sample,ldl$f.eid),]
caddat2$LDL<-t2$LDL
#LPA
t2<-lpa[match(caddat2$sample,lpa$f.eid),]
caddat2$LPA<-t2$LPA
#Diabets
t2<-dia2[match(caddat2$sample,dia2$f.eid),]
caddat2$Diabetes<-t2$Diabetes
#Smoking
t2<-smoke[match(bcdat2$sample,smoke$f.eid),]
bcdat2$smoke<-t2$smoke
#medication
meds<-fread("6153_medS_QC.ukb")
medsa<-fread("6153_medS.sample")
#lipid-medication
lpm<-meds[meds$value==1,]
t2<-lpm[match(medsa$f.eid,lpm$f.eid),]
medsa$lpm<-t2$value
medsa$lpm[is.na(medsa$lpm)]<-0
#blood pressure medication
bpm<-meds[meds$value==2,]
t2<-bpm[match(medsa$f.eid,bpm$f.eid),]
medsa$bpm<-t2$value
medsa$bpm[is.na(medsa$bpm)]<-0
#diabetes medication
#insulin hormone oral
diam<-meds[(meds$value==3 | meds$value==4 | meds$value==5),]
diam$value<-3
t2<-diam[match(medsa$f.eid,diam$f.eid),]
medsa$diam<-t2$value
medsa$diam[is.na(medsa$diam)]<-0
t2<-medsa[match(caddat2$sample,medsa$f.eid),]
caddat2$lpm<-t2$lpm
caddat2$bpm<-t2$bpm
caddat2$diam<-t2$diam
#FH+
fh<-fread("ukb_fh.sample")
fh$fh<-1
fh<-as.data.frame(fh)
fha<-fread("ukb_wes.fam")
fha<-as.data.frame(fha)
t2<-fh[match(fha$V1,fh$sample),]
t2$fh[is.na(t2$fh)]<-0
fha$fh<-t2$fh
t2<-fha[match(caddat2$sample,fha$V1),]
caddat2$fh<-t2$fh
# convert mmol/L to mg/dl
caddat2$ldl2<-caddat2$LDL*18.1082
#Exercise MET:7.5;21
t2<-exmet[match(caddat2$sample,exmet$f.eid),]
caddat2$MET<-t2$MET
#BC~alchol 
t2<-alchol_freq[match(bcdat2$sample,alchol_freq$f.eid),]
bcdat2$alchol_freq<-t2$freq
bcdat2$alchol_dosage<-t2$dosage
bcdat2$alchol_freq_class<-t2$freq_class
#BC~age
t2<-age[match(bcdat2$sample,age$f.eid),]
bcdat2$age<-t2$age
# statistics deciles
statgrs<-function(dat2l,stat){
  stat<-data.frame("age"=0,"agesd"=0,"sex"=0,"bmi"=0,"bmisd"=0,"bp"=0,"bpsd"=0,"ldl"=0,"ldlsd"=0,"lpa"=0,"lpasd"=0,"dia"=0,"smok"=0,"lpm"=0,"bpm"=0,"diam"=0,"fh"=0,"met"=0,"metsd"=0)
  for (i in 1:10){
    t<-dat2l[[i]]
    age<-mean(t$age,na.rm=TRUE)
    agesd<-sd(t$age,na.rm=TRUE)
    sex<-sum(t$sex)-nrow(t)#1male 2female - female
    bmi<-mean(t$bmi,na.rm=TRUE)
    bmisd<-sd(t$bmi,na.rm=TRUE)
    bp<-mean(t$BloodPressure,na.rm=TRUE)
    bpsd<-sd(t$BloodPressure,na.rm=TRUE)
    ldl<-mean(t$LDL,na.rm=TRUE)
    ldlsd<-sd(t$LDL,na.rm=TRUE)
    #ldl<-mean(t$ldl2,na.rm=TRUE)
    #ldlsd<-sd(t$ldl2,na.rm=TRUE)
    lpa<-mean(t$LPA,na.rm=TRUE)
    lpasd<-sd(t$LPA,na.rm=TRUE)
    dia<-sum(t$Diabetes,na.rm=TRUE)
    smok<-sum(t$smoke,na.rm=TRUE)
    lpm<-sum(t$lpm,na.rm=TRUE)
    bpm<-sum(t$bpm,na.rm=TRUE)
    diam<-sum(t$diam,na.rm=TRUE)
    fh<-sum(t$fh,na.rm=TRUE)
    met<-mean(t$MET,na.rm=TRUE)
    metsd<-sd(t$MET,na.rm=TRUE)
    stat[i,]<-c(age,agesd,sex,bmi,bmisd,bp,bpsd,ldl,ldlsd,lpa,lpasd,dia,smok,lpm,bpm,diam,fh,met,metsd)
  }
  return(stat)
}

# uGRS
dat2l<-split(caddat2,caddat2$uorder)
ustat<-statgrs(dat2l,ustat)
# wGRS
dat2l<-split(caddat2,caddat2$worder)
wstat<-statgrs(dat2l,wstat)

# smoking increased prevelence percentage
#prevstat<-function(dat2l,adres){
#  adres<-data.frame("smpre"=0,"nsmpre"=0,"diff_sm"=0,"bmihp"=0,"bmilp"=0,"diff_bmi"=0,"tn1p"=0,"tn2p"=0,"tn3p"=0,"diapre"=0,"dianpre"=0,"alc0pre"=0,"alc1pre"=0,"alc2pre"=0)
#  for (i in 1:10){
#    t<-dat2l[[i]]
#    tn<-t[is.na(t$smoke)==FALSE,]
#    sm<-tn[tn$smoke==1,]
#    nsm<-tn[tn$smoke==0,]
#    smpre<-nrow(sm[sm$CC==1,])/nrow(sm)
#    nsmpre<-nrow(nsm[nsm$CC==1,])/nrow(nsm)
#    diff_sm<-smpre-nsmpre
#    tn<-t[is.na(t$bmi)==FALSE,]
#    bmih<-tn[tn$bmi>=30,]
#    bmil<-tn[tn$bmi<30,]
#    bmihp<-nrow(bmih[bmih$CC==1,])/nrow(bmih)
#    bmilp<-nrow(bmil[bmil$CC==1,])/nrow(bmil)
#    diff_bmi<-bmihp-bmilp
#    tn<-t[is.na(t$MET)==FALSE,]
#    tn1<-tn[tn$MET<=7.5,]
#    tn2<-tn[(tn$MET>7.5 & tn$MET<=21),]
#    tn3<-tn[tn$MET>21,]
#    tn1p<-nrow(tn1[tn1$CC==1,])/nrow(tn1)
#    tn2p<-nrow(tn2[tn2$CC==1,])/nrow(tn2)
#    tn3p<-nrow(tn3[tn3$CC==1,])/nrow(tn3)
#    tn<-t[is.na(t$Diabetes)==FALSE,]
#    diay<-tn[tn$Diabetes==1,]
#    dian<-tn[tn$Diabetes==0,]
#    diapre<-nrow(diay[diay$CC==1,])/nrow(diay)
#    dianpre<-nrow(dian[dian$CC==1,])/nrow(dian)
#    diff_dia<-diapre-dianpre
#    tn<-t[is.na(t$alchol_freq_class)==FALSE,]
#    alc0<-tn[tn$alchol_freq_class==0,]
#    alc1<-tn[tn$alchol_freq_class==1,]
#    alc2<-tn[tn$alchol_freq_class==2,]
#    alc0pre<-nrow(alc0[alc0$CC==1,])/nrow(alc0)
#    alc1pre<-nrow(alc1[alc1$CC==1,])/nrow(alc1)
#    alc2pre<-nrow(alc2[alc2$CC==1,])/nrow(alc2)
#    adres[i,]<-c(smpre,nsmpre,diff_sm,bmihp,bmilp,diff_bmi,tn1p,tn2p,tn3p,diapre,dianpre,alc0pre,alc1pre,alc2pre)
#  }
#  return(adres)
#}

#Figure 2 version trait ~ pheno
# stat case, all, nra
statsub<-function(x,pheno){
  ressub<-data.frame("case"=0,"all"=0,"nra"=0,"pheno"=0)
  case<-sum(x$CC)
  all<-nrow(x)
  nra<-mean(x$NRA)
  ressub[1,]<-c(case,all,nra,pheno)
  return(ressub)
}
statsub<-function(x,pheno){
  ressub<-data.frame("case"=0,"all"=0,"nra"=0,"pheno"=0)
  case<-sum(x$CC)
  all<-nrow(x)
  nra<-mean(x$wgrs)
  ressub[1,]<-c(case,all,nra,pheno)
  return(ressub)
}
prevstatn<-function(dat2l,adres){
  adres<-data.frame("case"=0,"all"=0,"nra"=0,"pheno"=0,"decile"=0)
  for (i in 1:10){
    t<-dat2l[[i]]
    tn<-t[is.na(t$smoke)==FALSE,]
    sm<-tn[tn$smoke==1,]
    nsm<-tn[tn$smoke==0,]
    sm0<-statsub(sm,"smoke")
    nsm0<-statsub(nsm,"non-smoke")
    tn<-t[is.na(t$bmi)==FALSE,]
    bmih<-tn[tn$bmi>=30,]
    bmil<-tn[tn$bmi<30,]
    bmih0<-statsub(bmih,"bmi>=30")
    bmil0<-statsub(bmil,"bmi<30")
    tn<-t[is.na(t$sex)==FALSE,]
    sex1<-tn[tn$sex==1,]
    sex2<-tn[tn$sex==2,]
    sex10<-statsub(sex1,"male")
    sex20<-statsub(sex2,"female")
    tn<-t[is.na(t$MET)==FALSE,]
    tn1<-tn[tn$MET<=7.5,]
    tn2<-tn[(tn$MET>7.5 & tn$MET<=21),]
    tn3<-tn[tn$MET>21,]
    met1<-statsub(tn1,"MET<=7.5")
    met2<-statsub(tn2,"7.5<MET<=21")
    met3<-statsub(tn3,"21<MET")
    #tn1<-tn[tn$MET<3.75,]
    #tn2<-tn[tn$MET>=3.75,]
    #met1<-statsub(tn1,"low exercise frequancy")
    #met2<-statsub(tn2,"high exercise frequancy")
    tn<-t[is.na(t$Diabetes)==FALSE,]
    diay<-tn[tn$Diabetes==1,]
    dian<-tn[tn$Diabetes==0,]
    diay0<-statsub(diay,"diabetes")
    dian0<-statsub(dian,"non-diabetes")
    adres0<-rbind(sm0,nsm0,sex10,sex20,bmih0,bmil0,met1,met2,met3,diay0,dian0)
    adres0$decile<-i
    adres<-rbind(adres,adres0)
  }
  return(adres)
}
# uGRS
dat2l<-split(caddat2,caddat2$uorder)
#uadres<-prevstat(dat2l,uadres)
uadres<-prevstatn(dat2l,uadres)
t<-caddat2
adres0$decile<-"average"
uadres<-rbind(uadres,adres0)
# wGRS
dat2l<-split(caddat2,caddat2$worder)
#wadres<-prevstat(dat2l,wadres)
wadres<-prevstatn(dat2l,wadres)
# breast cancer
prevstatn<-function(dat2l,adres){
  adres<-data.frame("case"=0,"all"=0,"nra"=0,"pheno"=0,"decile"=0)
  for (i in 1:10){
    t<-dat2l[[i]]
    tn<-t[is.na(t$bmi)==FALSE,]
    bmih<-tn[tn$bmi>=30,]
    bmil<-tn[tn$bmi<30,]
    bmih0<-statsub(bmih,"bmi>=30")
    bmil0<-statsub(bmil,"bmi<30") #26.32
    tn<-t[is.na(t$alchol_freq_class)==FALSE,]
    alc0<-tn[tn$alchol_freq_class==0,]
    alc1<-tn[tn$alchol_freq_class==1,]
    alc2<-tn[tn$alchol_freq_class==2,]
    alc00<-statsub(alc0,"never")
    alc10<-statsub(alc1,"once")
    alc20<-statsub(alc2,"frequent") #9
    #tn<-t[is.na(t$fh)==FALSE,]
    #fhm<-tn[tn$fh==1,]
    #fhnm<-tn[tn$fh==0,]
    #fhm0<-statsub(fhm,"FH mutation")
    #fhnm0<-statsub(fhnm,"Non FH mutation")
    adres0<-rbind(bmih0,bmil0,alc00,alc10,alc20)
    adres0$decile<-i
    adres<-rbind(adres,adres0)
  }
  return(adres)
}
# uGRS
dat2l<-split(bcdat2,bcdat2$uorder)
#uadres<-prevstat(dat2l,uadres)
uadres<-prevstatn(dat2l,uadres)
t<-bcdat2
adres0$decile<-"average"
uadres<-rbind(uadres,adres0)
# wGRS
dat2l<-split(bcdat2,bcdat2$worder)
#wadres<-prevstat(dat2l,wadres)
wadres<-prevstatn(dat2l,wadres)


library(varhandle)
#barplot smoking
smot<-matrix(c("Never_Smoked","Low risk",1.63,"Ever_Smoked","Low risk",3.65,"Never_Smoked","Medium risk",3.21,"Ever_Smoked","Medium risk",6.61,"Never_Smoked","High risk",6.08,"Ever_Smoked","High risk",12.11),ncol=3,byrow=TRUE)
smot<-as.data.frame(smot)
smot<-unfactor(smot)
smot$V2 <- factor(smot$V2, levels = c("Low risk","Medium risk","High risk"))
smot$V1 <- factor(smot$V1, levels = c("Never_Smoked","Ever_Smoked"))
diff_df = smot %>%
  group_by(V2) %>%
  spread(V1, V3) %>%
  mutate(diff = Ever_Smoked - Never_Smoked,
         max_y = max(Ever_Smoked, Never_Smoked))
ggplot(smot,aes(V2,V3))+
  geom_bar(aes(y = max_y), data = diff_df, stat = "identity", fill = "grey80", width = 0.4)+
  geom_bar(aes(fill=V1),stat="identity",colour="black",position="dodge")+
  geom_text(aes(label = diff, y = max_y), vjust=-0.5, data = diff_df,hjust = 1, colour = scales::muted("red")) +
  theme_classic()+xlab("uGRS categories")+ylab("Adjusted CAD prevalence (%)")+ theme(legend.title =element_blank(),legend.position="top")+
  scale_fill_manual(values=c("#33A2C9","#CD375D"))

#barplot FH
smot<-matrix(c("Non_mutated","Low risk",1.88,"Mutated","Low risk",4.17,"Non_mutated","Medium risk",4.01,"Mutated","Medium risk",8.29,"Non_mutated","High risk",7.33,"Mutated","High risk",24),ncol=3,byrow=TRUE)
smot<-as.data.frame(smot)
smot<-unfactor(smot)
smot$V2 <- factor(smot$V2, levels = c("Low risk","Medium risk","High risk"))
smot$V1 <- factor(smot$V1, levels = c("Non_mutated","Mutated"))
diff_df = smot %>%
  group_by(V2) %>%
  spread(V1, V3) %>%
  mutate(diff = Mutated - Non_mutated,
         max_y = max(Mutated, Non_mutated))
ggplot(smot,aes(V2,V3))+
  geom_bar(aes(y = max_y), data = diff_df, stat = "identity", fill = "grey80", width = 0.4)+
  geom_bar(aes(fill=V1),stat="identity",colour="black",position="dodge")+
  geom_text(aes(label = diff, y = max_y), vjust=-0.5, data = diff_df,hjust = 1, colour = scales::muted("red")) +
  theme_classic()+xlab("uGRS categories")+ylab("Adjusted CAD prevalence (%)")+ theme(legend.title =element_blank(),legend.position="top")+
  scale_fill_manual(values=c("#33A2C9","#CD375D"))

#barplot BMI
smot<-matrix(c("BMI<30","Low risk",3.59,"BMI>=30","Low risk",3.66,"BMI<30","Medium risk",6.37,"BMI>=30","Medium risk",6.65,"BMI<30","High risk",10.25,"BMI>=30","High risk",11.09),ncol=3,byrow=TRUE)
smot<-as.data.frame(smot)
smot<-unfactor(smot)
smot$V2 <- factor(smot$V2, levels = c("Low risk","Medium risk","High risk"))
smot$V1 <- factor(smot$V1, levels = c("BMI<30","BMI>=30"))
diff_df = smot %>%
  group_by(V2) %>%
  spread(V1, V3) %>%
  mutate(diff = `BMI>=30` - `BMI<30`,
         max_y = max(`BMI>=30`, `BMI<30`))
ggplot(smot,aes(V2,V3))+
  geom_bar(aes(y = max_y), data = diff_df, stat = "identity", fill = "grey80", width = 0.4)+
  geom_bar(aes(fill=V1),stat="identity",colour="black",position="dodge")+
  geom_text(aes(label = round(diff,digits=2), y = max_y), vjust=-0.5, data = diff_df,hjust = 1, colour = scales::muted("red")) +
    theme_classic()+xlab("GRS categories")+ylab("Adjusted BC prevalence (%)")+ theme(legend.title =element_blank(),legend.position="top")+
  scale_fill_manual(values=c("#33A2C9","#CD375D"))

#barplot MET
smot<-matrix(c("MET<7.5","Low risk",2.20,"7.5<=MET<21","Low risk",1.82,"MET>=21","Low risk",1.60,"MET<7.5","Medium risk",5.10,"7.5<=MET<21","Medium risk",3.95,"MET>=21","Medium risk",3.31,"MET<7.5","High risk",8.74,"7.5<=MET<21","High risk",7.67,"MET>=21","High risk",6.45),ncol=3,byrow=TRUE)
smot<-as.data.frame(smot)
smot<-unfactor(smot)
smot$V2 <- factor(smot$V2, levels = c("Low risk","Medium risk","High risk"))
smot$V1 <- factor(smot$V1, levels = c("MET>=21","7.5<=MET<21","MET<7.5"))
diff_df = smot %>%
  group_by(V2) %>%
  spread(V1, V3) %>%
  mutate(diff = `MET<7.5` - `MET>=21`,
         max_y = max(`MET<7.5`, `MET>=21`))
ggplot(smot,aes(V2,V3))+
  geom_bar(aes(y = max_y), data = diff_df, stat = "identity", fill = "grey80", width = 0.4)+
  geom_bar(aes(fill=V1),stat="identity",colour="black",position="dodge")+
  geom_text(aes(label = diff, y = max_y), vjust=-0.5, data = diff_df,hjust = 1, colour = scales::muted("red")) +
  theme_classic()+xlab("uGRS categories")+ylab("Adjusted CAD prevalence (%)")+ theme(legend.title =element_blank(),legend.position="top")+
  scale_fill_manual(values=c("#33A2C9","#CD375D","#EEB422"))

smot<-matrix(c("MET<7.5","Low risk",2.81,"7.5<=MET<21","Low risk",2.23,"MET>=21","Low risk",1.69,"MET<7.5","Medium risk",5.11,"7.5<=MET<21","Medium risk",3.90,"MET>=21","Medium risk",3.30,"MET<7.5","High risk",8.95,"7.5<=MET<21","High risk",7.73,"MET>=21","High risk",6.43),ncol=3,byrow=TRUE)
smot<-as.data.frame(smot)
smot<-unfactor(smot)
smot$V2 <- factor(smot$V2, levels = c("Low risk","Medium risk","High risk"))
smot$V1 <- factor(smot$V1, levels = c("MET>=21","7.5<=MET<21","MET<7.5"))
diff_df = smot %>%
  group_by(V2) %>%
  spread(V1, V3) %>%
  mutate(diff = `MET<7.5` - `MET>=21`,
         max_y = max(`MET<7.5`, `MET>=21`))
ggplot(smot,aes(V2,V3))+
  geom_bar(aes(y = max_y), data = diff_df, stat = "identity", fill = "grey80", width = 0.4)+
  geom_bar(aes(fill=V1),stat="identity",colour="black",position="dodge")+
  geom_text(aes(label = diff, y = max_y), vjust=-0.5, data = diff_df,hjust = 1, colour = scales::muted("red")) +
  theme_classic()+xlab("wGRS categories")+ylab("Adjusted CAD prevalence (%)")+ theme(legend.title =element_blank(),legend.position="top")+
  scale_fill_manual(values=c("#33A2C9","#CD375D","#EEB422"))

# CAD physical activity
smot<-matrix(c("MET<7.5","Low risk",1,"7.5<=MET<21","Low risk",0.84,"MET>=21","Low risk",0.62,"MET<7.5","Medium risk",1.98,"7.5<=MET<21","Medium risk",1.61,"MET>=21","Medium risk",1.46,"MET<7.5","High risk",3.96,"7.5<=MET<21","High risk",3.23,"MET>=21","High risk",2.88),ncol=3,byrow=TRUE)
smot<-matrix(c("MET<7.5","Low risk",1,"7.5<=MET<21","Low risk",0.75,"MET>=21","Low risk",0.60,"MET<7.5","Medium risk",1.99,"7.5<=MET<21","Medium risk",1.69,"MET>=21","Medium risk",1.26,"MET<7.5","High risk",3.98,"7.5<=MET<21","High risk",3.27,"MET>=21","High risk",2.78),ncol=3,byrow=TRUE)
smot<-matrix(c("MET<7.5","Low risk",2.2,"7.5<=MET<21","Low risk",1.82,"MET>=21","Low risk",1.60,"MET<7.5","Medium risk",5.1,"7.5<=MET<21","Medium risk",3.95,"MET>=21","Medium risk",3.31,"MET<7.5","High risk",8.74,"7.5<=MET<21","High risk",7.67,"MET>=21","High risk",6.45),ncol=3,byrow=TRUE)
smot<-matrix(c("MET<7.5","Low risk",2.81,"7.5<=MET<21","Low risk",2.23,"MET>=21","Low risk",1.69,"MET<7.5","Medium risk",5.11,"7.5<=MET<21","Medium risk",3.90,"MET>=21","Medium risk",3.30,"MET<7.5","High risk",8.95,"7.5<=MET<21","High risk",7.73,"MET>=21","High risk",6.43),ncol=3,byrow=TRUE)
#barplot BMI T2D
library(varhandle)
smot<-matrix(c("BMI<30","Low risk",2.87,"BMI>=30","Low risk",16.29,"BMI<30","Medium risk",6.32,"BMI>=30","Medium risk",22.41,"BMI<30","High risk",11.69,"BMI>=30","High risk",33.12),ncol=3,byrow=TRUE)
smot<-matrix(c("BMI<30","Low risk",3.41,"BMI>=30","Low risk",3.79,"BMI<30","Medium risk",6.18,"BMI>=30","Medium risk",6.57,"BMI<30","High risk",10.42,"BMI>=30","High risk",11.21),ncol=3,byrow=TRUE)
#BC ALCOHOL
smot<-matrix(c("Class I","Low risk",3.63,"Class II","Low risk",3.50,"Class III","Low risk",3.66,"Class I","Medium risk",5.97,"Class II","Medium risk",6.25,"Class III","Medium risk",6.51,"Class I","High risk",9.61,"Class II","High risk",10.83,"Class III","High risk",11.10),ncol=3,byrow=TRUE)
smot<-as.data.frame(smot)
smot<-unfactor(smot)
smot$V2 <- factor(smot$V2, levels = c("Low risk","Medium risk","High risk"))
smot$V1 <- factor(smot$V1, levels = c("Class I","Class II","Class III"))
diff_df = smot %>%
  group_by(V2) %>%
  spread(V1, V3) %>%
  mutate(diff = `Class III` - `Class I`,
         max_y = max(`Class III`, `Class I`))
ggplot(smot,aes(V2,V3))+
  geom_bar(aes(y = max_y), data = diff_df, stat = "identity", fill = "grey80", width = 0.4)+
  geom_bar(aes(fill=V1),stat="identity",colour="black",position="dodge")+
  geom_text(aes(label = round(diff,digits=2), y = max_y), vjust=-0.5, data = diff_df,hjust = 1, colour = scales::muted("red")) +
  theme_classic()+xlab("uGRS categories")+ylab("Adjusted BC prevalence (%)")+ theme(legend.title =element_blank(),legend.position="top")+
  scale_fill_manual(values=c("#33A2C9","#CD375D","#EEB422"))

