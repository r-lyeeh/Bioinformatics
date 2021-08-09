## PRS calculation
## plink --bfile ${name} --score 2 3 4 sum --out ${name}
## PRsice2.sh

# number of allele count 
analysis<-function(dat,pheno,samplesize){
  if(missing(samplesize)){
    samplesize=200
  } else {
    samplesize
  }
  
  modComp  <- summary( glm(CC~NRA,family="binomial",data=dat) )$coefficients
  dat$grp  <- round(dat$NRA) #cut(dat$NRA,breaks = seq(min(dat$NRA),max(dat$NRA),by=2),include.lowest=TRUE)
  prev_p   <- aggregate(dat$CC,by=list(dat$grp),FUN=mean,na.rm=T)
  prev_n   <- aggregate(dat$CC,by=list(dat$grp),FUN=sum,na.rm=T)
  prev_N   <- aggregate(dat$CC,by=list(dat$grp),FUN=function(x) length(!is.na(x)))
  allele_c <- aggregate(dat$NRA,by=list(dat$grp),FUN=mean)
  d        <- cbind.data.frame(N=prev_N$x,n=prev_n$x,p=prev_p$x,numAllele=allele_c$x)
  d$p_se   <- sqrt(d$p*(1-d$p)/d$N)
  d        <- d[which(d$N>=samplesize),]
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
wanalysis<-function(dat,pheno,samplesize){
  if(missing(samplesize)){
    samplesize=200
  } else {
    samplesize
  }
  
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


## Produce figures
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
cairo_pdf(filename='BC_wGRS_per_01.pdf', width=11.06, height=12.39,pointsize = 12,family="Univers")
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
#558x306
#Figure 1 histogram case controls
#nf <- layout(matrix(c(1,3,2,4),ncol=2), width=c(558,558,558,558),height=c(306,306,306,306),TRUE)
cairo_pdf(filename='F1B_BC.pdf', width=5.58, height=3.06,pointsize = 12,family="Univers")
dev.off()

# normalized
ggplot(data=caddat2,aes(x=NRA,group=CC,fill=CC))+geom_histogram(aes(y=..density..*406996),binwidth = 1,color="white",size=0.2,alpha=0.7,position="identity")+xlab("Number of risk alleles")+ylab("Number of samples")+theme_classic()+labs(fill="Phenotype")+ annotate(geom="label", x = mean(caddat2[caddat2$CC==0,]$NRA)-10, y = 18000, label = round(mean(caddat2[caddat2$CC==0,]$NRA),digits = 1),fill="#42AD60")+ annotate(geom="label", x = mean(caddat2[caddat2$CC==1,]$NRA)+20, y = 4000, label = round(mean(caddat2[caddat2$CC==1,]$NRA),digits = 1),fill="#6DAFD8")+scale_fill_manual(values=c("#42AD60","#6DAFD8"),labels = c("control","case"))+ geom_vline(xintercept = mean(caddat2[caddat2$CC==1,]$NRA), linetype="longdash",color = "black", alpha=0.8)+geom_vline(xintercept = mean(caddat2[caddat2$CC==0,]$NRA), linetype="longdash",color = "black",alpha=0.8)+ ggtitle("UKBB CAD")+scale_y_continuous(name="Number of samples (control)",sec.axis=sec_axis(~./20.03919,name="Number of samples (case)"),expand = c(0, 0))
#raw data
# uGRS
ggplot(data=caddat2)+geom_histogram(aes(x=NRA,group=CC,fill=CC),binwidth = 1,color="white",size=0.2,position="identity")+xlab("Number of risk alleles")+ylab("Number of samples")+theme_classic()+labs(fill="Phenotype")+ annotate(geom="label", x = mean(caddat2[caddat2$CC==0,]$NRA)-10, y = 18000, label = round(mean(caddat2[caddat2$CC==0,]$NRA),digits = 2),fill="#42AD60")+ annotate(geom="label", x = mean(caddat2[caddat2$CC==1,]$NRA)+10, y = 3000, label = round(mean(caddat2[caddat2$CC==1,]$NRA),digits = 2),fill="#6DAFD8")+scale_fill_manual(values=c("#42AD60","#6DAFD8"),labels = c("control","case"))+ geom_vline(xintercept = mean(caddat2[caddat2$CC==1,]$NRA), linetype="longdash",color = "black", alpha=0.8)+geom_vline(xintercept = mean(caddat2[caddat2$CC==0,]$NRA), linetype="longdash",color = "black",alpha=0.8)+scale_y_continuous(expand = c(0, 0))+ ggtitle("UKBB_CAD")
#wGRS
ggplot(data=caddat2)+geom_histogram(aes(x=wgrs,group=CC,fill=CC),binwidth = 0.1,color="white",size=0.2,position="identity")+xlab("Number of weighted GRS")+ylab("Number of samples")+theme_classic()+labs(fill="Phenotype")+ annotate(geom="label", x = mean(caddat2[caddat2$CC==0,]$wgrs)-1, y = 18000, label = round(mean(caddat2[caddat2$CC==0,]$wgrs),digits = 2),fill="#42AD60")+ annotate(geom="label", x = mean(caddat2[caddat2$CC==1,]$wgrs)+1, y = 3000, label = round(mean(caddat2[caddat2$CC==1,]$wgrs),digits = 2),fill="#6DAFD8")+scale_fill_manual(values=c("#42AD60","#6DAFD8"),labels = c("control","case"))+ geom_vline(xintercept = mean(caddat2[caddat2$CC==1,]$wgrs), linetype="longdash",color = "black", alpha=0.8)+geom_vline(xintercept = mean(caddat2[caddat2$CC==0,]$wgrs), linetype="longdash",color = "black",alpha=0.8)+scale_y_continuous(expand = c(0, 0))+ ggtitle("UKBB_CAD")

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

library(varhandle)
#barplot smoking
#smot<-matrix(c("Never_Smoked","Low risk",1.63,"Ever_Smoked","Low risk",3.65,"Never_Smoked","Medium risk",3.21,"Ever_Smoked","Medium risk",6.61,"Never_Smoked","High risk",6.08,"Ever_Smoked","High risk",12.11),ncol=3,byrow=TRUE)
#smot<-as.data.frame(smot)
#smot<-unfactor(smot)
#smot$V2 <- factor(smot$V2, levels = c("Low risk","Medium risk","High risk"))
#smot$V1 <- factor(smot$V1, levels = c("Never_Smoked","Ever_Smoked"))
#diff_df = smot %>%
#  group_by(V2) %>%
#  spread(V1, V3) %>%
#  mutate(diff = Ever_Smoked - Never_Smoked,
#         max_y = max(Ever_Smoked, Never_Smoked))
#ggplot(smot,aes(V2,V3))+
#  geom_bar(aes(y = max_y), data = diff_df, stat = "identity", fill = "grey80", width = 0.4)+
#  geom_bar(aes(fill=V1),stat="identity",colour="black",position="dodge")+
#  geom_text(aes(label = diff, y = max_y), vjust=-0.5, data = diff_df,hjust = 1, colour = scales::muted("red")) +
#  theme_classic()+xlab("uGRS categories")+ylab("Adjusted CAD prevalence (%)")+ theme(legend.title =element_blank(),legend.position="top")+
#  scale_fill_manual(values=c("#33A2C9","#CD375D"))

# Figure 3
library(gridExtra)
library(egg)

plotor<-function(lipids_ugrs,title0){
  ggplot(lipids_ugrs, aes(x=as.numeric(as.character(lipids_ugrs$row.names)), y=Odds.Ratio,grpup=1)) + 
    geom_line() +
    geom_point()+
    geom_errorbar(aes(ymin=Odds.Ratio.CI.L, ymax=Odds.Ratio.CI.U), width=.2,
                  position=position_dodge(0.05))+theme_classic()+labs(title=title0,y="Odds Ratio")+scale_x_discrete(name="Decile", limits=c("1","2","3","4","5","6","7","8","9","10"))
}
p1<-plotor(lipids_ugrs,"Lipids risk allele count")
p2<-plotor(lipids_wgrs,"Lipids wGRS")
p3<-plotor(nonlipids_ugrs,"Non-lipids risk allele count")
p4<-plotor(nonlipids_wgrs,"Non-lipids wGRS")
ggarrange(p1, p2, p3,p4, 
          labels = c("A", "B", "C","D"),
          ncol = 2, nrow = 2)

uGRS<-read.table("uGRS_raw",header=TRUE,row.names = NULL,sep="\t")
uGRS$CIL<-uGRS$Prevalence-100*(uGRS$p_se)
uGRS$CIU<-uGRS$Prevalence+100*(uGRS$p_se)
t1<-uGRS[uGRS$Align=="Main score",]

uGRSdecile<-read.table("uGRS_decile",sep="\t",header=TRUE)
uGRSdecile$CIL<-uGRSdecile$p-uGRSdecile$se
uGRSdecile$CIU<-uGRSdecile$p+uGRSdecile$se
ggplot(uGRSdecile, aes(x=nra, y=p*100,grpup=Class,color=Class)) + 
  geom_line(size=1.8) +
  geom_point(size=4)+theme_classic()+labs(title="10th Decile",x="number of risk alleles",y="Prevalence (%)")+theme(text = element_text(size=20)) + labs(color='Class')+geom_errorbar(aes(ymin=CIL*100, ymax=CIU*100), size=1,width=.2,position=position_dodge(0.05),color="steelblue")+scale_y_continuous(limits=c(0,10))

uGRSs<-read.table("uGRS_3",sep="\t",header=TRUE)
uGRSs$CIL<-uGRSs$p-uGRSs$SE
uGRSs$CIU<-uGRSs$p+uGRSs$SE
p2<-ggplot(uGRSs, aes(x=class, y=p*100,group=class,color=class)) + 
  geom_point(size=4)+theme_classic()+theme(text = element_text(size=15))+geom_errorbar(aes(ymin=CIL*100, ymax=CIU*100), size=1,width=.1,position=position_dodge(0.05),color="steelblue")+scale_x_discrete(name="10th Decile",labels=c("Non-lipid related SNPs","Lipid related SNPs"))+scale_y_continuous(limits=c(0,10),position="right")+labs(x="number of risk alleles",y="Prevalence (%)")+theme(legend.position = "none")+
  theme(text = element_text(size=20)) 


p1<-ggplot(t1, aes(x=as.numeric(as.character(Decile)), y=Prevalence,grpup=group,color=group)) + 
  geom_line(size=1) +
  geom_point(size=2)+theme_classic()+labs(y="Prevalence (%)")+scale_x_discrete(name="Decile", limits=c("1","2","3","4","5","6","7","8","9","10"))+theme(text = element_text(size=16)) + labs(color='Group')+geom_errorbar(aes(ymin=CIL, ymax=CIU), size=1,width=.2,position=position_dodge(0.05),color="steelblue")+scale_y_continuous(limits=c(0,10))+theme(legend.position = "none")

ggarrange(p1, p2, 
          labels = c("A", "B"),
          ncol = 2, nrow = 1)

grid.arrange(p1,p2,ncol=2)


# Decile; PRS line shift figure
caddat2<-fread("caddat2.res")
decilep<-function(dat){
  modComp  <- summary( glm(CC~NRA,family="binomial",data=dat) )$coefficients
  dat$grp  <- round(dat$NRA) #cut(dat$NRA,breaks = seq(min(dat$NRA),max(dat$NRA),by=2),include.lowest=TRUE)
  prev_p   <- aggregate(dat$CC,by=list(dat$grp),FUN=mean,na.rm=T)
  prev_n   <- aggregate(dat$CC,by=list(dat$grp),FUN=sum,na.rm=T)
  prev_N   <- aggregate(dat$CC,by=list(dat$grp),FUN=function(x) length(!is.na(x)))
  allele_c <- aggregate(dat$NRA,by=list(dat$grp),FUN=mean)
  d        <- cbind.data.frame(N=prev_N$x,n=prev_n$x,p=prev_p$x,numAllele=allele_c$x)
  d$p_se   <- sqrt(d$p*(1-d$p)/d$N)
  d        <- d[which(d$N>=100),]
  d$px     <- d$p / d$numAllele
  d$px_se  <- d$p_se / d$numAllele
  t<-dat[dat$NRA>=161]
  d<-rbind(c(0,0,0,161,0,0,0),d)
  
  modelLog    <- glm(cbind(n,N)~numAllele,data=d,family = binomial(link="log"))
  predLog     <- predict(modelLog,type = "response")
  d$predLog<-predLog
  return(d)
}

decilepoint<-function(dat){
  quant<-quantile(dat$NRA,probs = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
  dat$uorder<-10
  dat[dat$NRA<=quant[1],5]<-1
  dat[quant[1]<dat$NRA & dat$NRA<=quant[2],5]<-2
  dat[quant[2]<dat$NRA & dat$NRA<=quant[3],5]<-3
  dat[quant[3]<dat$NRA & dat$NRA<=quant[4],5]<-4
  dat[quant[4]<dat$NRA & dat$NRA<=quant[5],5]<-5
  dat[quant[5]<dat$NRA & dat$NRA<=quant[6],5]<-6
  dat[quant[6]<dat$NRA & dat$NRA<=quant[7],5]<-7
  dat[quant[7]<dat$NRA & dat$NRA<=quant[8],5]<-8
  dat[quant[8]<dat$NRA & dat$NRA<=quant[9],5]<-9
  dat2l<-split(dat,dat$uorder)
  q<-data.frame("N"=0,"n"=0,"p"=0,"numAllele"=0,"p_se"=0)
  for (i in 1:10){
    t<-dat2l[[i]]
    q[i,1:4]<-c(nrow(t),sum(t$CC),sum(t$CC)/nrow(t),mean(t$NRA))
    q$p_se[i]<-sqrt(q$p[i]*(1-q$p[i])/q$N[i])
  }
  modelLog    <- glm(cbind(n,N)~numAllele,data=q,family = binomial(link="log"))
  predLog     <- predict(modelLog,type = "response")
  q$predLog<-predLog
  return(q)
}


[100-250]
caddatm<-caddat2[caddat2$sex==1,]
caddatm<-caddatm[caddatm$smoke==0,]
caddatm<-caddatm[caddatm$Diabetes==0,]
caddatm<-caddatm[caddatm$bmi<30,]

caddatfm<-caddat2[caddat2$sex==2,]
caddatfm<-caddatfm[caddatfm$smoke==0,]
caddatfm<-caddatfm[caddatfm$Diabetes==0,]
caddatfm<-caddatfm[caddatfm$bmi<30,]

caddatms<-caddatm[caddatm$smoke==1,]
#caddatmns<-caddatm[caddatm$smoke==0,]
caddatms<-caddatms[caddatms$Diabetes==0,]
caddatms<-caddatms[caddatms$bmi<30,]

caddatmns<-caddatms[caddatms$bmi>=30,]
caddatmns<-caddatmns[caddatmns$Diabetes==1,]

caddatfmnsdia<-caddat2[caddat2$sex==2,]
caddatfmnsdia<-caddatfmnsdia[caddatfmnsdia$smoke==0,]
caddatfmnsdia<-caddatfmnsdia[caddatfmnsdia$Diabetes==1,]
caddatfmnsdia

d0<-decilep(caddat2)
dm<-decilep(caddatm)
dfm<-decilep(caddatfm)
dms<-decilep(caddatms)
dmns<-decilep(caddatmns)
dfmnsdia<-decilep(caddatfmnsdia)

q0<-decilepoint(caddat2)
qm<-decilepoint(caddatm)
qfm<-decilepoint(caddatfm)
qms<-decilepoint(caddatms)
qmns<-decilepoint(caddatmns)
qfmnsdia<-decilepoint(caddatfmnsdia)


#simulation
ggplot()+ stat_smooth(data=d0,aes(x=numAllele,y=p),method = "glm", method.args = list(family = binomial(link = 'log')), se = FALSE, fullrange = T,color="grey",size=2)       + stat_smooth(data=dm,aes(x=numAllele,y=p),method = "glm", method.args = list(family = binomial(link = 'log')), se = FALSE, fullrange = T,color="#d3e0ea",size=2)+ stat_smooth(data=dfm,aes(x=numAllele,y=p),method = "glm", method.args = list(family = binomial(link = 'log')), se = FALSE, fullrange = T,color="#f5c0c0",size=2)+ stat_smooth(data=dms,aes(x=numAllele,y=p),method = "glm", method.args = list(family = binomial(link = 'log')), se = FALSE, fullrange = T,color="#d8ebe4",size=2)+ stat_smooth(data=dmns,aes(x=numAllele,y=p),method = "glm", method.args = list(family = binomial(link = 'log')), se = FALSE, fullrange = T,color="#b8b5ff",size=2)+ stat_smooth(data=dmns,aes(x=numAllele,y=p),method = "glm", method.args = list(family = binomial(link = 'log')), se = FALSE,color="#7868e6",size=2)   +stat_smooth(data=d0,aes(x=numAllele,y=p),method = "glm", method.args = list(family = binomial(link = 'log')), se = FALSE, color="black",size=2)                        + stat_smooth(data=dm,aes(x=numAllele,y=p),method = "glm", method.args = list(family = binomial(link = 'log')), se = FALSE,color="#1687a7",size=2)                       + stat_smooth(data=dfm,aes(x=numAllele,y=p),method = "glm", method.args = list(family = binomial(link = 'log')), se = FALSE, color="#f14668",size=2)              + stat_smooth(data=dms,aes(x=numAllele,y=p),method = "glm", method.args = list(family = binomial(link = 'log')), se = FALSE, color="#3a6351",size=2)      +theme_classic()+ylab("Prevalence (%)")+xlab("Mean number of risk alleles")+xlim(100,300)+labs(title="CAD")+ scale_y_continuous(labels = scales::percent)+ylim(0,1)

# dedcile 
ggplot() +geom_point(data=q0,aes(x=numAllele,y=p), color="black")+ stat_smooth(data=q0,aes(x=numAllele,y=p),method = "glm", method.args = list(family = binomial(link = 'log')), se = FALSE, color="black",size=2)+geom_errorbar(data=q0,aes(x=numAllele,ymin=p-p_se*10, ymax=p+p_se*10,fill="black"), size=0.3,width=1)   +geom_point(data=qms,aes(x=numAllele,y=p), color="#32a02c")+ stat_smooth(data=qms,aes(x=numAllele,y=p),method = "glm", method.args = list(family = binomial(link = 'log')), se = FALSE, color="#32a02c",size=2)+geom_errorbar(data=qms,aes(x=numAllele,ymin=p-p_se*5, ymax=p+p_se*5,fill="#32a02c"), size=0.3,width=1)                 +geom_point(data=qm,aes(x=numAllele,y=p), color="#ef9b2b")+ stat_smooth(data=qm,aes(x=numAllele,y=p),method = "glm", method.args = list(family = binomial(link = 'log')), se = FALSE, color="#ef9b2b",size=2)+geom_errorbar(data=qm,aes(x=numAllele,ymin=p-p_se*5, ymax=p+p_se*5,fill="#ef9b2b"), size=0.3,width=1)                          +geom_point(data=qfm,aes(x=numAllele,y=p), color="#f14668")+ stat_smooth(data=qfm,aes(x=numAllele,y=p),method = "glm", method.args = list(family = binomial(link = 'log')), se = FALSE, color="#f14668",size=2)+geom_errorbar(data=qfm,aes(x=numAllele,ymin=p-p_se*5, ymax=p+p_se*5,fill="#f14668"), size=0.3,width=1)                 +geom_point(data=qmns,aes(x=numAllele,y=p), color="#367eb8")+ stat_smooth(data=qmns,aes(x=numAllele,y=p),method = "glm", method.args = list(family = binomial(link = 'log')), se = FALSE, color="#367eb8",size=2)+geom_errorbar(data=qmns,aes(x=numAllele,ymin=p-p_se, ymax=p+p_se,fill="#367eb8"), size=0.3,width=1)+ stat_smooth(data=qms,aes(x=numAllele,y=p),method = "glm", method.args = list(family = binomial(link = 'log')), se = FALSE, fullrange = T,color="#d8ebe4",size=2)                 + scale_y_continuous(labels = scales::percent)+theme_classic()+ylab("Prevalence (%)")+xlab("Mean number of risk alleles")+labs(title="CAD")


# orginal data
ggplot()+geom_point(data=dmns,aes(x=numAllele,y=p),color="#367eb8")+stat_smooth(data=dmns,aes(x=numAllele,y=p),method = "glm", method.args = list(family = binomial(link = 'log')), se = FALSE,color="#367eb8",size=2)+geom_errorbar(data=dmns,aes(x=numAllele,ymin=p-p_se, ymax=p+p_se,fill="#367eb8"), size=0.3,width=1)       +geom_point(data=d0,aes(x=numAllele,y=p), color="black")+  stat_smooth(data=d0,aes(x=numAllele,y=p),method = "glm", method.args = list(family = binomial(link = 'log')), se = FALSE, color="black",size=2)+geom_errorbar(data=d0,aes(x=numAllele,ymin=p-p_se, ymax=p+p_se,fill="black"), size=0.3,width=1)                       +geom_point(data=dm,aes(x=numAllele,y=p),color="#ef9b2b")+ stat_smooth(data=dm,aes(x=numAllele,y=p),method = "glm", method.args = list(family = binomial(link = 'log')), se = FALSE,color="#ef9b2b",size=2)+geom_errorbar(data=dm,aes(x=numAllele,ymin=p-p_se, ymax=p+p_se,fill="#ef9b2b"), size=0.3,width=1)                       +geom_point(data=dfm,aes(x=numAllele,y=p), color="#f14668")+stat_smooth(data=dfm,aes(x=numAllele,y=p),method = "glm", method.args = list(family = binomial(link = 'log')), se = FALSE, color="#f14668",size=2)+geom_errorbar(data=dfm,aes(x=numAllele,ymin=p-p_se, ymax=p+p_se,fill="#f14668"), size=0.3,width=1)                        +geom_point(data=dms,aes(x=numAllele,y=p),color="#32a02c")+stat_smooth(data=dms,aes(x=numAllele,y=p),method = "glm", method.args = list(family = binomial(link = 'log')), se = FALSE, color="#32a02c",size=2)+geom_errorbar(data=dms,aes(x=numAllele,ymin=p-p_se, ymax=p+p_se,fill="#32a02c"), size=0.3,width=1)                       +theme_classic()+ylab("Prevalence (%)")+xlab("Mean number of risk alleles")+labs(title="CAD")+ scale_y_continuous(labels = scales::percent)+ylim(0,1)

#combined decile
ggplot()+geom_point(data=qmns,aes(x=numAllele,y=p),color="#367eb8")+stat_smooth(data=dmns,aes(x=numAllele,y=p),method = "glm", method.args = list(family = binomial(link = 'log')), se = FALSE,color="#367eb8",size=2,fullrange = T)+geom_errorbar(data=qmns,aes(x=numAllele,ymin=p-p_se, ymax=p+p_se,fill="#367eb8"), size=0.3,width=1)       +geom_point(data=q0,aes(x=numAllele,y=p), color="black")+  stat_smooth(data=d0,aes(x=numAllele,y=p),method = "glm", method.args = list(family = binomial(link = 'log')), se = FALSE, color="black",size=2)+geom_errorbar(data=q0,aes(x=numAllele,ymin=p-p_se*10, ymax=p+p_se*10,fill="black"), size=0.3,width=1)                       +geom_point(data=qm,aes(x=numAllele,y=p),color="#ef9b2b")+ stat_smooth(data=dm,aes(x=numAllele,y=p),method = "glm", method.args = list(family = binomial(link = 'log')), se = FALSE,color="#ef9b2b",size=2)+geom_errorbar(data=qm,aes(x=numAllele,ymin=p-p_se*5, ymax=p+p_se*5,fill="#ef9b2b"), size=0.3,width=1)                       +geom_point(data=qfm,aes(x=numAllele,y=p), color="#f14668")+stat_smooth(data=dfm,aes(x=numAllele,y=p),method = "glm", method.args = list(family = binomial(link = 'log')), se = FALSE, color="#f14668",size=2)+geom_errorbar(data=qfm,aes(x=numAllele,ymin=p-p_se*5, ymax=p+p_se*5,fill="#f14668"), size=0.3,width=1)                        +geom_point(data=qms,aes(x=numAllele,y=p),color="#32a02c")+stat_smooth(data=dms,aes(x=numAllele,y=p),method = "glm", method.args = list(family = binomial(link = 'log')), se = FALSE, color="#32a02c",size=2)+geom_errorbar(data=qms,aes(x=numAllele,ymin=p-p_se*5, ymax=p+p_se*5,fill="#32a02c"), size=0.3,width=1)                        +theme_classic()+ylab("Prevalence (%)")+xlab("Mean number of risk alleles")+labs(title="CAD")+ scale_y_continuous(labels = scales::percent)+ylim(0,0.5)
# dedcile 
p1<-ggplot() +geom_point(data=q0,aes(x=numAllele,y=p), color="black")+ stat_smooth(data=q0,aes(x=numAllele,y=p),method = "glm", method.args = list(family = binomial(link = 'log')), se = FALSE, color="black",size=2)+geom_errorbar(data=q0,aes(x=numAllele,ymin=p-p_se*10, ymax=p+p_se*10,fill="black"), size=0.3,width=1)   
p2<-p1+geom_point(data=qms,aes(x=numAllele,y=p), color="#32a02c")+ stat_smooth(data=qms,aes(x=numAllele,y=p),method = "glm", method.args = list(family = binomial(link = 'log')), se = FALSE, color="#32a02c",size=2)+geom_errorbar(data=qms,aes(x=numAllele,ymin=p-p_se*5, ymax=p+p_se*5,fill="#32a02c"), size=0.3,width=1)  
p3<-p2+geom_point(data=qm,aes(x=numAllele,y=p), color="#ef9b2b")+ stat_smooth(data=qm,aes(x=numAllele,y=p),method = "glm", method.args = list(family = binomial(link = 'log')), se = FALSE, color="#ef9b2b",size=2)+geom_errorbar(data=qm,aes(x=numAllele,ymin=p-p_se*5, ymax=p+p_se*5,fill="#ef9b2b"), size=0.3,width=1)
p4<-p3+geom_point(data=qfm,aes(x=numAllele,y=p), color="#f14668")+ stat_smooth(data=qfm,aes(x=numAllele,y=p),method = "glm", method.args = list(family = binomial(link = 'log')), se = FALSE, color="#f14668",size=2)+geom_errorbar(data=qfm,aes(x=numAllele,ymin=p-p_se*5, ymax=p+p_se*5,fill="#f14668"), size=0.3,width=1) 
p5<-p4+geom_point(data=qmns,aes(x=numAllele,y=p), color="#367eb8")+ stat_smooth(data=qmns,aes(x=numAllele,y=p),method = "glm", method.args = list(family = binomial(link = 'log')), se = FALSE, color="#367eb8",size=2)+geom_errorbar(data=qmns,aes(x=numAllele,ymin=p-p_se, ymax=p+p_se,fill="#367eb8"), size=0.3,width=1)
p6<-p5+stat_smooth(data=q0,aes(x=numAllele,y=p),method = "glm", method.args = list(family = binomial(link = 'log')), se = FALSE, color="grey",size=2,fullrange = T)+ scale_y_continuous(labels = scales::percent)+theme_classic()+ylab("Prevalence (%)")+xlab("Mean number of risk alleles")+labs(title="CAD")+ylim(0,0.5)+xlim(130,270)
#size
cairo_pdf("d7.pdf",width=13.515748, height=5.708661,pointsize = 12,family="Univers")
dev.off()

# figure 4
#new version difference, quadriac model C-statistics
plot_diff_prev1<-function(dn,trait,pheno1,pheno2,figure,quant1,quant2){
  #parameter:
  #d-"n","N","p","numAllele","p_se","px","px_se","group","colors","p2" 
  #trait:Disease; pheno1:E abs;pheno2: E pre;figure:figure name; number: min(d$p2);seq: range (d$p2);quant1:0.1 quantile;quant2:0.9 quantile
  #d$p<-d$p2-min(d$p2)
  # create barplot
  #bar1<-as.data.frame(t(d[d$group==pheno1,3]))
  #bar2<-as.data.frame(t(d[d$group==pheno2,3]))
  #bar<-rbind(bar1,bar2)
  #rownames(bar)<-c(pheno1,pheno2)
  #colnames(bar)<-seq(1,10,1)
  dn$p_se   <- sqrt(dn$p*(1-dn$p)/dn$N)
  dn$px     <- dn$p / dn$numAllele
  dn$px_se  <- dn$p_se / dn$numAllele
  d1<-dn[dn$group==pheno1,]
  d2<-dn[dn$group==pheno2,]
  
  d1$colors<-"#33A2C9"
  d2$colors<-"#CD375D"
  
  bar<-rbind(mean(d1$p),mean(d2$p))
  
  bar<-as.matrix(bar)
  
  d<-rbind(d1,d2)
  
  l    <- which(d[,"p_se"]>0)
  #modelLogit1  <- glm(cbind(n,N)~numAllele,data=d1,family = binomial(link="logit"))
  modelLogit1  <- glm(p~numAllele,data=d1,family = binomial(link="logit"),weights=N)
  modelProbit1 <- glm(cbind(n,N)~numAllele,data=d1,family = binomial(link="probit"))
  modelLog1    <- glm(cbind(n,N)~numAllele,data=d1,family = binomial(link="log"))
  modelLinear1 <- lm(p~numAllele,data=d1,weights = N)
  
  predLogit1   <- predict(modelLogit1,type = "response")
  predProbit1  <- predict(modelProbit1,type = "response")
  predLog1     <- predict(modelLog1,type = "response")
  predLin1     <- predict(modelLinear1,type="response")
  
  Rlogit1      <- round(unlist(cor.test(d1$p,predLogit1)[c("estimate","conf.int")]),2)
  Rprobit1     <- round(unlist(cor.test(d1$p,predProbit1)[c("estimate","conf.int")]),2)
  Rlog1        <- round(unlist(cor.test(d1$p,predLog1)[c("estimate","conf.int")]),2)
  Rlin1        <- round(unlist(cor.test(d1$p,predLin1)[c("estimate","conf.int")]),2)
  
  Cols1 <- c("white","#33A2C9","white","white","white","white")
  
  #modelLogit2  <- glm(cbind(n,N)~numAllele,data=d2,family = binomial(link="logit"))
  modelLogit2  <- glm(p~numAllele,data=d2,family = binomial(link="logit"),weights=N)
  modelProbit2 <- glm(cbind(n,N)~numAllele,data=d2,family = binomial(link="probit"))
  modelLog2    <- glm(cbind(n,N)~numAllele,data=d2,family = binomial(link="log"))
  modelLinear2 <- lm(p~numAllele,data=d2,weights = N)
  
  predLogit2   <- predict(modelLogit2,type = "response")
  predProbit2  <- predict(modelProbit2,type = "response")
  predLog2     <- predict(modelLog2,type = "response")
  predLin2     <- predict(modelLinear2,type="response")
  
  Rlogit2      <- round(unlist(cor.test(d2$p,predLogit2)[c("estimate","conf.int")]),2)
  Rprobit2     <- round(unlist(cor.test(d2$p,predProbit2)[c("estimate","conf.int")]),2)
  Rlog2        <- round(unlist(cor.test(d2$p,predLog2)[c("estimate","conf.int")]),2)
  Rlin2        <- round(unlist(cor.test(d2$p,predLin2)[c("estimate","conf.int")]),2)
  
  Cols2 <- c("white","#CD375D","white","white","white","white")
  
  d1$co<-d1$N-d1$n
  d2$co<-d2$N-d2$n
  obs1<-d1$co+d1$n
  obs2<-d2$co+d2$n  
  d10<-d1[rep(1:length(obs1),each=2),]
  d10$survived<-rep(c(0L,1L),times=length(obs1))
  d10$weight<-c(rbind(d1$n,d1$co))
  d10<-d10[rep(1:nrow(d10),times=d10$weight),]
  model<-glm(survived~numAllele,data=d10,family=binomial(link="logit"))
  ci1<-auc(model$fitted.values,d10$survived)
  
  d20<-d2[rep(1:length(obs1),each=2),]
  d20$survived<-rep(c(0L,1L),times=length(obs1))
  d20$weight<-c(rbind(d2$n,d2$co))
  d20<-d20[rep(1:nrow(d20),times=d20$weight),]
  model<-glm(survived~numAllele,data=d20,family=binomial(link="logit"))
  ci2<-auc(model$fitted.values,d20$survived)
  
  #pdf("CAD_smoke_uGRS_ori_dif.pdf",width=4.63,height=5.64)
  cairo_pdf(paste(figure,"_panelAB.pdf",sep=""),width=8.80,height=7.00)
  par(mai=c(0.5,0.5,0.5,0.5))
  layout(matrix(c(1,2,0,0), 1, 2, byrow = TRUE),widths=c(1,3))
  #bar plot
  par(mai=c(2,1,0.5,0))
  barplot(bar,beside=TRUE,col=c("#33A2C9","#CD375D"),ylab=paste("Mean ",trait," prevalence",sep=""),border = NA)
  #barplot(bar,beside=TRUE,col=c("#33A2C9","#CD375D","#EEB422"),ylab=paste("Mean ",trait," prevalence",sep=""),border = NA)
  legend(x=0,y=xy[3]-yinch(0.7),xpd=TRUE,legend=c(paste0(pheno1),paste0(pheno2)),col=c("#33A2C9","#CD375D"),cex=0.7,box.lty=0,lty=c(1,1),lwd=10,bg="transparent")
  #legend(x=0,y=xy[3]-yinch(0.7),xpd=TRUE,legend=c(paste0(pheno1),paste0(pheno2),paste0(pheno3)),col=c("#33A2C9","#CD375D","#EEB422"),cex=0.7,box.lty=0,lty=c(1,1),lwd=10,bg="transparent")
  # line plot
  par(mai=c(2,0.5,0.5,1))
  plot(d[l,"numAllele"],d[l,"p"],pch=19,axes=FALSE,col=d$colors,ylim=c(0,max(d$p+d$p_se)),xlim=c(min(d$numAllele-sd(d$numAllele)),max(d$numAllele+sd(d$numAllele))),xlab="",ylab="");axis(1);axis(4);mtext(paste("Prevalence differences in each decile of \n",pheno1," and ",pheno2," in ",trait,sep=""),side=4,line=3);mtext("Mean number of risk alleles",side=1,line=2.5)
  
  for(i in 1:nrow(d1)){
    if(d1[i,"p_se"]>0){
      arrows(d1[i,"numAllele"],d1[i,"p"]-d1[i,"p_se"],d1[i,"numAllele"],d1[i,"p"]+d1[i,"p_se"],col="#A5DFEF",angle=90,len=0.03)
      arrows(d1[i,"numAllele"],d1[i,"p"]+d1[i,"p_se"],d1[i,"numAllele"],d1[i,"p"]-d1[i,"p_se"],col="#A5DFEF",angle=90,len=0.03)
    }
  }
  
  for(i in 1:nrow(d2)){
    if(d2[i,"p_se"]>0){
      arrows(d2[i,"numAllele"],d2[i,"p"]-d2[i,"p_se"],d2[i,"numAllele"],d2[i,"p"]+d2[i,"p_se"],col="#EF9AB2",angle=90,len=0.04)
      arrows(d2[i,"numAllele"],d2[i,"p"]+d2[i,"p_se"],d2[i,"numAllele"],d2[i,"p"]-d2[i,"p_se"],col="#EF9AB2",angle=90,len=0.04)
    }
  }
  
  matlines(d1$numAllele,predLogit1,col="#33A2C9",lwd=2)
  matlines(d2$numAllele,predLogit2,col="#CD375D",lwd=2)
  #CAD uGRS 181.685,203.873
  polygon(x=c(quant1,quant1,quant2,quant2),y=c(0,max(d$p),max(d$p),0),col=rgb(0.2,0.2,0.2,0.2),border = NA)
  
  xy<-par("usr")
  #legend(x=162,y=xy[3]-yinch(0.7),xpd=TRUE,
  legend(x=mean(d$numAllele)-sd(d$numAllele),y=xy[3]-yinch(0.7),xpd=TRUE,
         legend=c(
           paste0(pheno1),
           #paste0("Logit: R=",Rlogit1[1]," (95%CI: [",Rlogit1[2],"-",Rlogit1[3],"])","; Steepness: ",round(modelLogit1$coefficients[2],digits=2)),
           paste0("Logit: R=",Rlogit1[1]," (95%CI: [",Rlogit1[2],"-",Rlogit1[3],"])"),
           paste0("Log: R=",Rlog1[1]," (95%CI: [",Rlog1[2],"-",Rlog1[3],"])"),
           paste0("Probit: R=",Rprobit1[1]," (95%CI: [",Rprobit1[2],"-",Rprobit1[3],"])"),
           paste0("Linear: R=",Rlin1[1]," (95%CI: [",Rlin1[2],"-",Rlin1[3],"])"),
           paste0("C-index: ",round(ci1[1],digits=3),"; Steepness: ",round(modelLogit1$coefficients[2],digits=3))),
         col=c(Cols1),cex=0.7,box.lty=0,lty=c(1,1,1,3),lwd=2,bg="transparent")
  
  #legend(x=197,y=xy[3]-yinch(0.7),xpd=TRUE,
  legend(x=mean(d$numAllele)+sd(d$numAllele),y=xy[3]-yinch(0.7),xpd=TRUE,
         legend=c(
           #paste0("Ever_smoked"),
           #paste0("BMI>=30"),
           paste0(pheno2),
           paste0("Logit: R=",Rlogit2[1]," (95%CI: [",Rlogit2[2],"-",Rlogit2[3],"])"),
           paste0("Log: R=",Rlog2[1]," (95%CI: [",Rlog2[2],"-",Rlog2[3],"])"),
           paste0("Probit: R=",Rprobit2[1]," (95%CI: [",Rprobit2[2],"-",Rprobit2[3],"])"),
           paste0("Linear: R=",Rlin2[1]," (95%CI: [",Rlin2[2],"-",Rlin2[3],"])"),
           paste0("C-index: ",round(ci2[1],digits=3),"; Steepness: ",round(modelLogit2$coefficients[2],digits=3))),
         col=c(Cols2),cex=0.7,box.lty=0,lty=c(1,1,1,3),lwd=2,bg="transparent")
  
  
  dev.off()
  
  d<-d1
  d$p<-d2$p-d1$p
  d$N<-d1$N+d2$N
  d$p_se<-sqrt(d1$p*(1-d1$p)/d1$N + d2$p*(1-d2$p)/d2$N)
  d$numAllele2<-d$numAllele^2
  
  l    <- which(d[,"p_se"]>0)
  
  modelLinear <- lm(p~1+numAllele,data=d,weights = 1/p_se)
  modelLinear2 <- lm(p~1+numAllele+numAllele2,data=d,weights = 1/p_se)
  
  predLin     <- predict(modelLinear,type="response")
  predLin2     <- predict(modelLinear2,type="response")
  
  Rlin        <- round(unlist(cor.test(d$p,predLin)[c("estimate","conf.int")]),2)
  Rlin2        <- round(unlist(cor.test(d$p,predLin2)[c("estimate","conf.int")]),2)
  
  Cols <- c("seagreen","coral")
  
  cairo_pdf(paste(figure,"_panelc.pdf",sep=""),width=5.746,height=7)
  plot(d[l,"numAllele"],d[l,"p"],pch=19,axes=FALSE,ylim=c(min(d$p-d$p_se),max(d$p+d$p_se)),xlim=c(min(d$numAllele-sd(d$numAllele)),max(d$numAllele+sd(d$numAllele))),xlab="Mean number of risk alleles",ylab="");axis(1);axis(2);mtext(paste("Prevalence differences in each decile between\n",pheno1," and ",pheno2, " groups in ",trait," (%)",sep=""),side=2,line=2)
  
  for(i in 1:nrow(d)){
    if(d[i,"p_se"]>0){
      arrows(d[i,"numAllele"],d[i,"p"]-d[i,"p_se"],d[i,"numAllele"],d[i,"p"]+d[i,"p_se"],col="grey",angle=90,len=0.05)
      arrows(d[i,"numAllele"],d[i,"p"]+d[i,"p_se"],d[i,"numAllele"],d[i,"p"]-d[i,"p_se"],col="grey",angle=90,len=0.05)
    }
  }
  
  matlines(d$numAllele,cbind(predLin,predLin),col="seagreen",lty=3)
  matlines(d$numAllele,cbind(predLin2,predLin2),col="coral1")
  
  polygon(x=c(quant1,quant1,quant2,quant2),y=c(0,max(d$p),max(d$p),0),col=rgb(0.2,0.2,0.2,0.2),border = NA)
  
  legend("topleft",
         legend=c(paste0("Linear: R=",Rlin[1]," (95%CI: [",Rlin[2],"-",Rlin[3],"])"),
                  paste0("Quadratic: R=",Rlin2[1]," (95%CI: [",Rlin2[2],"-",Rlin2[3],"])")),
         col=Cols,cex=0.9,box.lty=0,lty=c(3,1),lwd=2,bg="transparent")
  
  dev.off()
}


plot_diff_prev2<-function(dn,trait,pheno1,pheno2,figure,number,seq,quant1,quant2){
  #parameter:
  #d-"n","N","p","numAllele","p_se","px","px_se","group","colors","p2" 
  #trait:Disease; pheno1:E abs;pheno2: E pre;figure:figure name; number: min(d$p2);seq: range (d$p2);quant1:0.1 quantile;quant2:0.9 quantile
  #d$p<-d$p2-min(d$p2)
  # create barplot
  #bar1<-as.data.frame(t(d[d$group==pheno1,3]))
  #bar2<-as.data.frame(t(d[d$group==pheno2,3]))
  #bar<-rbind(bar1,bar2)
  #rownames(bar)<-c(pheno1,pheno2)
  #colnames(bar)<-seq(1,10,1)
  dn$p_se   <- sqrt(dn$p*(1-dn$p)/dn$N)
  dn$px     <- dn$p / dn$numAllele
  dn$px_se  <- dn$p_se / dn$numAllele
  d1<-dn[dn$group==pheno1,]
  d2<-dn[dn$group==pheno2,]
  d1$colors<-"#33A2C9"
  d2$colors<-"#CD375D"
  bar<-rbind(mean(d1$p),mean(d2$p))
  bar<-as.matrix(bar)
  d1$p2<-d1$p-mean(d1$p)
  d2$p2<-d2$p-mean(d2$p)
  d<-rbind(d1,d2)
  d$p<-d$p2+number
  d$N<-14500
  d$n<-ceiling(d$p*14500)
  d$p_se   <- sqrt(d$p*(1-d$p)/d$N)*5
  d$px     <- d$p / d$numAllele
  d$px_se  <- d$p_se / d$numAllele
  d$numAllele2<-d$numAllele^2
  d1<-d[d$group==pheno1,]
  d2<-d[d$group==pheno2,]
  
  
  
  l    <- which(d[,"p_se"]>0)
  modelLogit1  <- glm(cbind(n,N)~numAllele,data=d1,family = binomial(link="logit"))
  modelLinear11 <- lm(p~numAllele,data=d1,weights = N)
  modelLinear12 <- lm(p~1+numAllele+numAllele2,data=d1,weights = 1/p_se)
  
  predLogit1   <- predict(modelLogit1,type = "response")  
  predLin11     <- predict(modelLinear11,type="response")
  predLin12     <- predict(modelLinear12,type="response")
  
  Rlogit1      <- round(unlist(cor.test(d1$p,predLogit1)[c("estimate","conf.int")]),2)
  Rlin11        <- round(unlist(cor.test(d1$p,predLin11)[c("estimate","conf.int")]),2)
  Rlin12        <- round(unlist(cor.test(d1$p,predLin12)[c("estimate","conf.int")]),2)
  
  Cols1 <- c("white","white","#33A2C9","white")
  
  modelLogit2  <- glm(cbind(n,N)~numAllele,data=d2,family = binomial(link="logit"))
  modelLinear21 <- lm(p~numAllele,data=d2,weights = N)
  modelLinear22 <- lm(p~1+numAllele+numAllele2,data=d2,weights = 1/p_se)
  
  predLogit2   <- predict(modelLogit2,type = "response") 
  predLin21     <- predict(modelLinear21,type="response")
  predLin22     <- predict(modelLinear22,type="response")
  
  Rlogit2      <- round(unlist(cor.test(d2$p,predLogit2)[c("estimate","conf.int")]),2)
  Rlin21        <- round(unlist(cor.test(d2$p,predLin21)[c("estimate","conf.int")]),2)
  Rlin22        <- round(unlist(cor.test(d2$p,predLin22)[c("estimate","conf.int")]),2)
  
  Cols2 <- c("white","white","#CD375D","white")
  
  #pdf("CAD_smoke_uGRS_ori_dif.pdf",width=4.63,height=5.64)
  cairo_pdf(paste(figure,"_Q_panelAB.pdf",sep=""),width=8.80,height=7.00)
  par(mai=c(0.5,0.5,0.5,0.5))
  layout(matrix(c(1,2,0,0), 1, 2, byrow = TRUE),widths=c(1,3))
  #bar plot
  par(mai=c(2,1,0.5,0))
  barplot(bar,beside=TRUE,col=c("#33A2C9","#CD375D"),ylab=paste("Mean ",trait," prevalence",sep=""),border = NA)
  #barplot(bar,beside=TRUE,col=c("#33A2C9","#CD375D","#EEB422"),ylab=paste("Mean ",trait," prevalence",sep=""),border = NA)
  legend(x=0,y=xy[3]-yinch(0.7),xpd=TRUE,legend=c(paste0(pheno1),paste0(pheno2)),col=c("#33A2C9","#CD375D"),cex=0.7,box.lty=0,lty=c(1,1),lwd=10,bg="transparent")
  #legend(x=0,y=xy[3]-yinch(0.7),xpd=TRUE,legend=c(paste0(pheno1),paste0(pheno2),paste0(pheno3)),col=c("#33A2C9","#CD375D","#EEB422"),cex=0.7,box.lty=0,lty=c(1,1),lwd=10,bg="transparent")
  # line plot
  par(mai=c(2,0.5,0.5,1))
  plot(d[l,"numAllele"],d[l,"p"],pch=19,axes=FALSE,col=d$colors,ylim=c(min(d$p-d$p_se),max(d$p+d$p_se)),xlim=c(min(d$numAllele-sd(d$numAllele)),max(d$numAllele+sd(d$numAllele))),xlab="",ylab="");axis(1);axis(4,at=(seq+number),label=seq,cex.axis=0.8);mtext(paste("Prevalence differences in each decile compared\n to the average prevalence in ",trait,sep=""),side=4,line=3);mtext("Mean number of risk alleles",side=1,line=2.5)
  
  for(i in 1:nrow(d1)){
    if(d1[i,"p_se"]>0){
      arrows(d1[i,"numAllele"],d1[i,"p"]-d1[i,"p_se"],d1[i,"numAllele"],d1[i,"p"]+d1[i,"p_se"],col="#A5DFEF",angle=90,len=0.03)
      arrows(d1[i,"numAllele"],d1[i,"p"]+d1[i,"p_se"],d1[i,"numAllele"],d1[i,"p"]-d1[i,"p_se"],col="#A5DFEF",angle=90,len=0.03)
    }
  }
  
  for(i in 1:nrow(d2)){
    if(d2[i,"p_se"]>0){
      arrows(d2[i,"numAllele"],d2[i,"p"]-d2[i,"p_se"],d2[i,"numAllele"],d2[i,"p"]+d2[i,"p_se"],col="#EF9AB2",angle=90,len=0.04)
      arrows(d2[i,"numAllele"],d2[i,"p"]+d2[i,"p_se"],d2[i,"numAllele"],d2[i,"p"]-d2[i,"p_se"],col="#EF9AB2",angle=90,len=0.04)
    }
  }
  
  matlines(d1$numAllele,predLin12,col="#33A2C9",lwd=2)
  matlines(d2$numAllele,predLin22,col="#CD375D",lwd=2)
  #CAD uGRS 181.685,203.873
  polygon(x=c(quant1,quant1,quant2,quant2),y=c(0,max(d$p),max(d$p),0),col=rgb(0.2,0.2,0.2,0.2),border = NA)
  
  xy<-par("usr")
  #legend(x=162,y=xy[3]-yinch(0.7),xpd=TRUE,
  legend(x=mean(d$numAllele)-sd(d$numAllele),y=xy[3]-yinch(0.7),xpd=TRUE,
         legend=c(
           paste0(pheno1),
           #paste0("Logit: R=",Rlogit1[1]," (95%CI: [",Rlogit1[2],"-",Rlogit1[3],"])","; Steepness: ",round(modelLogit1$coefficients[2],digits=2)),
           paste0("Linear: R=",Rlin11[1]," (95%CI: [",Rlin11[2],"-",Rlin11[3],"])"),
           paste0("Quadratic: R=",Rlin12[1]," (95%CI: [",Rlin12[2],"-",Rlin12[3],"])"),
           paste0("Logit: R=",Rlogit1[1]," (95%CI: [",Rlogit1[2],"-",Rlogit1[3],"])","; Steepness: ",round(modelLogit1$coefficients[2],digits=2))),
         col=c(Cols1),cex=0.7,box.lty=0,lty=c(1,3),lwd=2,bg="transparent")
  
  #legend(x=197,y=xy[3]-yinch(0.7),xpd=TRUE,
  legend(x=mean(d$numAllele)+sd(d$numAllele),y=xy[3]-yinch(0.7),xpd=TRUE,
         legend=c(
           #paste0("Ever_smoked"),
           #paste0("BMI>=30"),
           paste0(pheno2),
           paste0("Linear: R=",Rlin21[1]," (95%CI: [",Rlin21[2],"-",Rlin21[3],"])"),
           paste0("Quadratic: R=",Rlin22[1]," (95%CI: [",Rlin22[2],"-",Rlin22[3],"])"),
           paste0("Logit: R=",Rlogit1[1]," (95%CI: [",Rlogit1[2],"-",Rlogit1[3],"])","; Steepness: ",round(modelLogit2$coefficients[2],digits=2))),
         col=c(Cols2),cex=0.7,box.lty=0,lty=c(1,3),lwd=2,bg="transparent")
  
  
  dev.off()
}

plot_diff_prev3<-function(dn,trait,pheno1,pheno2,figure,quant1,quant2){
  #parameter:
  #d-"n","N","p","numAllele","p_se","px","px_se","group","colors","p2" 
  #trait:Disease; pheno1:E abs;pheno2: E pre;figure:figure name; number: min(d$p2);seq: range (d$p2);quant1:0.1 quantile;quant2:0.9 quantile
  #d$p<-d$p2-min(d$p2)
  # create barplot
  #bar1<-as.data.frame(t(d[d$group==pheno1,3]))
  #bar2<-as.data.frame(t(d[d$group==pheno2,3]))
  #bar<-rbind(bar1,bar2)
  #rownames(bar)<-c(pheno1,pheno2)
  #colnames(bar)<-seq(1,10,1)
  dn$numAllele2<-dn$numAllele^2
  d1<-dn[dn$group==pheno1,]
  d2<-dn[dn$group==pheno2,]
  d1$colors<-"#33A2C9"
  d2$colors<-"#CD375D"
  bar<-as.data.frame(rbind(mean(d1$p),mean(d2$p)))
  bar$pheno<-factor(c(pheno1,pheno2),levels=c(pheno1,pheno2))
  
  
  d1$p2<-d1$p-mean(d1$p)
  d2$p2<-d2$p-mean(d2$p)
  d<-rbind(d1,d2)
  d$p<-d$p2*100
  d$p_se<- (sqrt(d1$p*(1-d1$p)/d1$N + d2$p*(1-d2$p)/d2$N))*100
  d1<-d[d$group==pheno1,]
  d2<-d[d$group==pheno2,]
  
  l    <- which(d[,"p_se"]>0)
  modelLinear11 <- lm(p~numAllele,data=d1,weights = N)
  modelLinear12 <- lm(p~1+numAllele+numAllele2,data=d1,weights = 1/p_se)
  
  predLin11     <- predict(modelLinear11,type="response")
  predLin12     <- predict(modelLinear12,type="response")
  
  Rlin11        <- round(unlist(cor.test(d1$p,predLin11)[c("estimate","conf.int")]),2)
  Rlin12        <- round(unlist(cor.test(d1$p,predLin12)[c("estimate","conf.int")]),2)
  
  Cols1 <- c("white","white","#33A2C9","white")
  
  modelLinear21 <- lm(p~numAllele,data=d2,weights = N)
  modelLinear22 <- lm(p~1+numAllele+numAllele2,data=d2,weights = 1/p_se)
  
  predLin21     <- predict(modelLinear21,type="response")
  predLin22     <- predict(modelLinear22,type="response")
  
  Rlin21        <- round(unlist(cor.test(d2$p,predLin21)[c("estimate","conf.int")]),2)
  Rlin22        <- round(unlist(cor.test(d2$p,predLin22)[c("estimate","conf.int")]),2)
  
  Cols2 <- c("white","white","#CD375D","white")
  
  #pdf("CAD_smoke_uGRS_ori_dif.pdf",width=4.63,height=5.64)
  cairo_pdf(paste(figure,"_Q.pdf",sep=""),width=11.29,height=14.57)
  par(mai=c(0.5,0.5,0.5,0.5),cex.lab=1.5, cex.axis=2.5)
  layout(matrix(c(1,2,2,2,3,3,4,4), 2, 4, byrow = TRUE))
  #bar plot
  par(mai=c(0.5,1,1,0.5))
  barplot(V1*100~pheno,data=bar,beside=TRUE,space=0,col=c("#33A2C9","#CD375D"),ylim=c(0,25),ylab=paste("Mean ",trait," prevalence (%)",sep=""),border = NA,args.legend = list(x = "topright"))
  #smot
  smot<-matrix(c(pheno1,"Low risk",2.05,pheno2,"Low risk",8.62,pheno1,"Medium risk",3.99,pheno2,"Medium risk",15.73,pheno1,"High risk",7.49,pheno2,"High risk",25.39),ncol=3,byrow=TRUE)
  smot<-as.data.frame(smot)
  smot<-unfactor(smot)
  smot$V2 <- factor(smot$V2, levels = c("Low risk","Medium risk","High risk"))
  smot$V1 <- factor(smot$V1, levels = c(pheno1,pheno2))
  smot$V3<-as.numeric(smot$V3)
  par(mai=c(0.5,2,1,0.5))
  barplot(V3~V1+V2,data=smot,beside=TRUE,col=c("#33A2C9","#CD375D"),ylim=c(0,25),ylab=paste(trait," prevalence (%)",sep=""),border = NA,args.legend = list(x = "topright"), xlab="GRS categories",)
  #print(ggplot(data=smot)+geom_bar(aes(y=V3,fill=V2)))
  # line plot
  #par(mai=c(2,0.5,0.5,1))
  par(mai=c(1,1,1,0))
  plot(d[l,"numAllele"],d[l,"p"],pch=19,axes=FALSE,col=d$colors,ylim=c(min(d$p-d$p_se),max(d$p+d$p_se)),xlim=c(min(d$numAllele-sd(d$numAllele)),max(d$numAllele+sd(d$numAllele))),xlab="Mean number of risk alleles",ylab=paste("Prevalence differences (%) in each decile compared\n to the average prevalence in ",trait,sep=""));axis(1);axis(2,pos=170)
  
  for(i in 1:nrow(d1)){
    if(d1[i,"p_se"]>0){
      arrows(d1[i,"numAllele"],d1[i,"p"]-d1[i,"p_se"],d1[i,"numAllele"],d1[i,"p"]+d1[i,"p_se"],col="#A5DFEF",angle=90,len=0.03)
      arrows(d1[i,"numAllele"],d1[i,"p"]+d1[i,"p_se"],d1[i,"numAllele"],d1[i,"p"]-d1[i,"p_se"],col="#A5DFEF",angle=90,len=0.03)
    }
  }
  
  for(i in 1:nrow(d2)){
    if(d2[i,"p_se"]>0){
      arrows(d2[i,"numAllele"],d2[i,"p"]-d2[i,"p_se"],d2[i,"numAllele"],d2[i,"p"]+d2[i,"p_se"],col="#EF9AB2",angle=90,len=0.04)
      arrows(d2[i,"numAllele"],d2[i,"p"]+d2[i,"p_se"],d2[i,"numAllele"],d2[i,"p"]-d2[i,"p_se"],col="#EF9AB2",angle=90,len=0.04)
    }
  }
  
  matlines(d1$numAllele,predLin12,col="#33A2C9",lwd=2)
  matlines(d2$numAllele,predLin22,col="#CD375D",lwd=2)
  #CAD uGRS 181.685,203.873
  polygon(x=c(quant1,quant1,quant2,quant2),y=c(min(d$p),max(d$p),max(d$p),min(d$p)),col=rgb(0.2,0.2,0.2,0.2),border = NA)
  
  xy<-par("usr")
  #legend(x=162,y=xy[3]-yinch(0.7),xpd=TRUE,
  legend(x=mean(d$numAllele)-sd(d$numAllele),y=xy[3]-yinch(0.7),xpd=TRUE,
         legend=c(
           paste0(pheno1),
           paste0("Quadratic: R=",Rlin12[1]," (95%CI: [",Rlin12[2],"-",Rlin12[3],"])")),
         col=c(Cols1),cex=0.7,box.lty=0,lty=c(1,3),lwd=2,bg="transparent")
  
  #legend(x=197,y=xy[3]-yinch(0.7),xpd=TRUE,
  legend(x=mean(d$numAllele)+sd(d$numAllele),y=xy[3]-yinch(0.7),xpd=TRUE,
         legend=c(
           paste0(pheno2),
           paste0("Quadratic: R=",Rlin22[1]," (95%CI: [",Rlin22[2],"-",Rlin22[3],"])")),
         col=c(Cols2),cex=0.7,box.lty=0,lty=c(1,3),lwd=2,bg="transparent")
  
  d1<-dn[dn$group==pheno1,]
  d2<-dn[dn$group==pheno2,]
  d<-d1
  d$p<-(d2$p-d1$p)*100
  d$N<-d1$N+d2$N
  d$p_se<-(sqrt(d1$p*(1-d1$p)/d1$N + d2$p*(1-d2$p)/d2$N))*100
  d$numAllele2<-d$numAllele^2
  
  l    <- which(d[,"p_se"]>0)
  
  modelLinear2 <- lm(p~1+numAllele+numAllele2,data=d,weights = 1/p_se)
  predLin2     <- predict(modelLinear2,type="response")
  Rlin2        <- round(unlist(cor.test(d$p,predLin2)[c("estimate","conf.int")]),2)
  
  Cols <- c("coral")
  
  #cairo_pdf(paste(figure,"_panelc.pdf",sep=""),width=5.746,height=7)
  par(mai=c(1,0,1,0.5))
  plot(d[l,"numAllele"],d[l,"p"],pch=19,axes=FALSE,ylim=c(min(d$p-d$p_se),max(d$p+d$p_se)),xlim=c(min(d$numAllele-sd(d$numAllele)),max(d$numAllele+sd(d$numAllele))),xlab="Mean number of risk alleles",ylab="");axis(1);axis(4,pos=210);mtext(paste("Prevalence differences in each decile compared\n to the average prevalence in ",trait,sep=""),side=4,line=4,cex=1.5)
  
  for(i in 1:nrow(d)){
    if(d[i,"p_se"]>0){
      arrows(d[i,"numAllele"],d[i,"p"]-d[i,"p_se"],d[i,"numAllele"],d[i,"p"]+d[i,"p_se"],col="grey",angle=90,len=0.05)
      arrows(d[i,"numAllele"],d[i,"p"]+d[i,"p_se"],d[i,"numAllele"],d[i,"p"]-d[i,"p_se"],col="grey",angle=90,len=0.05)
    }
  }
  
  matlines(d$numAllele,cbind(predLin2,predLin2),col="coral1")
  polygon(x=c(quant1,quant1,quant2,quant2),y=c(0,max(d$p),max(d$p),0),col=rgb(0.2,0.2,0.2,0.2),border = NA)
  
  legend("topleft",
         legend=c(paste0("Quadratic: R=",Rlin2[1]," (95%CI: [",Rlin2[2],"-",Rlin2[3],"])")),
         col=Cols,cex=0.9,box.lty=0,lty=c(1),lwd=2,bg="transparent")
  dev.off()
}

## example ggplot2
library(ggplot2)
library(tidyverse)
library(data.table)
library(varhandle)
library(ggsignif)
fite<-function(x,y){
  bar2<-as.data.frame(c(x,y))
  colnames(bar2)<-"V1"
  bar2$V2<-100-bar$V1
  return(bar2)
}

plot_diff_prev6<-function(dn,trait,pheno1,pheno2,figure,quant1,quant2,number){
  dn$numAllele2<-dn$numAllele^2
  d1<-dn[dn$group==pheno1,]
  d2<-dn[dn$group==pheno2,]
  d1$colors<-"#33A2C9"
  d2$colors<-"#CD375D"
  dn$p_se<- sqrt(d1$p*(1-d1$p)/d1$N + d2$p*(1-d2$p)/d2$N)
  d1<-dn[dn$group==pheno1,]
  d2<-dn[dn$group==pheno2,]
  bar<-as.data.frame(rbind(mean(d1$p)*100,mean(d2$p)*100))
  bar$pheno<-factor(c(pheno1,pheno2),levels=c(pheno1,pheno2))
  smot<-matrix(c(pheno1,"Low risk",min(d1$p),pheno2,"Low risk",min(d2$p),pheno1,"Medium risk",mean(d1$p[2:9]),pheno2,"Medium risk",mean(d2$p[2:9]),pheno1,"High risk",max(d1$p),pheno2,"High risk",max(d2$p)),ncol=3,byrow=TRUE)
  
  modelLog11 <- glm(p~numAllele,data=d1,family = binomial(link="logit"),weights=N)
  predLog11     <- predict(modelLog11,type="response")
  Rlog11        <- round(unlist(cor.test(d1$p,predLog11)[c("estimate","conf.int")]),2)
  modelLog21 <- glm(p~numAllele,data=d2,family = binomial(link="logit"),weights=N)
  predLog21     <- predict(modelLog21,type="response")
  Rlog21        <- round(unlist(cor.test(d2$p,predLog21)[c("estimate","conf.int")]),2)
  d1$predLog<-predLog11
  d2$predLog<-predLog21
  dn<-rbind(d1,d2)
  pol<-data.frame(x=c(quant1,quant1,quant2,quant2),y=c(min(dn$p),max(dn$p),max(dn$p),min(dn$p)))
  cols <- c(pheno1="#33A2C9",pheno2="#CD375D","Difference in prevalence (%)"="#cccccc")
  names(cols)[1:2]<-c(pheno1,pheno2)
  p0<-ggplot()+geom_polygon(data=pol, mapping=aes(x=x, y=y),fill="grey",alpha=0.4)+geom_point(data=dn,aes(x=numAllele,y=p,color=group),pch=19)+geom_line(data=dn,aes(x=numAllele,y=predLog,group=group,color=group),size=1.5)+theme_classic()+ylab(paste(trait," prevalence (%)",sep=""))+geom_errorbar(data=dn,aes(x=numAllele,ymin=p-p_se, ymax=p+p_se,color=group), size=0.3,width=0.1)+xlab("wGRS")+ theme(legend.title =element_blank(),legend.position="top",legend.text=element_text(size=13),plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),text = element_text(size=18))+ scale_color_manual(values=cols[1:2],breaks = names(cols[1:2]),labels=c(paste0(pheno1,": Logit: R=",Rlog11[1]," (95%CI: [",Rlog11[2],"-",Rlog11[3],"])"),paste0(pheno2,": Logit: R=",Rlog21[1]," (95%CI: [",Rlog21[2],"-",Rlog21[3],"])")))+guides(color=guide_legend(nrow=2))
  
  d1$p2<-d1$p-mean(d1$p)
  d2$p2<-d2$p-mean(d2$p)
  d<-rbind(d1,d2)
  d$p<-d$p2*100
  d$p_se<- (sqrt(d1$p*(1-d1$p)/d1$N + d2$p*(1-d2$p)/d2$N))*100
  d1<-d[d$group==pheno1,]
  d2<-d[d$group==pheno2,]
  #p.values <- fisher.test(fite(bar[1,1],bar[2,1]))$p.value
  p.values <- fisher.test(data.frame("V1"=c(sum(d1$n),(sum(d1$N)-sum(d1$n))),"V2"=c(sum(d2$n),(sum(d2$N)-sum(d2$n)))))$p.value
  labels <- symnum(p.values, corr = FALSE, cutpoints = c(0,  .001,.01,.05, 1), symbols = c("***","**","*","p=NS"))
  # bar plot 1
  #p1<-ggplot(data=bar, aes(x=pheno, y=V1,fill=pheno))+geom_bar(stat="identity")+theme_classic()+ylab(paste("Mean ",trait," prevalence (%)",sep=""))+ theme(legend.position="none",plot.margin=unit(c(0.5,0,0.5,0.5),"cm"),text = element_text(size=18))+scale_fill_manual(values=c("#33A2C9","#CD375D"))+geom_signif(comparisons = list(c(pheno1,pheno2)), annotations =labels,textsize = 10,vjust=0.2)+ylim(0,number)+geom_text(x=pheno1, y=max(bar$V1), label=fisher.test(fite(bar[1,1],bar[2,1]))$p.value)
  p1<-ggplot(data=bar, aes(x=pheno, y=V1,fill=pheno))+geom_bar(stat="identity",colour="black")+
    theme_classic()+ylab(paste(trait," prevalence (%)",sep=""))+ theme(legend.position="none",plot.margin=unit(c(0.5,0,0.5,3),"cm"),text = element_text(size=18))+scale_fill_manual(values=c("#33A2C9","#CD375D"))+geom_signif(y_position = max(bar$V1)+0.5,xmin=0.955,xmax=1.975,annotations =labels,textsize = 10,vjust=0.2)+ylim(0,number)
  # bar plot 2
  smot<-as.data.frame(smot)
  smot<-unfactor(smot)
  smot$V2 <- factor(smot$V2, levels = c("Low risk","Medium risk","High risk"))
  smot$V1 <- factor(smot$V1, levels = c(pheno1,pheno2))
  smot$V3<-as.numeric(smot$V3)*100
  diff_df = smot %>%
    group_by(V2) %>%
    spread(V1, V3) %>%
    mutate(diff = !!as.name(pheno2) - !!as.name(pheno1),
           max_y = max(!!as.name(pheno2), !!as.name(pheno1)))
  diff_df$diff<-round(diff_df$diff,digits = 2)
  p.values <- c(fisher.test(data.frame("V1"=c(sum(d1$n[1]),(sum(d1$N[1])-sum(d1$n[1]))),"V2"=c(sum(d2$n[1]),(sum(d2$N[1])-sum(d2$n[1])))))$p.value,fisher.test(data.frame("V1"=c(sum(d1$n[2:9]),(sum(d1$N[2:9])-sum(d1$n[2:9]))),"V2"=c(sum(d2$n[2:9]),(sum(d2$N[2:9])-sum(d2$n[2:9])))))$p.value,fisher.test(data.frame("V1"=c(sum(d1$n[10]),(sum(d1$N[10])-sum(d1$n[10]))),"V2"=c(sum(d2$n[10]),(sum(d2$N[10])-sum(d2$n[10])))))$p.value)
  #p.values <- fisher.test(data.frame("V1"=c(sum(d1$n),(sum(d1$N)-sum(d1$n))),"V2"=c(sum(d2$n),(sum(d2$N)-sum(d2$n)))))$p.value
  labels <- symnum(p.values, corr = FALSE, cutpoints = c(0,  .001,.01,.05, 1), symbols = c("***","**","*","p=NS"))
  y.values <- sapply(split(smot, smot$V2), function(x){max(sapply(split(x, x$V1), function(xx){boxplot(x$V3, plot=F)$stats[5, ]}))})+1.5
  sigpos=data.frame(x=c(0.775, 1.775,2.775), xend=c(1.225, 2.225,3.225),y=y.values, annotation=labels)
  cols <- c(pheno1="#33A2C9",pheno2="#CD375D","Difference in prevalence (%)"="#cccccc")
  names(cols)[1:2]<-c(pheno1,pheno2)
  p2<-ggplot(smot,aes(V2,V3))+geom_bar(aes(y = max_y, fill = "Difference in prevalence (%)"), data = diff_df, stat = "identity", width = 0.4)+geom_bar(aes(fill=V1),stat="identity",colour="black",position="dodge")+geom_text(aes(label = diff, y = max_y), vjust=-0.5, data = diff_df,hjust = 1, colour = scales::muted("red")) +theme_classic()+xlab("wGRS")+ylab(paste(trait," prevalence (%)",sep=""))+theme(legend.title = element_blank(),legend.position = c(0.1,1),plot.margin=unit(c(0,7,0.5,1.5),"cm"),legend.direction = "horizontal",text = element_text(size=18),legend.text=element_text(size=14))+scale_fill_manual(values=cols,breaks = names(cols))+ scale_y_continuous(breaks = seq(0, number, by = 5))+geom_signif(y_position = y.values,xmin=sigpos$x,xmax=sigpos$xend,annotations = sigpos$annotation,textsize=10,vjust=0.2) 
  # top, right, bottom, and left margin 1,2,0.5,2
  # par bottom, left, top, right
  
  modelLinear11 <- lm(p~numAllele,data=d1,weights = N)
  modelLinear12 <- lm(p~1+numAllele+numAllele2,data=d1,weights = 1/p_se)
  predLin11     <- predict(modelLinear11,type="response")
  predLin12     <- predict(modelLinear12,type="response")
  Rlin11        <- round(unlist(cor.test(d1$p,predLin11)[c("estimate","conf.int")]),2)
  Rlin12        <- round(unlist(cor.test(d1$p,predLin12)[c("estimate","conf.int")]),2)
  modelLinear21 <- lm(p~numAllele,data=d2,weights = N)
  modelLinear22 <- lm(p~1+numAllele+numAllele2,data=d2,weights = 1/p_se)
  predLin21     <- predict(modelLinear21,type="response")
  predLin22     <- predict(modelLinear22,type="response")
  Rlin21        <- round(unlist(cor.test(d2$p,predLin21)[c("estimate","conf.int")]),2)
  Rlin22        <- round(unlist(cor.test(d2$p,predLin22)[c("estimate","conf.int")]),2)
  d1$predLin<-predLin12
  d2$predLin<-predLin22
  d<-rbind(d1,d2)
  pol<-data.frame(x=c(quant1,quant1,quant2,quant2),y=c(min(d$p),max(d$p),max(d$p),min(d$p)))
  p3<-ggplot()+geom_polygon(data=pol, mapping=aes(x=x, y=y),fill="grey",alpha=0.4)+geom_point(data=d,aes(x=numAllele,y=p,color=group),pch=19)+geom_line(data=d,aes(x=numAllele,y=predLin,group=group,color=group),size=1.5)+theme_classic()+ylab(paste("Difference prevalence (%) in each decile compared\n to the average ",trait," prevalence of respective group ",sep=""))+geom_errorbar(data=d,aes(x=numAllele,ymin=p-p_se, ymax=p+p_se,color=group), size=0.3,width=0.1)+xlab("wGRS")+ theme(legend.title =element_blank(),legend.position="top",legend.text=element_text(size=13),plot.margin=unit(c(0.5,0.5,0.5,0.5),"cm"),text = element_text(size=18))+ scale_color_manual(values=cols[1:2],breaks = names(cols[1:2]),labels=c(paste0(pheno1,": Quadratic: R=",Rlin12[1]," (95%CI: [",Rlin12[2],"-",Rlin12[3],"])"),paste0(pheno2,": Quadratic: R=",Rlin22[1]," (95%CI: [",Rlin22[2],"-",Rlin22[3],"])")))+guides(color=guide_legend(nrow=2))
  
  d1<-dn[dn$group==pheno1,]
  d2<-dn[dn$group==pheno2,]
  d<-d1
  d$p<-(d2$p-d1$p)*100
  d$N<-d1$N+d2$N
  d$p_se<-(sqrt(d1$p*(1-d1$p)/d1$N + d2$p*(1-d2$p)/d2$N))*100
  d$numAllele2<-d$numAllele^2
  modelLinear2 <- lm(p~1+numAllele+numAllele2,data=d,weights = 1/p_se)
  predLin2     <- predict(modelLinear2,type="response")
  Rlin2        <- round(unlist(cor.test(d$p,predLin2)[c("estimate","conf.int")]),2)
  d$predLin<-predLin2
  pol<-data.frame(x=c(quant1,quant1,quant2,quant2),y=c(min(d$p),max(d$p),max(d$p),min(d$p)))
  p4<-ggplot()+geom_polygon(data=pol, mapping=aes(x=x, y=y),fill="grey",alpha=0.4)+geom_point(data=d,aes(x=numAllele,y=p),pch=19)+geom_line(data=d,aes(x=numAllele,y=predLin,color=group),size=1.5)+theme_classic()+ylab(paste("Difference in ",trait, "prevalence (%) in each decile\n between subjects ",pheno1," and ",pheno2,sep=""))+geom_errorbar(data=d,aes(x=numAllele,ymin=p-p_se, ymax=p+p_se), size=0.3,width=0.1)+xlab("wGRS")+ scale_color_manual(values="coral",labels=c(paste0("Quadratic: R=",Rlin2[1]," (95%CI: [",Rlin2[2],"-",Rlin2[3],"])")))+theme(legend.title =element_blank(),legend.position="top",legend.text=element_text(size=16),plot.margin=unit(c(0.5,1,0.5,0.5),"cm"),text = element_text(size=18))+scale_y_continuous(position = "left",limits=c(min(d$p-d$p_se),max(d$p+d$p_se)))+xlim(c(min(d$numAllele-sd(d$numAllele)),max(d$numAllele+sd(d$numAllele))))
  myplot1 <- arrangeGrob(p1, top = textGrob("A", x = unit(0, "npc"), y   = unit(1, "npc"), just=c("left","top"),gp=gpar(col="black", fontsize=18)))
  myplot2 <- arrangeGrob(p2, top = textGrob("B", x = unit(1, "npc"), y = unit(1, "npc"), just=c("right","top"),gp=gpar(col="black", fontsize=18)))
  myplot3 <- arrangeGrob(p3, top = textGrob("C", x = unit(0, "npc"), y  = unit(1, "npc"), just=c("left","top"),gp=gpar(col="black", fontsize=18)))
  myplot4 <- arrangeGrob(p0, top = textGrob("D", x = unit(0, "npc"), y = unit(1, "npc"), just=c("left","top"),gp=gpar(col="black",  fontsize=18)))
  myplot5 <- arrangeGrob(p4, top = textGrob("E", x = unit(0, "npc"), y = unit(1, "npc"), just=c("left","top"),gp=gpar(col="black",  fontsize=18)))
  cairo_pdf(paste(figure,"_T.pdf",sep=""),width=17,height=16)
  grid.arrange(arrangeGrob(myplot1,myplot2, ncol=2, nrow=1,widths=c(1,3)),arrangeGrob(myplot4,myplot3,myplot5, ncol=3, nrow=1),nrow=2)
  dev.off()
}

## example ggplot2
pol<-data.frame(x=c(181.685,181.685,203.873,203.873),y=c(0,max(d$p),max(d$p),0))
p<-ggplot()+geom_polygon(data=pol, mapping=aes(x=x, y=y),fill="grey",alpha=0.4)+geom_point(data=d,aes(x=numAllele,y=p,color=group),pch=19)+geom_line(data=d,aes(x=numAllele,y=predLogit,group=group,color=group),size=1.5)+theme_classic()+ylab("Difference prevalence compared to the average in CAD (%)")+geom_errorbar(data=d,aes(x=numAllele,ymin=p-p_se, ymax=p+p_se,color=group), size=0.3,width=1)+xlab("Mean number of risk alleles")+ theme(legend.title =element_blank(),legend.position="top",legend.text=element_text(size=15))+ scale_color_manual(values=c("#CD375D","#33A2C9"))
q<-ggdraw(add_sub(p, paste("Never_smoked\n",
                           paste0("Logit: R=",Rlogit1[1]," (95%CI: [",Rlogit1[2],"-",Rlogit1[3],"])"),"\n",
                           paste0("Log: R=",Rlog1[1]," (95%CI: [",Rlog1[2],"-",Rlog1[3],"])"),"\n",
                           paste0("Probit: R=",Rprobit1[1]," (95%CI: [",Rprobit1[2],"-",Rprobit1[3],"])"),"\n",
                           paste0("Linear: R=",Rlin1[1]," (95%CI: [",Rlin1[2],"-",Rlin1[3],"])")),size=10,x=0,y=0.5,hjust=0))

ggdraw(add_sub(q, paste("Ever_smoked\n",
                        paste0("Logit: R=",Rlogit2[1]," (95%CI: [",Rlogit2[2],"-",Rlogit2[3],"])"),"\n",
                        paste0("Log: R=",Rlog2[1]," (95%CI: [",Rlog2[2],"-",Rlog2[3],"])"),"\n",
                        paste0("Probit: R=",Rprobit2[1]," (95%CI: [",Rprobit2[2],"-",Rprobit2[3],"])"),"\n",
                        paste0("Linear: R=",Rlin2[1]," (95%CI: [",Rlin2[2],"-",Rlin2[3],"])")),size=10,x=0.5,y=1.58,hjust=0))


cairo_ps(filename='CAD_smoke_uGRS_gg.pdf', width=5.78, height=7.13,pointsize = 12,family="Univers")
dev.off()
