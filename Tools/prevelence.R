setwd("~/Desktop/Papers/LRT_Heri_paper/DataFeb2020/")
load("prs.RData")
load("../New_Data/diseases_ukb.RData")
Diseases <- d[,c("IID","BreastCancer","CAD","T2D","IBD")]
#pc <- read.table("RunID_37403_Diagnosis_C61_RunID_31139_death_C61.pheno",h=T,stringsAsFactors = F)[,c(2,6)]
#colnames(pc) <- c("IID","ProstateCancer")
#Diseases <- merge(Diseases,pc,by="IID")
m  <- merge(Diseases,prs,by="IID")

analyse <- function(ds){
  dat      <- m[,c(ds,paste0("PRS_",ds))]; colnames(dat) <- c("CC","NRA")
  modComp  <- summary( glm(CC~NRA,family="binomial",data=dat) )$coefficients
  dat$grp  <- round(dat$NRA) #cut(dat$NRA,breaks = seq(min(dat$NRA),max(dat$NRA),by=2),include.lowest=TRUE)
  prev_p   <- aggregate(dat$CC,by=list(dat$grp),FUN=mean,na.rm=T)
  prev_n   <- aggregate(dat$CC,by=list(dat$grp),FUN=sum,na.rm=T)
  prev_N   <- aggregate(dat$CC,by=list(dat$grp),FUN=function(x) length(!is.na(x)))
  allele_c <- aggregate(dat$NRA,by=list(dat$grp),FUN=mean)
  d        <- cbind.data.frame(N=prev_N$x,n=prev_n$x,p=prev_p$x,numAllele=allele_c$x)
  d$p_se   <- sqrt(d$p*(1-d$p)/d$N)
  d        <- d[which(d$N>=10),]
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
  
  Cols <- c("seagreen","coral1","goldenrod")
  l    <- which(d[,"p_se"]>0)
  
  png(paste0("plots/",ds,".prevalence.fit.png"),width=1500,height=1250,res=250)
  plot(d[l,"numAllele"],d[l,"p"],pch=19,axes=FALSE,ylim=c(0,max(d$p+d$p_se)),main=ds,
       xlab="Mean number of risk alleles",ylab="Prevalence");axis(1);axis(2)
  for(i in 1:nrow(d)){
    if(d[i,"p_se"]>0){
      arrows(d[i,"numAllele"],d[i,"p"]-d[i,"p_se"],d[i,"numAllele"],d[i,"p"]+d[i,"p_se"],col="grey",angle=90,len=0.05)
      arrows(d[i,"numAllele"],d[i,"p"]+d[i,"p_se"],d[i,"numAllele"],d[i,"p"]-d[i,"p_se"],col="grey",angle=90,len=0.05)
    }
  }
  matlines(d$numAllele,cbind(predLogit,predProbit,predLog),col=Cols,lwd=2)
  abline(a=modelLinear$coefficients[1],b=modelLinear$coefficients[2],col=4,lty=3)
  legend("topleft",
         legend=c(paste0("Logit: R=",Rlogit[1]," (95%CI: [",Rlogit[2],"-",Rlogit[3],"])"),
                  paste0("Probit: R=",Rprobit[1]," (95%CI: [",Rprobit[2],"-",Rprobit[3],"])"),
                  paste0("Log: R=",Rlog[1]," (95%CI: [",Rlog[2],"-",Rlog[3],"])"),
                  paste0("Linear: R=",Rlin[1]," (95%CI: [",Rlin[2],"-",Rlin[3],"])")),
         col=c(Cols,"blue"),cex=0.9,box.lty=0,lty=c(1,1,1,3),lwd=2)
  dev.off()
  
  ## Second  plot
  png(paste0("plots/",ds,".prevalence_per_allele.png"),width=1500,height=1250,res=250)
  plot(d[l,"numAllele"],d[l,"px"],pch=19,axes=FALSE,ylim=c(0,max(d$px+d$px_se)),main=ds,
       xlab="Mean number of risk alleles",ylab="Prevalence/Allele");axis(1);axis(2)
  for(i in 1:nrow(d)){
    if(d[i,"px_se"]>0){
      arrows(d[i,"numAllele"],d[i,"px"]-d[i,"px_se"],d[i,"numAllele"],d[i,"px"]+d[i,"px_se"],col="grey",angle=90,len=0.05)
      arrows(d[i,"numAllele"],d[i,"px"]+d[i,"px_se"],d[i,"numAllele"],d[i,"px"]-d[i,"px_se"],col="grey",angle=90,len=0.05)
    }
  }
  md <- lm(d[l,"px"]~d[l,"numAllele"],weights = d[l,"N"])
  mx <- coef(md)
  abline(a=mx[1],b=mx[2],col="grey",lty=2)
  
  modelLinear <- lm(px~numAllele,data=d,weights = N)
  predLin     <- predict(modelLinear,type="response")
  Rlin        <- round(unlist(cor.test(d$p,predLin)[c("estimate","conf.int")]),2)
  legend("topleft",legend=paste0("Linear: R=",Rlin[1]," (95%CI: [",Rlin[2],"-",Rlin[3],"])"),box.lty=0,lty=2,col="grey")
  dev.off()
  
  results <- c(modComp[2,c(1,4)],cor_prev_nra=cor(d$p,d$numAllele))
  names(results) <- c("beta_allele_count","p_allele_count","corr_prev_allele_count")
  return(results)
}
diseases <- colnames(Diseases)[-1]
A        <- do.call("rbind",lapply(diseases,analyse)); rownames(A) <- diseases

