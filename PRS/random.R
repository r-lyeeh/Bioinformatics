library(data.table)
library(tidyverse)
library(rsq)
randomR1<-function(trait,n,i){
  PRS<-fread(paste0("UKBB_",trait,"_",n,"_",i,".profile"))
  dat<-PRS %>% select(1,6,3)
  dat$CC<-dat$PHENO-1
  dat<-dat[,-3]
  colnames(dat)<-c("sample","NRA","CC")
  modComp  <- summary( glm(CC~NRA,family="binomial",data=dat) )$coefficients
  dat$grp  <- round(dat$NRA) #cut(dat$NRA,breaks = seq(min(dat$NRA),max(dat$NRA),by=2),include.lowest=TRUE)
  prev_p   <- aggregate(dat$CC,by=list(dat$grp),FUN=mean,na.rm=T)
  prev_n   <- aggregate(dat$CC,by=list(dat$grp),FUN=sum,na.rm=T)
  prev_N   <- aggregate(dat$CC,by=list(dat$grp),FUN=function(x) length(!is.na(x)))
  allele_c <- aggregate(dat$NRA,by=list(dat$grp),FUN=mean)
  d        <- cbind.data.frame(N=prev_N$x,n=prev_n$x,p=prev_p$x,numAllele=allele_c$x)
  #d$p_se   <- sqrt(d$p*(1-d$p)/d$N)
  d        <- d[which(d$N>=200),]
  #d$px     <- d$p / d$numAllele
  #d$px_se  <- d$p_se / d$numAllele
  return(d)
}
randomR2<-function(trait,n,i,d){
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
  
  #AIC(), BIC() result weird 
  rsqLogit<-rsq(modelLogit,type="lr")
  rsqProbit<-rsq(modelProbit,type="lr")
  rsqLog<-rsq(modelLog,type="lr")
  rsqLinear<-rsq(modelLinear,type="lr") 

  cor<-data.frame("R"=0,"Rdown"=0,"Rup"=0,"slope"=0,"Function"=0,"trait"=0,"nSNP"=0,"rsq"=0)
  cor[1,]<-c(Rlogit[1],Rlogit[2],Rlogit[3],modelLogit$coefficients[2],"Logit",trait,n,rsqLogit)
  cor[2,]<-c(Rprobit[1],Rprobit[2],Rprobit[3],modelProbit$coefficients[2],"Probit",trait,n,rsqProbit)
  cor[3,]<-c(Rlog[1],Rlog[2],Rlog[3],modelLog$coefficients[2],"Log",trait,n,rsqLog)
  cor[4,]<-c(Rlin[1],Rlin[2],Rlin[3],modelLinear$coefficients[2],"Lin",trait,n,rsqLinear)
  fwrite(cor,paste0(trait,"_",n,"_",i,".random"),col.names = FALSE)
}

args<-commandArgs(TRUE)
d<-randomR1(args[1],args[2],args[3])
randomR2(args[1],args[2],args[3],d)
