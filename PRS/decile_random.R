library(data.table)
library(tidyverse)
decilepoint<-function(dat){
  quant<-quantile(dat$NRA,probs = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
  dat$uorder<-10
  dat$uorder[dat$NRA<=quant[1]]<-1
  dat$uorder[quant[1]<dat$NRA & dat$NRA<=quant[2]]<-2
  dat$uorder[quant[2]<dat$NRA & dat$NRA<=quant[3]]<-3
  dat$uorder[quant[3]<dat$NRA & dat$NRA<=quant[4]]<-4
  dat$uorder[quant[4]<dat$NRA & dat$NRA<=quant[5]]<-5
  dat$uorder[quant[5]<dat$NRA & dat$NRA<=quant[6]]<-6
  dat$uorder[quant[6]<dat$NRA & dat$NRA<=quant[7]]<-7
  dat$uorder[quant[7]<dat$NRA & dat$NRA<=quant[8]]<-8
  dat$uorder[quant[8]<dat$NRA & dat$NRA<=quant[9]]<-9
  dat2l<-split(dat,dat$uorder)
  q<-data.frame("N"=0,"n"=0,"p"=0,"numAllele"=0)
  for (i in names(table(dat$uorder))){
    t<-dat2l[[i]]
    q[i,1:4]<-c(nrow(t),sum(t$CC),sum(t$CC)/nrow(t),mean(t$NRA))
  }
  return(q)
}

randomR<-function(trait,n,i){
  PRS<-fread(paste0("UKBB_",trait,"_",n,"_",i,".profile"))
  dat<-PRS %>% select(1,6,3)
  dat$CC<-dat$PHENO-1
  dat<-dat[,-3]
  colnames(dat)<-c("sample","NRA","CC")
  d<-decilepoint(dat)
 
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
  
  cor<-data.frame("R"=0,"Rdown"=0,"Rup"=0,"slope"=0,"Function"=0,"trait"=0,"nSNP"=0)
  cor[1,]<-c(Rlogit[1],Rlogit[2],Rlogit[3],modelLogit$coefficients[2],"Logit",trait,n)
  cor[2,]<-c(Rprobit[1],Rprobit[2],Rprobit[3],modelProbit$coefficients[2],"Probit",trait,n)
  cor[3,]<-c(Rlog[1],Rlog[2],Rlog[3],modelLog$coefficients[2],"Log",trait,n)
  cor[4,]<-c(Rlin[1],Rlin[2],Rlin[3],modelLinear$coefficients[2],"Lin",trait,n)
  fwrite(cor,paste0(trait,"_",n,"_",i,"_decile.random"),col.names = FALSE)
}

args<-commandArgs(TRUE)
randomR(args[1],args[2],args[3])
