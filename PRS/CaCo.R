library(data.table)
library(tidyverse)
library(pROC)
#library(ROSE)
#library(DMwR)
randomR<-function(trait,n,i){
  PRS<-fread(paste0("UKBB_",trait,"_",n,"_",i,".profile"))
  dat<-PRS %>% select(1,6,3)
  dat$CC<-dat$PHENO-1
  dat<-dat[,-3]
  colnames(dat)<-c("sample","NRA","CC")
  dat.case<-dat[dat$CC==1,]
  dat.control<-dat[dat$CC==0,]
  dat.control<-dat.control[sample(nrow(dat.control),nrow(dat.case)),]
  dat<-rbind(dat.case,dat.control)

  modelLogit  <- glm(CC~NRA,data=dat,family = binomial(link="logit"))
  modelProbit <- glm(CC~NRA,data=dat,family = binomial(link="probit"))
  #modelLog    <- glm(CC~NRA,data=dat,family = binomial(link="log"))
  modelLinear <- lm(CC~NRA,data=dat)
  
  predLogit   <- predict(modelLogit,type = "response")
  predProbit  <- predict(modelProbit,type = "response")
  #predLog     <- predict(modelLog,type = "response")
  predLin     <- predict(modelLinear,type="response")
  
  Rlogit      <- round(unlist(cor.test(dat$CC,predLogit)[c("estimate","conf.int")]),2)
  Rprobit     <- round(unlist(cor.test(dat$CC,predProbit)[c("estimate","conf.int")]),2)
  #Rlog        <- round(unlist(cor.test(dat$CC,predLog)[c("estimate","conf.int")]),2)
  Rlin        <- round(unlist(cor.test(dat$CC,predLin)[c("estimate","conf.int")]),2)
  
  roc.fitLogit<-roc(response = modelLogit$data$CC, predictor = modelLogit$fitted.values, direction = "<", ci = TRUE, ci.method = "boot", boot.n = 100)
  roc.fitProbit<-roc(response = modelProbit$data$CC, predictor = modelProbit$fitted.values, direction = "<", ci = TRUE, ci.method = "boot", boot.n = 100)
  roc.fitLin<-roc(response = dat$CC, predictor = modelLinear$fitted.values, direction = "<", ci = TRUE, ci.method = "boot", boot.n = 100)
  
  cor<-data.frame("R"=0,"Rdown"=0,"Rup"=0,"slope"=0,"AUC"=0,"Function"=0,"trait"=0,"nSNP"=0)
  cor[1,]<-c(Rlogit[1],Rlogit[2],Rlogit[3],modelLogit$coefficients[2],roc.fitLogit$auc[1],"Logit",trait,n)
  cor[2,]<-c(Rprobit[1],Rprobit[2],Rprobit[3],modelProbit$coefficients[2],roc.fitProbit$auc[1],"Probit",trait,n)
  #cor[3,]<-c(Rlog[1],Rlog[2],Rlog[3],modelLog$coefficients[2],"Log",trait,n)
  cor[3,]<-c(Rlin[1],Rlin[2],Rlin[3],modelLinear$coefficients[2],roc.fitProbit$auc[1],"Lin",trait,n)
  fwrite(cor,paste0(trait,"_",n,"_",i,"_CC.random"),col.names = FALSE)
}

args<-commandArgs(TRUE)
randomR(args[1],args[2],args[3])
