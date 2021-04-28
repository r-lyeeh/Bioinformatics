install.packages("ROSE")
library(ROSE)
## imbalanced data, 1 random selection of samples grenerate an AUC = 0.6033 which is poorly fitted
## try python script to deal with this problem
t<-caddat2 %>% select(1,2,3)
prop.table(table(t$CC))
t.temp<-t[1:nrow(t),]
#over sampling 0.18
data_balanced_over <- ovun.sample(CC~NRA, data = t.temp, method = "over",N = 406996*2,seed=1)$data
table(data_balanced_over$CC)
modelLogit  <- glm(CC~NRA,data=data_balanced_over,family = binomial(link="logit"))
predLogit   <- predict(modelLogit,type = "response")
Rlogit1      <- round(unlist(cor.test(data_balanced_over$CC,predLogit)[c("estimate","conf.int")]),2)
(roc.fit1 <- roc(response = modelLogit$data$CC, predictor = modelLogit$fitted.values, direction = "<", ci = TRUE, ci.method = "boot", boot.n = 10))
rocobj1 <- roc(data_balanced_over$CC, predLogit)

# raw 0.18
t.case<-t[t$CC==1,]
t.control<-t[t$CC==0,]
t.control<-t.control[sample(nrow(t.control), nrow(t.case)), ]
data_raw<-rbind(t.case,t.control)
modelLogit  <- glm(CC~NRA,data=data_raw,family = binomial(link="logit"))
predLogit   <- predict(modelLogit,type = "response")
Rlogit      <- round(unlist(cor.test(data_raw$CC,predLogit)[c("estimate","conf.int")]),2)
(roc.fit1 <- roc(response = modelLogit$data$CC, predictor = modelLogit$fitted.values, direction = "<", ci = TRUE, ci.method = "boot", boot.n = 10))


# under sampling 0.18
data_balanced_under <- ovun.sample(CC~NRA, data = t.temp, method = "under", N = 20310*2, seed = 1)$data
table(data_balanced_under$CC)
modelLogit  <- glm(CC~NRA,data=data_balanced_under,family = binomial(link="logit"))
predLogit   <- predict(modelLogit,type = "response")
Rlogit2      <- round(unlist(cor.test(data_balanced_under$CC,predLogit)[c("estimate","conf.int")]),2)
(roc.fit1 <- roc(response = modelLogit$data$CC, predictor = modelLogit$fitted.values, direction = "<", ci = TRUE, ci.method = "boot", boot.n = 10))
rocobj2 <- roc(data_balanced_under$CC, predLogit)

# both sampling 0.18
data_balanced_both <- ovun.sample(CC ~ NRA, data = t.temp, method = "both", p=0.5, N=213653, seed = 1)$data
table(data_balanced_both$CC)
modelLogit  <- glm(CC~NRA,data=data_balanced_both,family = binomial(link="logit"))
predLogit   <- predict(modelLogit,type = "response")
Rlogit3      <- round(unlist(cor.test(data_balanced_both$CC,predLogit)[c("estimate","conf.int")]),2)
(roc.fit1 <- roc(response = modelLogit$data$CC, predictor = modelLogit$fitted.values, direction = "<", ci = TRUE, ci.method = "boot", boot.n = 10))
rocobj3 <- roc(data_balanced_both$CC, predLogit)

# rose sampling 0.18
data.rose <- ROSE(CC~NRA, data = t.temp, seed = 1)$data
table(data.rose$CC)
modelLogit  <- glm(CC~NRA,data=data.rose,family = binomial(link="logit"))
predLogit   <- predict(modelLogit,type = "response")
Rlogit4      <- round(unlist(cor.test(data.rose$CC,predLogit)[c("estimate","conf.int")]),2)
(roc.fit1 <- roc(response = modelLogit$data$CC, predictor = modelLogit$fitted.values, direction = "<", ci = TRUE, ci.method = "boot", boot.n = 10))
rocobj4 <- roc(data.rose$CC, predLogit)

# SMOTE sampling 0.20
library(DMwR)
t.temp0<-t.temp
t.temp0$CC<-as.factor(t.temp0$CC)
#trainSplit<-SMOTE(CC~NRA,t.temp0,perc.over=20310,perc.under=406996)
data.smote <- SMOTE(CC~NRA, data = t.temp0)
table(data.smote$CC)
modelLogit  <- glm(CC~NRA,data=data.smote,family = binomial(link="logit"))
predLogit   <- predict(modelLogit,type = "response")
data.smote$CC<-as.numeric(data.smote$CC)
Rlogit5      <- round(unlist(cor.test(data.smote$CC,predLogit)[c("estimate","conf.int")]),2)
(roc.fit1 <- roc(response = modelLogit$data$CC, predictor = modelLogit$fitted.values, direction = "<", ci = TRUE, ci.method = "boot", boot.n = 10))
rocobj5 <- roc(data.smote$CC, predLogit)

g2 <- ggroc(list(under=rocobj1, over=rocobj2, both=rocobj3,ROSE=rocobj4,SMOTE=rocobj5))+theme_classic()


# ggplot Logistic regression
ggplot(data=data.smote,aes(x=NRA,y=CC))+geom_point()+geom_smooth(method = "glm",  formula = y~x, method.args = list(family = "binomial"), se = FALSE)
ggplot(data=data.smote,aes(x=NRA,y=CC))+geom_point()+geom_smooth(method = "glm",  formula = y~log(x), method.args = list(family = "binomial"), se = FALSE)
ggplot(data=data.smote,aes(x=NRA,y=CC))+geom_point()+geom_smooth(method = "glm",  formula = y~x, method.args = list(family = binomial(link="logit")), se = FALSE)+theme_classic()+xlab("Number of risk alleles")+ylab("Case/Control")+ggtitle("CAD 198 SNPs")
cbl<-cbind(data.smote$NRA,predLogit)
ggplot(data=data.smote,aes(x=NRA,y=CC))+geom_point()+geom_line(data=cbl,aes(x=NRA,y=predLogit))

# decile random
test<-fread("../decile/CAD_150_decile_sum.random")
colnames(test)<-c("R","Rdown","Rup","Slope","Function","Trait","nSNP")
test$Function<-factor(test$Function,levels=c("Logit","Log","Probit","Lin"))
ggplot(test)+geom_boxplot(aes(x=Function,y=R,fill=Function))+theme_classic()+ylab("Correlation R")+ggtitle("CAD 150 SNPs (decile)")
ggplot(test)+geom_boxplot(aes(x=Function,y=R,fill=Function),outlier.shape = NA)+theme_classic()+ylab("Correlation R")+ggtitle("CAD 150 SNPs (decile)")+ylim(0.97,1)

# CC random
test<-fread("../CC-R/CAD_150_CC_sum.random")
colnames(test)<-c("R","Rdown","Rup","Slope","AUC","Function","Trait","nSNP")
test$Function<-factor(test$Function,levels=c("Logit","Log","Probit","Lin"))
ggplot(test)+geom_boxplot(aes(x=Function,y=R,fill=Function))+theme_classic()+ylab("Correlation R")+ggtitle("CAD 150 SNPs (Case/Control)")
ggplot(test)+geom_boxplot(aes(x=Function,y=AUC,fill=Function),outlier.shape = NA)+theme_classic()+ylab("AUC")+ggtitle("CAD 150 SNPs (Case/Control)")

# SMOTE random
test<-fread("../SMOTE/SMOTE_sum.random")
colnames(test)<-c("R","Rdown","Rup","Slope","AUC","Function","Trait","nSNP")
test$Function<-factor(test$Function,levels=c("Logit","Log","Probit","Lin"))
ggplot(test)+geom_boxplot(aes(x=Function,y=R,fill=Function))+theme_classic()+ylab("Correlation R")+ggtitle("CAD 10 SNPs (Case/Control, SMOTE)")
ggplot(test)+geom_boxplot(aes(x=Function,y=AUC,fill=Function),outlier.shape = NA)+theme_classic()+ylab("AUC")+ggtitle("CAD 10 SNPs (Case/Control, SMOTE)")

## figures
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

## 10 SNPs -> seperate to equal decile
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
    print(i)
    t<-dat2l[[i]]
    q[i,1:4]<-c(nrow(t),sum(t$CC),sum(t$CC)/nrow(t),mean(t$NRA))
  }
  return(q)
}

## 10 SNPs -> seperate to equal sample size
decilepoint<-function(dat){
  dat$uorder<-10
  dat$uorder[1:(nrow(dat)/10)]<-1
  dat$uorder[(nrow(dat)/10):(2*(nrow(dat)/10))]<-2
  dat$uorder[(2*(nrow(dat)/10)):(3*(nrow(dat)/10))]<-3
  dat$uorder[(3*(nrow(dat)/10)):(4*(nrow(dat)/10))]<-4
  dat$uorder[(4*(nrow(dat)/10)):(5*(nrow(dat)/10))]<-5
  dat$uorder[(5*(nrow(dat)/10)):(6*(nrow(dat)/10))]<-6
  dat$uorder[(6*(nrow(dat)/10)):(7*(nrow(dat)/10))]<-7
  dat$uorder[(7*(nrow(dat)/10)):(8*(nrow(dat)/10))]<-8
  dat$uorder[(8*(nrow(dat)/10)):(9*(nrow(dat)/10))]<-9
  dat2l<-split(dat,dat$uorder)
  q<-data.frame("N"=0,"n"=0,"p"=0,"numAllele"=0)
  for (i in names(table(dat$uorder))){
    print(i)
    t<-dat2l[[i]]
    q[i,1:4]<-c(nrow(t),sum(t$CC),sum(t$CC)/nrow(t),mean(t$NRA))
  }
  return(q)
}


test<-fread("t1")
colnames(test)<-c("R","Rdown","Rup","Function","group")
test$group<-factor(test$group,levels=c("CAD_10","CAD_50","CAD_99","CAD_150","BC_10","BC_50","BC_100","PC_10","PC_50","PC_100"))
test$Function<-factor(test$Function,levels=c("Logit","Log","Probit","Lin"))
ggplot(data=test,aes(x=group,y=R,fill=Function))+geom_boxplot(outlier.shape = NA)+theme_classic()+ylim(0.85,1)+ylab("Correlation R")

test1<-test %>%select("R","Function","group")
test2<-test %>%select("Rdown","Function","group")
colnames(test2)[1]<-"R"
test1<-rbind(test1,test2)
test2<-test %>%select("Rup","Function","group")
colnames(test2)[1]<-"R"
test1<-rbind(test1,test2)
ggplot(data=test1,aes(x=group,y=R,fill=Function))+geom_boxplot(outlier.shape = NA)+theme_classic()+ylim(0.85,1)+ylab("Correlation R")

ggplot(data=test00,aes(x=Function,y=R,fill=Function))+geom_boxplot(outlier.shape = NA)+theme_classic()+ylim(0.85,1)+ylab("Correlation R")

test$group<-factor(test$group,levels=c("CAD_10","CAD_50","CAD_99","CAD_150"))
test$Function<-factor(test$Function,levels=c("Logit","Log","Probit","Lin"))
ggplot(data=test,aes(x=group,y=R,fill=Function))+geom_boxplot(outlier.shape = NA)+theme_classic()+ylim(0.85,1)+ylab("Correlation R")

test<-fread("Version_50/50sample.random")
test<-fread("Version_200/Function_R.random")

test1<-fread("profile/BC_tf.random")
test2<-fread("profile/PC_tf.random")
test3<-fread("profile/CAD_tf.random")

normat<-function(test){
  colnames(test)<-c("R","Rdown","Rup","slope","Function","group")
  t1<-test[test$Function=="Probit",]
  t1$slope<-scale(t1$slope)
  t2<-test[test$Function=="Log",]
  t2$slope<-scale(t2$slope)
  t3<-test[test$Function=="Logit",]
  t3$slope<-scale(t3$slope)
  t4<-test[test$Function=="Lin",]
  t4$slope<-scale(t4$slope)  
  t<-rbind(t1,t2,t3,t4)
  return(t)
}
test10<-normat(test1)
test20<-normat(test2)
test30<-normat(test3)
test<-rbind(test10,test20,test30)
test$group<-factor(test$group,levels=c("CAD_10","CAD_50","CAD_150","BC_10","BC_50","BC_100","PC_10","PC_50","PC_100"))
test$Function<-factor(test$Function,levels=c("Logit","Log","Probit","Lin"))
ggplot(data=test,aes(x=group,y=slope,fill=Function))+geom_boxplot(outlier.shape = NA)+theme_classic()+ylab("Slope")

test<-fread("profile/testf.random")
colnames(test)<-c("R","Rdown","Rup","slope","Function","group")
test<-test[test$Function=="Lin",]
test$group<-factor(test$group,levels=c("CAD_10","CAD_50","CAD_150","BC_10","BC_50","BC_100","PC_10","PC_50","PC_100"))
test$Function<-factor(test$Function,levels=c("Logit","Log","Probit","Lin"))
ggplot(data=test,aes(x=group,y=slope))+geom_boxplot(outlier.shape = NA)+theme_classic()+ylab("Slope (Linear)")
