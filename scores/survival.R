library("survival")
library("survminer")
# read data
HES<-fread("~/Desktop/Fabian/pheno/42000_DateHESCAD.ukb")
t2<-date[match(caddat2$sample,date$f.eid),]
caddat2$date<-t2$f.42000.0.0
caddat2$datedia<-t2$f.42001.0.0
assdate<-fread("~/Desktop/Fabian/pheno/53_Assess_Date.ukb")
t2<-assdate[match(caddat2$sample,assdate$f.eid),]
caddat2$assdate<-t2$f.53.0.0
ukbpc<-fread("~/Desktop//Scattered/MO_Johan_20200221/genotypes/ukb_imp_20PCs.eigenvec")
t2<-ukbpc[match(caddat2$sample,ukbpc$V1),]
cadadt2$PC1<-t2$V3
cadadt2$PC2<-t2$V4
cadadt2$PC3<-t2$V5
cadadt2$PC4<-t2$V6
cadadt2$PC5<-t2$V7
cadadt2$PC6<-t2$V8
cadadt2$PC7<-t2$V9
cadadt2$PC8<-t2$V10
cadadt2$PC9<-t2$V11
cadadt2$PC10<-t2$V12
caddat2$time<-caddat2$date-caddat2$assdate
death<-fread("~/Desktop/Fabian/pheno/40000_DateOfDeath.ukb")
t2<-death[match(CIPsurv$sample,death$f.eid),]
CIPsurv$deathdate<-t2$f.40000.1.0


# re-sampling
cad0<-caddat2[caddat2$CC==0,]
cad1<-caddat2[caddat2$CC==1,]
cad0<-cad0[is.na(cad0$date)==TRUE,]
cad1<-cad1[is.na(cad1$date)==FALSE,]
test<-rbind(cad0,cad1)
#control group time NA -large value
test$time[is.na(test$time)==TRUE]<-max(test$time,na.rm=TRUE)
#max(test$time)=3578
test<-test[test$time>=0,]

#sample size estimating
Cox.Equality <- function (alpha, beta, loghr, p1, p2, d) 
{
  n = (qnorm(1-alpha/2) + qnorm(1 - beta))^2/(log(loghr)^2 * p1 * p2 * d)
  n
}

#Cox.Equality(0.05,0.2,2,0.5,0.5,0.8)
# using equal number
cad01<-cad0[sample(nrow(cad0), 5722), ]
test<-rbind(cad01,cad1)
test$time[is.na(test$time)==TRUE]<-max(test$time,na.rm = TRUE)
test<-test[test$time>=0,]

#survival analysis
#PFS=status PFS_T = time
#单因素Cox回归分析
res.cox <- coxph(Surv(time, CC) ~ wgrs+sex+age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = test)
res.cox <- coxph(Surv(time, CC) ~ NRA+sex+age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = test)
res.cox <- coxph(Surv(time, CC) ~ qrisk3+sex+age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = test)
summary(res.cox)
survminer::ggforest(res.cox, data = test)
ggsurvplot(survfit(res.cox),data=test, color = "#2E9FDF",
           ggtheme = theme_minimal(),ylim=c(0.99,1))+ggtitle("CAD")+xlab("Time (year)")

fit <- survfit(Surv(time, CC) ~ NRA, data = test)
ggsurvplot(fit, fun="event",palette =c("#E7B800", "#2E9FDF"),data = test, ylim=c(0.99,1),risk.table = TRUE,pval=TRUE,ggtheme = theme_classic(),censor.shape="|", censor.size = 4,size=1,conf.int.alpha=c(0.1),conf.int=F,ncensor.plot = TRUE)

#kaplan-meier
#KM <- survfit(Surv(time, CC) ~ NRA+sex+age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = test, type = 'kaplan-meier', conf.type = 'log')
#summary(KM)

# Deciles
t1<-test[test$worder==1,]
t1<-test[test$worder==10,]
t1<-test[test$worder!=1 & test$worder!=10,]
t1<-test[test$qrisk3<=5,]
t1<-test[test$qrisk3>20,]
t1<-test[test$qrisk3>5 & test$qrisk3<=20,]
res.cox <- coxph(Surv(time, CC) ~ qrisk3+sex+age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = t1)

#C-index
res.cox <- coxph(Surv(time, CC) ~ qrisk3+sex+age+Diabetes+smoke+LDL+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = t1)
summary(res.cox)$concordance
res.cox <- coxph(Surv(time, CC) ~ NRA+qrisk3+sex+age+Diabetes+smoke+LDL+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = t1)
summary(res.cox)$concordance
res.cox <- coxph(Surv(time, CC) ~ wgrs+qrisk3+sex+age+Diabetes+smoke+LDL+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = t1)
summary(res.cox)$concordance
t1<-test[test$qrisk3<=5,]
t1<-test[test$qrisk3>20,]
t1<-test[test$qrisk3>5 & test$qrisk3<=20,]
t1<-test[test$sex==1,]
t1<-test[test$sex==0,]
t1<-test[test$Diabetes==1,]
t1<-test[test$Diabetes==0,]
t1<-test[test$smoke==1,]
t1<-test[test$smoke==0,]
t1<-test[test$LDL<=2.6,]
t1<-test[test$LDL>2.6 & test$LDL<=4.2,]
t1<-test[test$LDL>4.2,]
t1<-test[test$age<=44,]
t1<-test[test$age>44 & test$age<=54,]
t1<-test[test$age>54 & test$age<=64,]
t1<-test[test$age>64,]

res.cox <- coxph(Surv(time, CC) ~ qrisk3+sex+age+Diabetes+smoke+LDL+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = test)
summary(res.cox)

# C-index test
library(survcomp)
cindexp<-function(t1){
  t1<-t1%>%select("sample","NRA","CC","wgrs","uorder","worder","sex","LDL","Diabetes","smoke","age","qrisk3","time","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
  t1<-t1[complete.cases(t1)==TRUE,]
  c1 <- concordance.index(x=t1$qrisk3+t1$sex+t1$age+t1$Diabetes+t1$smoke+t1$LDL+t1$PC1+t1$PC2+t1$PC3+t1$PC4+t1$PC5+t1$PC6+t1$PC7+t1$PC8+t1$PC9+t1$PC10, surv.time=t1$time, surv.event=t1$CC, method="noether")
  c2 <- concordance.index(x=t1$wgrs+t1$qrisk3+t1$sex+t1$age+t1$Diabetes+t1$smoke+t1$LDL+t1$PC1+t1$PC2+t1$PC3+t1$PC4+t1$PC5+t1$PC6+t1$PC7+t1$PC8+t1$PC9+t1$PC10, surv.time=t1$time, surv.event=t1$CC, method="noether")
  p1<-cindex.comp(c1, c2)$p.value
  c1 <- concordance.index(x=t1$qrisk3+t1$sex+t1$age+t1$Diabetes+t1$smoke+t1$LDL+t1$PC1+t1$PC2+t1$PC3+t1$PC4+t1$PC5+t1$PC6+t1$PC7+t1$PC8+t1$PC9+t1$PC10, surv.time=t1$time, surv.event=t1$CC, method="noether")
  c2 <- concordance.index(x=t1$NRA+t1$qrisk3+t1$sex+t1$age+t1$Diabetes+t1$smoke+t1$LDL+t1$PC1+t1$PC2+t1$PC3+t1$PC4+t1$PC5+t1$PC6+t1$PC7+t1$PC8+t1$PC9+t1$PC10, surv.time=t1$time, surv.event=t1$CC, method="noether")
  p2<-cindex.comp(c1, c2)$p.value
  return(c(p1,p2))
}

cindp<-data.frame("p1"=0,"p2"=0)
t1<-test[test$qrisk3<=5,]
cindp[1,]<-cindexp(t1)
t1<-test[test$qrisk3>5 & test$qrisk3<=20,]
cindp[2,]<-cindexp(t1)
t1<-test[test$qrisk3>20,]
cindp[3,]<-cindexp(t1)
t1<-test[test$sex==1,]
cindp[4,]<-cindexp(t1)
t1<-test[test$sex==0,]
cindp[5,]<-cindexp(t1)
t1<-test[test$Diabetes==1,]
cindp[6,]<-cindexp(t1)
t1<-test[test$Diabetes==0,]
cindp[7,]<-cindexp(t1)
t1<-test[test$smoke==1,]
cindp[8,]<-cindexp(t1)
t1<-test[test$smoke==0,]
cindp[9,]<-cindexp(t1)
t1<-test[test$LDL<=2.6,]
cindp[10,]<-cindexp(t1)
t1<-test[test$LDL>2.6 & test$LDL<=4.2,]
cindp[11,]<-cindexp(t1)
t1<-test[test$LDL>4.2,]
cindp[12,]<-cindexp(t1)
t1<-test[test$age<=44,]
cindp[13,]<-cindexp(t1)
t1<-test[test$age>44 & test$age<=54,]
cindp[14,]<-cindexp(t1)
t1<-test[test$age>54 & test$age<=64,]
cindp[15,]<-cindexp(t1)
t1<-test[test$age>64,]
cindp[16,]<-cindexp(t1)


#批量进行单因素Cox生存分析
covariates <- c("NRA","sex","age","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
univ_formulas <- sapply(covariates,function(x) as.formula(paste('Surv(time, CC)~', x)))
#univ_formulas
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = test)})
#univ_results
#要整理输出，注意要一个因素一行才行，像treatmentArm这样子chr类型的，会被当成不同种类而有不同的p-value（见下文多因素的输出）
#chr类型会被堪称不同种类，所以要变成numeric（看age）
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2);#coeficient beta
                         HR <-signif(x$coef[2], digits=2);#exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"],2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR, wald.test, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)","wald.test", "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
as.data.frame(res)
#多因素Cox生存分析

coxR<-lapply(2:(dim(test)[2]),function(x){
  #coxR<-lapply(11:19,function(x){
  n=test[,x]
  res.cox <- coxph(Surv(time, CC) ~ qrisk3+sex+age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = test)  # complex Coxph, ori expr
  s <- summary(res.cox)
  gene<-colnames(test)[x]
  
  #也可以单独提取信息，但千万要注意选定好p-value 和HR，不同特征的位置不一样
  p.value<-signif(signif(s$coefficients)[1,5], digits=3)
  HR <-exp(signif(s$coefficients[1], digits=3))
  res<-c(gene,HR,p.value)
  return(res)
})
res <- t(as.data.frame(coxR, check.names = FALSE))
row.names(res)=res[,1]
res=as.data.frame(res[,-1])
res[] <- lapply(res, function(x) {
  if(is.factor(x)) as.numeric(as.character(x)) else x
}) 
colnames(res)=c("HR","pvalue")

res_p0=res[res[, "pvalue"]<0.05,]
res_p0=res_p0[order(res_p0[,"HR"],decreasing=TRUE),]
head(res_p0)



library(forestplot)
data <- read.csv("Desktop/Temp/ForestPlotData.csv", stringsAsFactors=FALSE)
data <- fread("Desktop/Temp/test.txt",sep="\t")
data <- read.csv("Desktop/Temp/Decile.csv", stringsAsFactors=FALSE)
data <- read.csv("Desktop/Temp/Cindex.csv", stringsAsFactors=FALSE)
#查看数据
head(data)

## 构建tabletext，更改列名称，展示更多信息
#np <- ifelse(!is.na(data$Count), paste(data$Count," (",data$Percent,")",sep=""), NA)

## The rest of the columns in the table.
tabletext <- cbind(c("Subgroup","\n",data$Variable),
                   c("Reference model overall C-index (95% CI)","\n",data$C1),
                   c("Plus NRA C-index changes (95% CI)","\n",data$dcNRA),
                   c("Plus PRS C-index changes (95% CI)","\n",data$dcPRS))
#c("P Value","\n",data$P))
##绘制森林图
forestplot(labeltext=tabletext, graph.pos=3,
           mean=c(NA,NA,data$dcNRAm),
           lower=c(NA,NA,data$dcNRAm-data$dcNRAs), upper=c(NA,NA,data$dcNRAm+data$dcNRAs),
           boxsize=0.5)

## 定义亚组，方便后面线条区分
subgps <- c(4,5,8,9,12,13,16,17,20,21,24,25,28,29,32,33)
data$Variable[subgps] <- paste("  ",data$Variable[subgps])

cairo_pdf("Figure3_surv.pdf",width =12.21,height=4.67)
cairo_pdf("Desktop/Temp/Figure3_surv.pdf",width =12.50,height=5.90)
forestplot(labeltext=tabletext,
           graph.pos=3, #为Pvalue箱线图所在的位置
           mean=c(NA,NA,data$dcNRAm),
           lower=c(NA,NA,data$dcNRAm-data$dcNRAs), upper=c(NA,NA,data$dcNRAm+data$dcNRAs),
           #定义标题
           title="C-index Plot",
           ##定义x轴
           xlab="Plus NRA C-index changes (95% CI) vs. reference model",
           ##根据亚组的位置，设置线型，宽度造成“区块感”
           #hrzl_lines=list("3" = gpar(lwd=1, col="#99999922"),
           #                 "7" = gpar(lwd=60, lineend="butt", columns=c(2:6), col="#99999922"),
           #                 "15" = gpar(lwd=60, lineend="butt", columns=c(2:6), col="#99999922"),
           #                 "23" = gpar(lwd=60, lineend="butt", columns=c(2:6), col="#99999922"),
           #                 "31" = gpar(lwd=60, lineend="butt", columns=c(2:6), col="#99999922")),
           #fpTxtGp函数中的cex参数设置各个组件的大小
           txt_gp=fpTxtGp(label=gpar(cex=0.85),#1.25
                          ticks=gpar(cex=1.1),
                          xlab=gpar(cex = 1.2),
                          title=gpar(cex = 1.2)),
           ##fpColors函数设置颜色
           col=fpColors(box="#1c61b6", lines="#1c61b6", zero = "gray50"),
           #箱线图中基准线的位置
           zero=0,#zero=1
           cex=0.9, lineheight = "auto",
           colgap=unit(8,"mm"),
           #箱子大小，线的宽度
           lwd.ci=2, boxsize=0.2)
#箱线图两端添加小竖线，高度
#ci.vertices=TRUE, ci.vertices.height = 0.4)

setwd("/home/pang/Desktop/CAD_IS_PAD_20200121/data/Final")
namematch<-fread("HES/match.txt")
test1<-fread("HES/hesin_diag.txt")
t2<-namematch[match(test1$eid,namematch$malik),]
test1$schunkert<-t2$schunkert
# use only incidence 0 
test1<-test1[test1$ins_index==0,]
test1<-test1[test1$level!=3,]
test1.1<-filter(test1,diag_icd9==4400 | diag_icd9==4402 | diag_icd9==4438 | diag_icd9==4439 | diag_icd10=="I7000" | diag_icd10=="I7001" |diag_icd10=="I702" |diag_icd10=="I7020" |diag_icd10=="I7021" |diag_icd10=="I708" |diag_icd10=="I7080" |diag_icd10=="I709" |diag_icd10=="I7090" |diag_icd10=="I738" |diag_icd10=="I739")
fwrite(test1.1,"survival/ICD910_PAD.txt")
test4<-fread("survival/test1")
test4.1<-filter(test4,f.20002.0.0==1067 | f.20002.0.0==1087 | f.20002.0.0==1088 | f.20004.0.0==1102 | f.20004.0.0==1108 | f.20004.0.0==1440)
fwrite(test4.1,"survival/self-report_PAD.txt")
test2<-fread("HES/hesin_oper.txt")
t2<-namematch[match(test2$eid,namematch$malik),]
test2$schunkert<-t2$schunkert
test2<-test2[test2$ins_index==0,]
test2<-test2[test2$level!=3,]
test2.1<-filter(test2,oper4=="X093" | oper4=="X094" | oper4=="X095" | oper4=="L216" | oper4=="L513" | oper4=="L516" | oper4=="L518" | oper4=="L521" | oper4=="L522" | oper4=="L541" | oper4=="L544" | oper4=="L548" | oper4=="L591" | oper4=="L592" | oper4=="L593" | oper4=="L594" | oper4=="L595" | oper4=="L596" | oper4=="L597" | oper4=="L598" | oper4=="L601" | oper4=="L602" | oper4=="L631" | oper4=="L635" | oper4=="L639" | oper4=="L667")
fwrite(test2.1,"survival/OPCS4_PAD.txt")
test3<-fread("HES/hesin.txt")
t2<-namematch[match(test3$eid,namematch$malik),]
test3$schunkert<-t2$schunkert
test3<-test3[test3$ins_index==0,]
# convert data format
library(lubridate)
test3$epistart<-dmy(test3$epistart)
### sample sum
test5<-fread("survival/PAD_all_rainer.txt")
t2<-test3[match(test5$eid,test3$eid),]
test5$schunkert<-t2$schunkert
test5$epistart<-t2$epistart
fwrite(test5, "survival/PAD_epistart.txt")
### survival data summary
CIPsurv<-Fabian
t2<-assdate[match(CIPsurv$sample,assdate$f.eid),]
CIPsurv$assdate<-t2$f.53.0.0
t2<-test5[match(CIPsurv$sample,test5$schunkert),]
CIPsurv$PADdate<-t2$epistart
CIPsurv$PADdate<-as.IDate(CIPsurv$PADdate)
# test5 2282 PAD cases
# 1415 in namematch, 867 lost
#t0<-subset(test5, !(eid %in% namematch$malik))
# 989 lost in CAD UKB
#t1<- subset(test5, !(schunkert %in% CIPsurv$sample))
# 122 exist in Rainer, not in CADUKB
#t2<-subset(t1,!(eid %in% t0$eid))
t2<-HES[match(CIPsurv$sample,HES$f.eid),]
CIPsurv$MIdate<-t2$f.42000.0.0
CIPsurv$ISdate<-t2$f.42008.0.0
CIPsurv$MIdate[CIPsurv$MIdate=="1900-01-01"]<-NA

#diagnose
CIPsurv$diagnose<-"Control"
CIPsurv$diagnose[is.na(CIPsurv$MIdate)==FALSE]<-"CAD"
CIPsurv$diagnose[is.na(CIPsurv$ISdate)==FALSE]<-"IS"
CIPsurv$diagnose[is.na(CIPsurv$PADdate)==FALSE]<-"PAD"
# state 0control 1CAD 2IS 3PAD
CIPsurv$state<-0
CIPsurv$state[CIPsurv$diagnose=="CAD"]<-1
CIPsurv$state[CIPsurv$diagnose=="IS"]<-2
CIPsurv$state[CIPsurv$diagnose=="PAD"]<-3
# time calculation
CIPsurv$time<-10000
CIPsurv$time0<-CIPsurv$MIdate-CIPsurv$assdate
CIPsurv$time0[CIPsurv$time0<0]<-"remove"
CIPsurv$time[CIPsurv$diagnose=="CAD"]<-CIPsurv$time0[CIPsurv$diagnose=="CAD"]
CIPsurv$time0<-CIPsurv$ISdate-CIPsurv$assdate
CIPsurv$time0[CIPsurv$time0<0]<-"remove"
CIPsurv$time[CIPsurv$diagnose=="IS"]<-CIPsurv$time0[CIPsurv$diagnose=="IS"]
CIPsurv$time0<-CIPsurv$PADdate-CIPsurv$assdate
CIPsurv$time0[CIPsurv$time0<0]<-"remove"
CIPsurv$time[CIPsurv$diagnose=="PAD"]<-CIPsurv$time0[CIPsurv$diagnose=="PAD"]
CIPsurv$status<-1
CIPsurv$status[CIPsurv$state==0]<-0
CIPsurv$deathtime<-CIPsurv$deathdate-CIPsurv$assdate
test<-CIPsurv[CIPsurv$time!="remove",]
test1<-test[as.numeric(test$time)>test$deathtime,]
t1<-subset(CIPsurv, sample %in% test1$sample)
t1$time<-t1$deathtime
t2<-subset(CIPsurv, !(sample %in% test1$sample))
CIPsurv0<-CIPsurv
CIPsurv<-rbind(t1,t2)
#test2<-test%>%mutate(time=replace(time,sample%in%test1$sample,deathtime))

#10 PCs
t2<-caddat2[match(CIPsurv$sample,caddat2$sample),]
CIPsurv$PC1<-t2$PC1
CIPsurv$PC2<-t2$PC2
CIPsurv$PC3<-t2$PC3
CIPsurv$PC4<-t2$PC4
CIPsurv$PC5<-t2$PC5
CIPsurv$PC6<-t2$PC6
CIPsurv$PC7<-t2$PC7
CIPsurv$PC8<-t2$PC8
CIPsurv$PC9<-t2$PC9
CIPsurv$PC10<-t2$PC10

test<-CIPsurv[CIPsurv$time!="remove"]
test$time<-as.numeric(test$time)
test<-test[test$state!=0,]
#test2<-test
#test2$state<-4
#test<-rbind(test,test2)
fit<- survfit(Surv(time, status) ~ state, data = test)
#ggsurvplot(fit, fun="event",palette =c("#E7B800", "#2E9FDF","#83C5BE"),data = test,risk.table = TRUE,pval=TRUE,ggtheme = theme_classic(),censor.shape="|", censor.size = 4,size=1,conf.int.alpha=c(0.1),conf.int=F,ncensor.plot = TRUE)
ggsurvplot(fit, 
           data=test,
           #fun = "cumhaz",
           fun="event",
           risk.table = FALSE, 
           cumevents = TRUE,
           #pval = TRUE, 
           #pval.coord = c(0, 0.25), conf.int = FALSE,
           legend.labs = c("CAD","IS","PAD","CADISPAD"),
           xlab = "Time in days",
           #cumevents.title = "Cumulative number of recurrences",
           size = rep(0.7, 4),
           #xlim = c(0, 100), ylim = c(0, 0.1),
           alpha = 0.7,
           #break.time.by = 10,
           ggtheme = theme_classic(),legend.title="Disease",
           # risk.table.y.text.col = TRUE,
           # risk.table.y.text = TRUE,
           palette = "Set1")

pdf("~/Desktop/CAD_IS_PAD_20200121/data/Final/figures/Survival_CADISPAD_4.pdf")

t0<-test[test$state==0,]
t0<-t0[sample(nrow(t0), nrow(t0)/20), ]
t1<-test[test$state!=0,]
t1<-t1[sample(nrow(t1),nrow(t1)/20),]
t2<-rbind(t0,t1)
library(dynpred)
t1<-t2%>%select("sample","sex","age","metagrs","status","time","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")
t1<-t1[complete.cases(t1)==TRUE,]
cindex(Surv(time, status) ~ metagrs+sex+age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=t1)

test0<-test

# always strange result???
t1<-test[test$qrisk3>20,]
t1<-test[test$qrisk3>5 & test$qrisk3<=20,]
t1<-test[test$qrisk3<=5,]
t1<-t1[is.na(t1$qrisk3)==FALSE,]
res.cox <- coxph(Surv(time2, CC) ~ qrisk3+sex+age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = t1)
summary(res.cox)
res.cox <- coxph(Surv(time2, CC) ~ qrisk3, data = t1)
t1<-test[test$SCORE<=0.05,]
res.cox <- coxph(Surv(time2, CC) ~ SCORE+sex+age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = t1)

#### Trait-specific PRS
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

analysis<-function(dat,pheno){
  
  modComp  <- summary( glm(CC~NRA,family="binomial",data=dat) )$coefficients
  dat$grp  <- round(dat$NRA) #cut(dat$NRA,breaks = seq(min(dat$NRA),max(dat$NRA),by=2),include.lowest=TRUE)
  #dat$grp<-cut(dat$NRA,breaks = seq(min(dat$NRA),max(dat$NRA),by=0.1),include.lowest=TRUE)
  prev_p   <- aggregate(dat$CC,by=list(dat$grp),FUN=mean,na.rm=T)
  prev_n   <- aggregate(dat$CC,by=list(dat$grp),FUN=sum,na.rm=T)
  prev_N   <- aggregate(dat$CC,by=list(dat$grp),FUN=function(x) length(!is.na(x)))
  allele_c <- aggregate(dat$NRA,by=list(dat$grp),FUN=mean)
  d        <- cbind.data.frame(N=prev_N$x,n=prev_n$x,p=prev_p$x,numAllele=allele_c$x)
  d$p_se   <- sqrt(d$p*(1-d$p)/d$N)
  d        <- d[which(d$N>=200),]
  d$normnum<-(d$numAllele-min(dat$NRA))/(max(dat$NRA)-min(dat$NRA))
  
  #modelLogit  <- glm(cbind(n,N)~numAllele,data=d,family = binomial(link="logit"))
  #modelLog    <- glm(cbind(n,N)~numAllele,data=d,family = binomial(link="log"))
  modelLogit  <- glm(p~numAllele,data=d,family = binomial(link="logit"))
  modelLog    <- glm(p~numAllele,data=d,family = binomial(link="log"))
  
  predLogit   <- predict(modelLogit,type = "response")
  predLog     <- predict(modelLog,type = "response")
  
  Rlogit      <- round(unlist(cor.test(d$p,predLogit)[c("estimate","conf.int")]),2)
  Rlog        <- round(unlist(cor.test(d$p,predLog)[c("estimate","conf.int")]),2)
  
  Cols <- c("seagreen","coral1","goldenrod","blue")
  l    <- which(d[,"p_se"]>0)
  
  plot(d[l,"normnum"],d[l,"p"],pch=19,axes=FALSE,ylim=c(0,max(d$p+d$p_se)),xlim=c(0,1),main=pheno,
       xlab="Polygenic risk score",ylab="Prevalence");axis(1);axis(2)
  for(i in 1:nrow(d)){
    if(d[i,"p_se"]>0){
      arrows(d[i,"normnum"],d[i,"p"]-d[i,"p_se"],d[i,"normnum"],d[i,"p"]+d[i,"p_se"],col="grey",angle=90,len=0.05)
      arrows(d[i,"normnum"],d[i,"p"]+d[i,"p_se"],d[i,"normnum"],d[i,"p"]-d[i,"p_se"],col="grey",angle=90,len=0.05)
    }
  }
  matlines(d$normnum,cbind(predLogit,predLog),col=Cols,lwd=2)
  quant<-quantile(dat$NRA,probs=c(0.1,0.9))
  polygon(x=c((quant[1]-min(dat$NRA))/(max(dat$NRA)-min(dat$NRA)),(quant[1]-min(dat$NRA))/(max(dat$NRA)-min(dat$NRA)),(quant[2]-min(dat$NRA))/(max(dat$NRA)-min(dat$NRA)),(quant[2]-min(dat$NRA))/(max(dat$NRA)-min(dat$NRA))),y=c(0,max(d$p),max(d$p),0),col=rgb(0.2,0.2,0.2,0.2),border = NA)
  legend(x=0,y=max(d$p+d$p_se),
         legend=c(paste0("Logit: R=",Rlogit[1]," (95%CI: [",Rlogit[2],"-",Rlogit[3],"])"),
                  paste0("Log: R=",Rlog[1]," (95%CI: [",Rlog[2],"-",Rlog[3],"])")),
         col=Cols,cex=0.9,box.lty=0,lwd=2,bg="transparent")
  #dev.off()
}

analysisde<-function(dat,pheno){
  dat<-dat%>%mutate(quantile=ntile(NRA,10))
  prev_p   <- aggregate(dat$CC,by=list(dat$quantile),FUN=mean,na.rm=T)
  prev_n   <- aggregate(dat$CC,by=list(dat$quantile),FUN=sum,na.rm=T)
  prev_N   <- aggregate(dat$CC,by=list(dat$quantile),FUN=function(x) length(!is.na(x)))
  allele_c <- aggregate(dat$NRA,by=list(dat$quantile),FUN=mean)
  d        <- cbind.data.frame(N=prev_N$x,n=prev_n$x,p=prev_p$x,numAllele=allele_c$x)
  d$p_se   <- sqrt(d$p*(1-d$p)/d$N)
  d        <- d[which(d$N>=400),]
  d$normnum<-(d$numAllele-min(dat$NRA))/(max(dat$NRA)-min(dat$NRA))
  
  #modelLogit  <- glm(cbind(n,N)~numAllele,data=d,family = binomial(link="logit"))
  #modelLog    <- glm(cbind(n,N)~numAllele,data=d,family = binomial(link="log"))
  modelLogit  <- glm(p~numAllele,data=d,family = binomial(link="logit"))
  modelLog    <- glm(p~numAllele,data=d,family = binomial(link="log"))
  
  predLogit   <- predict(modelLogit,type = "response")
  predLog     <- predict(modelLog,type = "response")
  
  Rlogit      <- round(unlist(cor.test(d$p,predLogit)[c("estimate","conf.int")]),2)
  Rlog        <- round(unlist(cor.test(d$p,predLog)[c("estimate","conf.int")]),2)
  
  Cols <- c("seagreen","coral1","goldenrod","blue")
  l    <- which(d[,"p_se"]>0)
  
  plot(d[l,"normnum"],d[l,"p"],pch=19,axes=FALSE,ylim=c(0,max(d$p+d$p_se)),xlim=c(0,1),main=pheno,
       xlab="Polygenic risk score",ylab="Prevalence");axis(1);axis(2)
  for(i in 1:nrow(d)){
    if(d[i,"p_se"]>0){
      arrows(d[i,"normnum"],d[i,"p"]-d[i,"p_se"],d[i,"normnum"],d[i,"p"]+d[i,"p_se"],col="grey",angle=90,len=0.05)
      arrows(d[i,"normnum"],d[i,"p"]+d[i,"p_se"],d[i,"normnum"],d[i,"p"]-d[i,"p_se"],col="grey",angle=90,len=0.05)
    }
  }
  matlines(d$normnum,cbind(predLogit,predLog),col=Cols,lwd=2)
  quant<-quantile(dat$NRA,probs=c(0.1,0.9))
  polygon(x=c((quant[1]-min(dat$NRA))/(max(dat$NRA)-min(dat$NRA)),(quant[1]-min(dat$NRA))/(max(dat$NRA)-min(dat$NRA)),(quant[2]-min(dat$NRA))/(max(dat$NRA)-min(dat$NRA)),(quant[2]-min(dat$NRA))/(max(dat$NRA)-min(dat$NRA))),y=c(0,max(d$p),max(d$p),0),col=rgb(0.2,0.2,0.2,0.2),border = NA)
  legend(x=0,y=max(d$p+d$p_se),
         legend=c(paste0("Logit: R=",Rlogit[1]," (95%CI: [",Rlogit[2],"-",Rlogit[3],"])"),
                  paste0("Log: R=",Rlog[1]," (95%CI: [",Rlog[2],"-",Rlog[3],"])")),
         col=Cols,cex=0.9,box.lty=0,lwd=2,bg="transparent")
  #dev.off()
}

#CAD
cadprs<-fread("~/Desktop/CAD_IS_PAD_20200121/data/Final/PRS/Specific/UKBB_CAD.all.score")
cadprs<-cadprs%>%select("FID","0.00050005")
t2<-CIPsurv[match(cadprs$FID,CIPsurv$sample),]
cadprs$CC<-t2$state
cadprs$time<-t2$time
cadprs$status<-t2$status
cadprs<-cadprs[cadprs$CC!=2 & cadprs$CC!=3,]
cadprs<-cadprs[cadprs$time!="remove",]
colnames(cadprs)[2]<-"NRA"
cadprs$CC<-cadprs$CC/2
cadprs$NRA2<-range01(cadprs$NRA)
pdf("UKBB_IS_specific.pdf",width=4.63,height=5.64)
analysis(cadprs,"UKBB IS specific")
dev.off()

test<-CIPsurv[CIPsurv$time!="remove"]
test$time<-as.numeric(test$time)
t1<-test[test$state!=3&test$state!=2,]
res.cox <- coxph(Surv(time, status) ~ CADsprs0+sex+age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = t1)
summary(res.cox)
#summary(res.cox)$concordance

CIPsurv$PRScom0<-1.2894*CIPsurv$CADsprs+1.0921*CIPsurv$ISsprs+1.1838*CIPsurv$PADsprs
CIPsurv$PRScom<-range01(CIPsurv$PRScom0)
res.cox <- coxph(Surv(time, status) ~ PRScom0+sex+age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = test)
cadprs<-CIPsurv%>%select(sample,PRScom0,status)
colnames(cadprs)<-c("sample","NRA","CC")

pdf("UKBB_CADISPAD_1.pdf",width=4.63,height=5.64)
analysis(cadprs,"UKBB CADISPAD")
dev.off()
pdf("UKBB_CADISPAD_2.pdf",width=4.63,height=5.64)
analysisde(cadprs,"UKBB CADISPAD")
dev.off()

t1<-t1[t1$time!=10000,]
fit <- survfit(Surv(time, status) ~ sex, data = t1)
ggsurvplot(fit, fun="event",palette =c("#E7B800", "#2E9FDF"),data = t1,risk.table = TRUE,pval=TRUE,ggtheme = theme_classic(),censor.shape="|", censor.size = 4,size=1,conf.int.alpha=c(0.1),conf.int=F,ncensor.plot = TRUE)
ggsurvplot(fit, fun="event",palette =c("#E7B800", "#2E9FDF"),data = t1,xlim=c(0,5000),risk.table = TRUE,pval=TRUE,ggtheme = theme_classic(),censor.shape="|", censor.size = 4,size=1,conf.int.alpha=c(0.1),conf.int=F,ncensor.plot = TRUE)
pdf("IS_survival_02.pdf",width=6.80,height=7.00)
dev.off()

t1$time[t1$time==10000]<-5000
res.cox <- coxph(Surv(time, status) ~ CADsprs+sex+age+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10, data = t1)
pdf("CAD_survival_04.pdf",width=6.80,height=7.00)
ggsurvplot(survfit(res.cox),data=t1, palette= "#2E9FDF",
           ggtheme = theme_minimal(),ylim=c(0.97,1))+ggtitle("CAD")+xlab("Time (year)")
dev.off()
