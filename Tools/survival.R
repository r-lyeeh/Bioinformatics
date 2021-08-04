###survival####
##Target: Rsearch on known molecule biomarkers(ER,HER2), pathological stage, sample type and lymph code status.
##Data Source:TCGA_BRCA
breast<-read.csv("/Users/magica/Desktop/clinical_data",sep="\t")
library(survival)
#breast1<-breast[!is.na(breast$X_EVENT),]#remove event na 
breast2<-data.frame(breast$sampleID,breast$X_EVENT,breast$X_OS,breast$AJCC_Stage_nature2012,breast$ER_Status_nature2012,breast$HER2_Final_Status_nature2012,breast$Node_Coded_nature2012,breast$er_level_cell_percentage_category,breast$sample_type)# get needed features and OS/EVENT
#1. Descripitive Statistical 
summary(breast2$breast.AJCC_Stage_nature2012)
summary(breast2$breast.ER_Status_nature2012)
summary(breast2$breast.HER2_Final_Status_nature2012)
summary(breast2$breast.Node_Coded_nature2012)
summary(breast2$breast.er_level_cell_percentage_category)
#2. Data Clean & univariate coxph
breast2<-breast2[!is.na(breast2$breast.X_EVENT),]#remove event na 
breast3<-data.frame(breast$X_EVENT,breast$X_OS,breast$Age_at_Initial_Pathologic_Diagnosis_nature2012)
kmfit1<-survfit(Surv(breast3$breast2.breast.X_OS,breast3$breast2.breast.X_EVENT)~breast3$breast2.breast.sample_type)
plot(kmfit1,col=c("red","blue"),lty=c(2,1))
legend("topright",legend=c("normal","tumor"),col=c("red","blue"),lty=c(2,1),title="Sample_type")
#HER2
breast4<-data.frame(breast2$breast.X_EVENT,breast2$breast.X_OS,breast2$breast.HER2_Final_Status_nature2012)
breast4<-breast4[breast4$breast2.breast.HER2_Final_Status_nature2012!="Equivocal",]
breast4<-breast4[breast4$breast2.breast.HER2_Final_Status_nature2012!="",]
kmfit2<-survfit(Surv(breast4$breast2.breast.X_OS,breast4$breast2.breast.X_EVENT)~breast4$breast2.breast.HER2_Final_Status_nature2012)
summary(kmfit2)
plot(kmfit2,col=c("red","blue"),lty=c(2,1))
legend("topright",legend=c("Negative","Positive"),col=c("red","blue"),lty=c(2,1),title="HER2")

breast5<-data.frame(breast2$breast.X_EVENT,breast2$breast.X_OS,breast2$breast.ER_Status_nature2012)
breast5<-breast5[breast5$breast2.breast.ER_Status_nature2012!="Indeterminate",]
breast5<-breast5[breast5$breast2.breast.ER_Status_nature2012!="",]
kmfit3<-survfit(Surv(breast5$breast2.breast.X_OS,breast5$breast2.breast.X_EVENT)~breast5$breast2.breast.ER_Status_nature2012)
print(kmfit3,show.rmean=True)
plot(kmfit3,col=c("red","blue"),lty=c(2,1))
legend("topright",legend=c("Negative","Positive"),col=c("red","blue"),lty=c(2,1),title="ER")

breast6<-data.frame(breast2$breast.X_EVENT,breast2$breast.X_OS,breast2$breast.Node_Coded_nature2012)
breast6<-breast6[breast6$breast2.breast.Node_Coded_nature2012!="",]
kmfit4<-survfit(Surv(breast6$breast2.breast.X_OS,breast6$breast2.breast.X_EVENT)~breast6$breast2.breast.Node_Coded_nature2012)
print(kmfit4,show.rmean=True)
plot(kmfit4,col=c("red","blue"),lty=c(2,1))
legend("topright",legend=c("Negative","Positive"),col=c("red","blue"),lty=c(2,1),title="Node_Coded")

breast7<-data.frame(breast2$breast.X_EVENT,breast2$breast.X_OS,breast2$breast.er_level_cell_percentage_category)
breast7<-breast7[breast7$breast2.breast.er_level_cell_percentage_category!="",]
new2<-as.character(breast7[,3])
new2[new2!="90-99%"]<-"under 90%"
kmfit5<-survfit(Surv(breast7$breast2.breast.X_OS,breast7$breast2.breast.X_EVENT)~new2)
print(kmfit5,show.rmean=True)
plot(kmfit5,col=c("red","blue"),lty=c(2,1))
legend("topright",legend=c("90-99%","under 90%"),col=c("red","blue"),lty=c(2,1),title="ER_LEVEL")

breast8<-data.frame(breast2$breast.X_EVENT,breast2$breast.X_OS,breast2$breast.AJCC_Stage_nature2012)
breast8<-breast8[breast8$breast2.breast.AJCC_Stage_nature2012!="",]
breast8<-breast8[breast8$breast2.breast.AJCC_Stage_nature2012!="Stage X",]
Stage<-as.character(breast8[,3])
Stage[Stage=="Stage IA"]<-"Stage I"
Stage[Stage=="Stage IB"]<-"Stage I"
Stage[Stage=="Stage IIA"]<-"Stage II"
Stage[Stage=="Stage IIB"]<-"Stage II"
Stage[Stage=="Stage IIIA"]<-"Stage III"
Stage[Stage=="Stage IIIB"]<-"Stage III"
Stage[Stage=="Stage IIIC"]<-"Stage III"
kmfit6<-survfit(Surv(breast8$breast2.breast.X_OS,breast8$breast2.breast.X_EVENT)~Stage)
print(kmfit6,show.rmean=True)
plot(kmfit6,col=c("red","blue","green","purple"))
legend("topright",legend=c("Stage I","Stage II","Stage III","Stage IV"),col=c("red","blue","green","purple"),lty=c(1),title="Stage")

breast9<-breast2
breast9<-breast9[breast9$breast.HER2_Final_Status_nature2012!="Equivocal",]
breast9<-breast9[breast9$breast.HER2_Final_Status_nature2012!="",]
breast9<-breast9[breast9$ER!="Indeterminate",]#
breast9<-breast9[breast9$breast.ER_Status_nature2012!="",]
breast9<-breast9[breast9$breast.er_level_cell_percentage_category!="",]
breast9<-breast9[breast9$breast.AJCC_Stage_nature2012!="",]
breast9<-breast9[breast9$breast.AJCC_Stage_nature2012!="Stage X",]
breast9<-breast9[breast9$breast.Node_Coded_nature2012!="",]
Stage0<-as.character(breast9$breast.AJCC_Stage_nature2012)
Stage0[Stage0=="Stage IA"]<-"Stage I"
Stage0[Stage0=="Stage IB"]<-"Stage I"
Stage0[Stage0=="Stage IIA"]<-"Stage II"
Stage0[Stage0=="Stage IIB"]<-"Stage II"
Stage0[Stage0=="Stage IIIA"]<-"Stage III"
Stage0[Stage0=="Stage IIIB"]<-"Stage III"
Stage0[Stage0=="Stage IIIC"]<-"Stage III"
ER_LEVEL0<-as.character(breast9$ER_LEVEL)
ER_LEVEL0[ER_LEVEL0!="90-99%"]<-"under 90%"
ER<-breast9$ER
HER2<-breast9$HER2
Node_Coded<-breast9$Node_Coded
OS<-breast9$OS
EVENT<-breast9$EVENT
Stage<-Stage0
ER_LEVEL<-ER_LEVEL0

survdiff(Surv(breast9$breast.X_OS,breast9$breast.X_EVENT)~breast9$breast.ER_Status_nature2012+strata(ER_LEVEL))

coxfit=coxph(Surv(OS,EVENT)~Node_Coded+ER+HER2+Stage0+ER_LEVEL0,data=breast9)
coxfit=coxph(Surv(breast.X_OS,breast.X_EVENT)~Stage0,data=breast9)
