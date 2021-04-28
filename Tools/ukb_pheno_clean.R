setwd("/home/pang/Desktop/Fabian")
library("data.table")
# medication
# creat QC.tab according to fam
med<-fread("43089.t3")
med2<-melt(data=med,id.vars = c("f.eid"))
med3<-na.omit(med2)
fwrite(med3,"20003_medication_QC.ukb",sep="\t")

# 27109 data
t27109<-fread("27109.t3")
# AGE
age<-cbind(t27109$f.eid,t27109$f.21022.0.0)
colnames(age)<-c("f.eid","age")
fwrite(age,"age_QC.ukb",sep="\t")

# function for three batches
merba<-function(data,pheno){
  data<-as.data.frame(data)
  data$Final<-data$V2
  for (i in 1:nrow(data)){
    if (is.na(data$V4[i])==FALSE){
      data$Final[i]=data$V4[i]
    } else if(is.na(data$V3[i])==FALSE){
      data$Final[i]=data$V3[i]
    } else{
      data$Final[i]==data$V2[i]
    }
  }
  data<-data[,-(2:4)]
  colnames(data)<-c("f.eid",pheno)
  return(data)
}

mermean<-function(data,pheno,n){
  data<-as.data.frame(data)
  data$mean<-rowMeans(data[,2:n],na.rm=TRUE)
  data<-data[,-(2:n)]
  colnames(data)<-c("f.eid",pheno)
  return(data)
}
# Ethic keep the latest result if conflict
ethic<-cbind(t27109$f.eid,t27109$f.21000.0.0,t27109$f.21000.1.0,t27109$f.21000.2.0)
ethic<-merba(ethic,"ethic")
# test inconsistent all.equal(a,b)
fwrite(ethic,"21000_ethic_QC.ukb",sep="\t")

# BMI 
bmi<-cbind(t27109$f.eid,t27109$f.21001.0.0,t27109$f.21001.1.0,t27109$f.21001.2.0)
bmi<-merba(bmi,"bmi")
fwrite(bmi,"21001_bmi_QC.ukb",sep="\t")
## Systolic blood pressure 
# use average
bp<-cbind(t27109$f.eid,t27109$f.4080.0.0,t27109$f.4080.0.1,t27109$f.4080.1.0,t27109$f.4080.1.1,t27109$f.4080.2.0,t27109$f.4080.2.1)
bp<-mermean(bp,"BloodPressure",7)
fwrite(bp,"4080_blopre_QC.ukb",sep="\t")

t40052<-fread("40052.t3")
#Cholesterol 
cho<-cbind(t40052$f.eid,t40052$f.30690.0.0,t40052$f.30690.1.0)
cho<-mermean(hdl,"Cholesterol",3)
fwrite(cho,"30690_chol_QC.ukb",sep="\t")

#HDL
hdl<-cbind(t40052$f.eid,t40052$f.30760.0.0,t40052$f.30760.1.0)
hdl<-mermean(hdl,"HDL",3)
fwrite(hdl,"30760_HDL_QC.ukb",sep="\t")

#LDL
ldl<-cbind(t40052$f.eid,t40052$f.30780.0.0,t40052$f.30780.1.0)
ldl<-mermean(ldl,"LDL",3)
fwrite(ldl,"30780_LDL_QC.ukb",sep="\t")

#TC
tc<-cbind(t40052$f.eid,t40052$f.30870.0.0,t40052$f.30870.1.0)
tc<-mermean(tc,"TC",3)
fwrite(tc,"30870_TC_QC.ukb",sep="\t")

t22580<-fread("22580.t3")
#Diabetes
dia<-cbind(t22580$f.eid,t22580$f.2443.0.0,t22580$f.2443.1.0,t22580$f.2443.2.0)
dia<-mermean(dia,"Diabetes",4)
dia2<-na.omit(dia)
dia2$new<-2
dia2[dia2$Diabetes==0,3]<-0
dia2[dia2$Diabetes==1,3]<-1
dia2$n2<-as.character(dia2$new)
dia2<-dia2[,-(2:3)]
colnames(dia2)[2]<-"Diabetes"
table(unlist(dia2$n2))
#0 no 1 yes 2 ever
fwrite(dia2,"2443_diabetes_QC.ukb",sep="\t")

#Smoking
smoke<-cbind(t22580$f.eid,t22580$f.20116.0.0,t22580$f.20116.1.0,t22580$f.20116.2.0)
smoke<-merba(smoke,"smoke")
# convert 2 to 1
# only keep 0 and 1
smoke2<-na.omit(smoke)
smoke2[smoke2$smoke==2,2]<-1
table(unlist(smoke2$smoke))
smoke21<-smoke2[smoke2$smoke==1,]
smoke22<-smoke2[smoke2$smoke==0,]
smoke2<-rbind(smoke21,smoke22)
fwrite(smoke2,"20116_smoking_QC.ukb",sep="\t")

#Ciagrettes mean
cia<-cbind(t22580$f.eid,t22580$f.3456.0.0,t22580$f.3456.1.0,t22580$f.3456.2.0)
cia<-mermean(cia,"Ciagrettes",4)
fwrite(cia,"3456_ciag_QC.ukb",sep="\t")

#Medication S
meds<-fread("meds.t3")
med2<-melt(data=meds,id.vars = c("f.eid"))
med3<-na.omit(med2)
med2<-med3[med3$value>0,]
fwrite(med2,"6153_medS_QC.ukb",sep="\t")

#Father illness
fai<-fread("fai.t3")
fai2<-melt(data=fai,id.vars = c("f.eid"))
fai3<-na.omit(fai2)
fai2<-fai3[fai3$value>0,]
fwrite(fai2,"20107_FaIl_QC.ukb",sep="\t")

#Mother illness
moi<-fread("moi.t3")
moi2<-melt(data=moi,id.vars = c("f.eid"))
moi3<-na.omit(moi2)
moi2<-moi3[moi3$value>0,]
fwrite(moi2,"20110_MoIl_QC.ukb",sep="\t")

#Sibling illness
sibi<-fread("sibi.t3")
sibi2<-melt(data=sibi,id.vars = c("f.eid"))
sibi3<-na.omit(sibi2)
sibi2<-sibi3[sibi3$value>0,]
fwrite(sibi2,"20111_SibIl_QC.ukb",sep="\t")

#lipoprotein A
lpa<-fread("lpa.t3")
lpa<-mermean(lpa,"LPA",3)
fwrite(lpa,"30790_lpa_QC.ukb",sep="\t")

#exercise
exe_freq<-fread("3637_exercise_freq.ukb")
colnames(exe_freq)<-c("V1","V2","V3","V4")
exfrq<-merba(exe_freq,"exercise_frequency")
exfrq2<-na.omit(exfrq)
exfrq3<-exfrq2[exfrq2$exercise_frequency>0,]
exfrq3[exfrq3$exercise_frequency==1,2]<-1/4
exfrq3[exfrq3$exercise_frequency==2,2]<-3/4
exfrq3[exfrq3$exercise_frequency==3,2]<-1
exfrq3[exfrq3$exercise_frequency==4,2]<-3
exfrq3[exfrq3$exercise_frequency==5,2]<-5
exfrq3[exfrq3$exercise_frequency==6,2]<-7
fwrite(exfrq3,"3637_exercise_freq_perweek.ukb",sep="\t")

exe_du<-fread("3647_exercise_duration.ukb")
colnames(exe_du)<-c("V1","V2","V3","V4")
exdu<-merba(exe_du,"exercise_duration")
exdu2<-na.omit(exdu)
exdu3<-exdu2[exdu2$exercise_duration>0,]
exdu3[exdu3$exercise_duration==2,2]<-15
exdu3[exdu3$exercise_duration==3,2]<-30
exdu3[exdu3$exercise_duration==4,2]<-60
exdu3[exdu3$exercise_duration==5,2]<-90
exdu3[exdu3$exercise_duration==6,2]<-120
exdu3[exdu3$exercise_duration==7,2]<-180
exdu3[exdu3$exercise_duration==1,2]<-5
fwrite(exdu3,"3647_exercise_duration_min.ukb",sep="\t")

ex0<-exdu3[match(exfrq3$f.eid,exdu3$f.eid),]
exfrq3$exercise_duration<-ex0$exercise_duration
exmet<-na.omit(exfrq3)
exmet$MET<-exmet$exercise_frequency*exmet$exercise_duration*0.05
fwrite(exmet,"3647_3637_exercise_MET.ukb",sep="\t")

#ICD10
test<-fread("ICD10.sample")
test<-as.data.frame(test)
rownames(test)<-test$V1
test<-test[,-1]
ind <- apply(test, 1, function(x) all(is.na(x)))
test <- test[ !ind, ]
test2<-fread("pheno/ukb_breast_cancer.txt")
test$sample<-rownames(test)
t2<-test2[match(test$sample,test2$ID),]
test$breastcancer<-t2$Female_breast_cancer
test0<-test %>% select(214,215)
test0[is.na(test0$breastcancer)==TRUE,1]<-0
fwrite(test0,"ICD10.sample")

#Alchol assumption
alchol<-fread("1558_Alchol")
alchol[alchol<0]<-NA
alchol_freq<-alchol[,1:4]
colnames(alchol_freq)<-c("V1","V2","V3","V4")
alchol_freq<-merba(alchol_freq,"freq")
alchol_freq$dosage<-rowSums(alchol[,5:19],na.rm=TRUE)
alchol2<-na.omit(alchol_freq)
alchol_freq<-alchol2
alchol_freq$freq_class<-alchol_freq$freq
alchol_freq[alchol_freq$freq==1,4]<-2
alchol_freq[alchol_freq$freq==2,4]<-2
alchol_freq[alchol_freq$freq==3,4]<-1
alchol_freq[alchol_freq$freq==4,4]<-0
alchol_freq[alchol_freq$freq==5,4]<-0
alchol_freq[alchol_freq$freq==6,4]<-0
fwrite(alchol_freq,"1558_Alchol_QC2.ukb")


#ICD code for diagnose
# SLE
SLE<-apply(icd10, 1, function(r) any(r %like% c("M32")))
test$SLE<-0
test$SLE[SLE]<-1

#CKD
CKD<-apply(icd10, 1, function(r) any(r %like% c("N183", "N184","N185")))
test$CKD<-0
test$CKD[CKD]<-1

#MIG
MIG<-apply(icd10, 1, function(r) any(r %like% c("G43", "G440")))
test$MIG<-0
test$MIG[MIG]<-1

#SMI
SMI<-apply(icd10, 1, function(r) any(r %like% c("F20", "F21","F22","F23","F24","F25","F26","F27","F28","F29","F31")))
test$SMI<-0
test$SMI[SMI]<-1

#ED N48.4, F52.2 (-I27.0, I27.2)
ED<-apply(icd10, 1, function(r) any(r %like% c("N484", "G522")))
test$ED<-0
test$ED[ED]<-1

#HIV
HIV<-apply(icd10, 1, function(r) any(r %like% c("B20", "B21","B22","B23","B24")))
test$HIV<-0
test$HIV[HIV]<-1

t2<-test[match(Fabian$sample,test$f.eid),]
t1<-test[is.na(test$bpm)==TRUE & is.na(test$bmp2)==TRUE,]
test$bmpa<-rowSums(cbind(test$bpm,test$bmp2),na.rm = TRUE)
test$bmpa[test$sample %in% t1$sample]<-NA
