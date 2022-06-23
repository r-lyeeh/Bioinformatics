setwd("/home/pang/Desktop/Fabian/pheno")
library("data.table")
library("tidyverse")
##Data input###
# age 21022_age_QC.ukb
# sex 31_sex 1male 2female
# chol 30690_chol_QC.ukb #38.941
# SBP 4080_blopre_QC.ukb
# smoker 20116_smoking_QC.ukb
# hdl 30760_HDL_QC.ukb #38.67
# ldl 30780_HDL_QC.ukb #38.67
# BP treatment 6153_medS_QC.ukb
# atrial fibrillation ICD10code
# atypicalantipsy 20003_medication_QC.ukb
# steroid tablets (corticosteroids) 20003_medication_QC.ukb
# erectile disfunction (impotence2) ICD10
# migraine ICD10
# rheumatoid arthritis ICD10
# chronic kidney disease ICD10
# severe mental illness ICD10
# systemic lupus erythematosis ICD10
# weight
# height
# bmi 21001_bmi_QC.ukb
# ethinicity 21000_ethic_QC.ukb  
# heart attack relative/family 20107_FaIl_QC.ukb, 20110_MoIl_QC.ukb, 20111_SibIl_QC.ukb
# townsend 189_townsend_QC.ukb
# trig 30870_TC_QC.ukb #88.496
# cpd 3456_ciag_QC.ukb  
# SIMDSC10 26427_SIMD.ukb 
# DM/diabetes1 2443_diabetes_QC.ukb
# cvd 1control 2case

#convert mmol to mgdl

### read data
Fabian<-fread("RCV_pheno.ukb")

#1.Score
#variable: age,chol,SBP,smoker,sex
#HDL,
#output:cvdr

#age,sex,chol,hdl,sbp,smoke

#falseifNA <- function(x){
#  ifelse(is.na(x), FALSE, x)
#}


#risk
riskf<-function(pheno,i){
  #Non CVD
    #step1
    if (pheno$sex[i]==0){
      #male
      s0age<-exp(-exp(-26.7)*(pheno$age[i]-20)^5.64)
      s0age10<-exp(-exp(-26.7)*(pheno$age[i]-10)^5.64)
      w<-0.02*(pheno$chol[i]-6)+0.022*(pheno$SBP[i]-120)+0.63*pheno$smoke[i]
      sage<-(s0age)^exp(w)
      sage10<-(s0age10)^exp(w)
      s10age<-sage10/sage
      risk00<-1-s10age
    }
    else if (pheno$sex[i]==1){
      #female
      s0age<-exp(-exp(-31)*(pheno$age[i]-20)^6.62)
      s0age10<-exp(-exp(-31)*(pheno$age[i]-10)^6.62)
      w<-0.02*(pheno$chol[i]-6)+0.022*(pheno$SBP[i]-120)+0.63*pheno$smoke[i]
      sage<-(s0age)^exp(w)
      sage10<-(s0age10)^exp(w)
      s10age<-sage10/sage
      risk00<-1-s10age
    }
    #step2 smoke
    
  #CVD
    #step1
    if (pheno$sex[i]==0){
      #male
      s0age<-exp(-exp(-22.1)*(pheno$age[i]-20)^4.71)
      s0age10<-exp(-exp(-22.1)*(pheno$age[i]-10)^4.71)
      w<-0.24*(pheno$chol[i]-6)+0.018*(pheno$SBP[i]-120)+0.71*pheno$smoke[i]
      sage<-(s0age)^exp(w)
      sage10<-(s0age10)^exp(w)
      s10age<-sage10/sage
      risk20<-1-s10age
    }
    else if (pheno$sex[i]==1){
      #female
      s0age<-exp(-exp(-29.8)*(pheno$age[i]-20)^6.36)
      s0age10<-exp(-exp(-29.8)*(pheno$age[i]-10)^6.36)
      w<-0.24*(pheno$chol[i]-6)+0.018*(pheno$SBP[i]-120)+0.71*pheno$smoke[i]
      sage<-(s0age)^exp(w)
      sage10<-(s0age10)^exp(w)
      s10age<-sage10/sage
      risk20<-1-s10age
    }
  risk10<-100*(risk00+risk20)
  #q<-c(risk00,risk20,risk10)
  return(risk10)
}

temp<-Fabian %>% select("sample","age","sex","sbp","smoking","chol")
colnames(temp)<-c("sample","age","sex","SBP","smoke","chol")
temp$smoke[temp$smoke<2]<-0
temp$smoke[temp$smoke>2]<-1
temp<-temp[temp$age>=40 & temp$age<=65,]
temp<-temp[complete.cases(temp), ]
temp$risk10<-0
for (i in 1:nrow(temp)){
  temp$risk10[i]<-riskf(temp,i)
}
#nrow(temp)=402960



#2. Framingham
#FRS <- 1 - S_0*(10y)^{exp(sum_(i=1)^p (beta_i*X_i) - sum_(i=1)^p (beta_i*mean(X_i)))}
#variable: age, chol,hdl,untreated,treated (BP treatment),smoking,diabetes(DM),sex
# annotate
#as.data.frame(table(pheno$dia))
#current smoker = 1

framscore<-function(pheno,i){
  #men
  if(pheno$sex[i]==0){
    if (falseifNA(pheno$BPtreatment[i]==1)){
      #treated
      eqLm<-3.06117*log(pheno$age[i])+1.12370*log(pheno$chol[i])+(-0.93263)*log(pheno$hdl[i])+1.99881*log(pheno$SBP[i])+0.65451*pheno$smoke[i]+0.57367*pheno$dia[i]
      eqAm<-eqLm-23.9802
      eqBm<-exp(eqAm)
      FRS<- (1 - 0.88936^eqBm)*100
    } else if (falseifNA(pheno$BPtreatment[i]==0)){
      #untreated
      eqLm<-3.06117*log(pheno$age[i])+1.12370*log(pheno$chol[i])+(-0.93263)*log(pheno$hdl[i])+1.93303*log(pheno$SBP[i])+0.65451*pheno$smoke[i]+0.57367*pheno$dia[i]
      eqAm<-eqLm-23.9802
      eqBm<-exp(eqAm)
      FRS<- (1 - 0.88936^eqBm)*100
    }
  } else if (pheno$sex[i]==1){
    #women
    if (falseifNA(pheno$BPtreatment[i]==1)){
      #treated
      eqLm<-2.32888*log(pheno$age[i])+1.20904*log(pheno$chol[i])+(-0.70833)*log(pheno$hdl[i])+2.82263*log(pheno$SBP[i])+0.52873*pheno$smoke[i]+0.69154*pheno$dia[i]
      eqAm<-eqLm-26.1931
      eqBm<-exp(eqAm)
      FRS<- (1 - 0.95012^eqBm)*100
    } else if (falseifNA(pheno$BPtreatment[i]==0)){
      #untreated
      eqLm<-2.32888*log(pheno$age[i])+1.20904*log(pheno$chol[i])+(-0.70833)*log(pheno$hdl[i])+2.76157*log(pheno$SBP[i])+0.52873*pheno$smoke[i]+0.69154*pheno$dia[i]
      eqAm<-eqLm-26.1931
      eqBm<-exp(eqAm)
      FRS<- (1 - 0.95012^eqBm)*100
    }
  }
  return(FRS)
}

#LWomen = 31.764001 x ln(Age) + 22.465206 x ln(Total cholesterol) + (-1.187731) x ln(HDL cholesterol) + 2.552905 x ln(Systolic BP) + 0.420251 x Treated for blood pressure + 13.07543 x Smoker + -5.060998 x ln(Age) x ln(Total cholesterol) + -2.996945 x ln(Age) x Smoker - 146.5933061

#PWomen = 1 - 0.98767^exp(LWomen)



temp<-Fabian %>% select("sample","age","sex","sbp","smoking","chol","HDL","bmpa","dia")
colnames(temp)<-c("sample","age","sex","SBP","smoke","chol","hdl","BPtreatment","dia")
#current smoker
temp$smoke[temp$smoke<=2]<-0
temp$smoke[temp$smoke>2]<-1
temp<-temp[temp$age>=30 & temp$age<=74,]
temp$chol<-temp$chol*38.67
temp$hdl<-temp$hdl*38.67
temp<-temp[complete.cases(temp), ]
temp$FRS<-0

for (i in 1:nrow(temp)){
  temp$FRS[i]<-framscore(temp,i)
}


#3.Pooled ACC/AHA
#PCE = 1 - S_0(10y)^{exp[sum_(i=1)^p (beta_i*X_i) - sum_(i=1)^p (beta_i*mean(X_i))]}
#variable: age, chol,hdl,untreated,treated,smoking,diabetes,sex,race

poolscore<-function(pheno,i){
  if(pheno$sex[i]==0){
    #male
    if (falseifNA(pheno$BPtreatment[i]==1)){
      #treated
      if (pheno$ethic[i]==1){
        eqLw<-12.344*log(pheno$age[i])+11.853*log(pheno$chol[i])+(-2.664)*log(pheno$age[i])*log(pheno$chol[i])+(-7.990)*log(pheno$hdl[i])+1.769*log(pheno$age[i])*log(pheno$hdl[i])+1.797*log(pheno$SBP[i])+7.837*pheno$smoke[i]+(-1.795)*log(pheno$age[i])*pheno$smoke[i]+0.658*pheno$dia[i]
        eqAw<-eqLw-61.18
        eqBw<-exp(eqAw)
        PCEw<-(1-0.9144^eqBw)*100
      } else if (pheno$ethic[i]==6){
        eqLw<-2.469*log(pheno$age[i])+0.302*log(pheno$chol[i])+(-0.307)*log(pheno$hdl[i])+1.809*log(pheno$SBP[i])+0.549*pheno$smoke[i]+0.645*pheno$dia[i]
        eqAw<-eqLw-19.54
        eqBw<-exp(eqAw)
        PCEw<-(1-0.8954^eqBw)*100
      }
    } else if (falseifNA(pheno$BPtreatment[i]==0)){
      #untreated
      if (pheno$ethic[i]==1){
        eqLw<-12.344*log(pheno$age[i])+11.853*log(pheno$chol[i])+(-2.664)*log(pheno$age[i])*log(pheno$chol[i])+(-7.990)*log(pheno$hdl[i])+1.769*log(pheno$age[i])*log(pheno$hdl[i])+1.764*log(pheno$SBP[i])+7.837*pheno$smoke[i]+(-1.795)*log(pheno$age[i])*pheno$smoke[i]+0.658*pheno$dia[i]
        eqAw<-eqLw-61.18
        eqBw<-exp(eqAw)
        PCEw<-(1-0.9144^eqBw)*100
      } else if (pheno$ethic[i]==6){
        eqLw<-2.469*log(pheno$age[i])+0.302*log(pheno$chol[i])+(-0.307)*log(pheno$hdl[i])+1.916*log(pheno$SBP[i])+0.549*pheno$smoke[i]+0.645*pheno$dia[i]
        eqAw<-eqLw-19.54
        eqBw<-exp(eqAw)
        PCEw<-(1-0.8954^eqBw)*100
      }
    }
  } else if (pheno$sex[i]==1){
    #female
    if (falseifNA(pheno$BPtreatment[i]==1)){
      #treated
      if (pheno$ethic[i]==1){
        eqLw<-(-29.799)*log(pheno$age[i])+4.884*(log(pheno$age[i])^2)+13.540*log(pheno$chol[i])+(-3.114)*log(pheno$age[i])*log(pheno$chol[i])+(-13.578)*log(pheno$hdl[i])+3.149*log(pheno$age[i])*log(pheno$hdl[i])+2.019*log(pheno$SBP[i])+7.574*pheno$smoke[i]+(-1.665)*log(pheno$age[i])*pheno$smoke[i]+0.661*pheno$dia[i]
        eqAw<-eqLw--29.18
        eqBw<-exp(eqAw)
        PCEw<-(1-0.9665^eqBw)*100
      } else if (pheno$ethic[i]==6){
        eqLw<-17.114*log(pheno$age[i])+0.940*log(pheno$chol[i])+(-18.920)*log(pheno$hdl[i])+4.475*log(pheno$age[i])*log(pheno$hdl[i])+29.291*log(pheno$SBP[i])+(-6.432)*log(pheno$age[i])*log(pheno$SBP[i])+0.691*pheno$smoke[i]+0.874*pheno$dia[i]
        eqAw<-eqLw-86.61
        eqBw<-exp(eqAw)
        PCEw<-(1-0.9533^eqBw)*100
      }
    } else if (falseifNA(pheno$BPtreatment[i]==0)){
      #untreated
      if (pheno$ethic[i]==1){
        eqLw<-(-29.799)*log(pheno$age[i])+4.884*(log(pheno$age[i])^2)+13.540*log(pheno$chol[i])+(-3.114)*log(pheno$age[i])*log(pheno$chol[i])+(-13.578)*log(pheno$hdl[i])+3.149*log(pheno$age[i])*log(pheno$hdl[i])+1.957*log(pheno$SBP[i])+7.574*pheno$smoke[i]+(-1.665)*log(pheno$age[i])*pheno$smoke[i]+0.661*pheno$dia[i]
        eqAw<-eqLw--29.18
        eqBw<-exp(eqAw)
        PCEw<-(1-0.9665^eqBw)*100
      } else if (pheno$ethic[i]==6){
        eqLw<-17.114*log(pheno$age[i])+0.940*log(pheno$chol[i])+(-18.920)*log(pheno$hdl[i])+4.475*log(pheno$age[i])*log(pheno$hdl[i])+27.820*log(pheno$SBP[i])+(-6.087)*log(pheno$age[i])*log(pheno$SBP[i])+0.691*pheno$smoke[i]+0.874*pheno$dia[i]
        eqAw<-eqLw-86.61
        eqBw<-exp(eqAw)
        PCEw<-(1-0.9533^eqBw)*100
      }
    }
  }
  return(PCEw)
}
#temp$chol<-temp$chol*38.941

temp<-Fabian %>% select("sample","age","sex","sbp","smoking","chol","HDL","bpm","dia","ethic")
colnames(temp)<-c("sample","age","sex","SBP","smoke","chol","hdl","BPtreatment","dia","ethic")
#current smoker
temp$smoke[temp$smoke<=2]<-0
temp$smoke[temp$smoke>2]<-1
temp<-temp[temp$age>=40 & temp$age<=79,]
temp$chol<-temp$chol*38.67
temp$hdl<-temp$hdl*38.67
temp<-temp[complete.cases(temp), ]
temp$PCE<-0
for (i in 1:nrow(temp)){
  temp$PCE[i]<-poolscore(temp,i)
}

pcr_probs <- with(
  subset(Fabian,age>=40& chol>3.37 & chol<8.27 &sbp>=90 &sbp<=200),
  predict_10yr_ascvd_risk(
    sex = sex,
    sex_levels = list(female = 1, male = 0),
    race = ethic,
    age_years = age,
    chol_total_mgdl = chol*38.67, ##*38.67
    chol_hdl_mgdl = HDL*38.67,
    bp_sys_mmhg = sbp,
    bp_meds = bpm,
    smoke_current = smoking,
    diabetes = dia,
    race_levels = list(black = 6, white = 1),
    smoke_current_levels = list(no =c(1,2), yes = 4),
    bp_meds_levels = list(no = 0, yes = 1),
    diabetes_levels = list(no = 0, yes =1)
  )
)
temp<-subset(Fabian,age>=40& chol>3.37 & chol<8.27 &sbp>=90 &sbp<=200)
temp$PCE<-pcr_probs
##1 female 0 male
##1 male 2 female

#4.QRISK3
#variable: age, sex, race, BMI, chol, hdl, SBP, BP treatment, DM, smoke, CPD, family, Townsend deprivation index, rheumatoid arthritis, atrial fibrillation, SBP variability, CKD stage, SLE, erectile dysfunction, migraine, severe mental illness, HIV/AIDS, Corticosteroids, Atyp. antypsychotics
#0men 1women kg cm 1 White or not stated
#2 Indian
#3 Pakistani
#4 Bangladeshi
#5 Other Asian
#6 Black Caribbean
#7 Black African
#8 Chinese
#9 Other ethnic group
#blood pressure mmHG
#1 non-smoker
#2 ex-smoker
#3 light smoker (less than 10)
#4 moderate smoker (10 to 19)
#5 heavy smoker (20 or over)
library(QRISK3)
test<-Fabian %>% select("age","sex","AF","atypical","corticosteroids","ED","MIG","RA","CKD","SMI","SLE","bmpa","dia1","dia2","weight","height","ethic","familyHeart","chol_hdl","sbp","stdsbp","smoking","Townsend","BMI","sample")
test<-test[complete.cases(test), ]
test<-test[test$age>=25 & test$age<=84,]
#test<-test[test$chol_hdl<=11,]
test_all_rst <- QRISK3_2017(data=test_all, patid="sample", gender="sex", age="age",atrial_fibrillation="AF", atypical_antipsy="atypical",regular_steroid_tablets="corticosteroids", erectile_disfunction="ED",migraine="MIG", rheumatoid_arthritis="RA",chronic_kidney_disease="CKD", severe_mental_illness="SMI",systemic_lupus_erythematosis="SLE",blood_pressure_treatment="bpm", diabetes1="dia1",diabetes2="dia2", weight="weight", height="height",ethiniciy="ethic", heart_attack_relative="familyHeart",cholesterol_HDL_ratio="chol_hdl", systolic_blood_pressure="sbp",std_systolic_blood_pressure="stdsbp", smoke="smoking",townsend="Townsend")

#bloodtreatment medicine ->new
#tset01<-subset(test0, !(ID %in% test$sample))

test_all <- QRISK3_2019_test
test_all<-test_all[,-24]
test_all<-test_all[,-1]
colnames(test)<-colnames(test_all)
test_all<-rbind(test_all,test)
test_all_rst <- QRISK3_2017(data=test_all, patid="ID", gender="gender", age="age",atrial_fibrillation="b_AF", atypical_antipsy="b_atypicalantipsy",regular_steroid_tablets="b_corticosteroids", erectile_disfunction="b_impotence2",migraine="b_migraine", rheumatoid_arthritis="b_ra",chronic_kidney_disease="b_renal", severe_mental_illness="b_semi",systemic_lupus_erythematosis="b_sle",blood_pressure_treatment="b_treatedhyp", diabetes1="b_type1",diabetes2="b_type2", weight="weight", height="height",ethiniciy="ethrisk", heart_attack_relative="fh_cvd",cholesterol_HDL_ratio="rati", systolic_blood_pressure="sbp",std_systolic_blood_pressure="sbps5", smoke="smoke_cat",townsend="town")

### strange need to merge my data with example data 
#test_all <- QRISK3_2019_test
#test_all<-rbind(test_all,test)



#5.ASSIGN
#Men
#variable:age,tc,hdlc,sbp,diabetes,family,cpd,SIMDSC10,sex
assscore<-function(pheno,i){
  if (pheno$sex[i]==0){
    #male
    L<-0.05698*pheno$age[i]+0.22286*pheno$chol[i]+(-0.53684)*pheno$hdl[i]+0.0118*pheno$SBP[i]+0.81558*pheno$dia[i]+0.27500*pheno$fam[i]+0.02005*pheno$ciag[i]+0.006296*pheno$simd[i]
    A<-L-5.447113
    B<-exp(A)
    assP<-100*(1-(0.8831^B))
  } else if (pheno$sex[i]==1){
    L<-0.07203*pheno$age[i]+0.12720*pheno$chol[i]+(-0.55836)*pheno$hdl[i]+0.01064*pheno$SBP[i]+0.97727*pheno$dia[i]+0.49159*pheno$fam[i]+0.02724*pheno$ciag[i]+0.009386*pheno$simd[i]
    A<-L-5.418379
    B<-exp(A)
    assP<-100*(1-(0.9365^B))
  }
  return(assP)
}


temp<-Fabian %>% select("sample","age","sex","sbp","chol","HDL","SIMD","ciag","familyHeart","dia")
colnames(temp)<-c("sample","age","sex","SBP","chol","hdl","simd","ciag","fam","dia")
temp<-temp[temp$age>=30 & temp$age<=74,]
#SMID=15.89
#temp$simd[is.na(temp$simd)==TRUE]<-15.89
temp<-temp[complete.cases(temp), ]
temp$assP<-0

for (i in 1:nrow(temp)){
  temp$assP[i]<-assscore(temp,i)
}

#6.Procam
#stoke
#sex male:1 female:0;DM>120=1;SBP>140
#variable: age, sex,hdl,LDL, TC, SBP, DM, smoke, family
sf<-function(x,start,end,step,endr){
  if(x<=start){
    return(0)
  }else if (x>start &x<=end){
    return(ceiling((x-start)/step))
  }else {
    return(endr)
  }
}

riskIS<-function(x){
  if(x<=42){
    return(0.1)
  } else if (x>42&x<=46){
    return(0.2)
  } else if (x>46&x<=49){
    return(0.3)
  } else if (x>49&x<=51){
    return(0.4)
  } else if (x>51&x<=53){
    return(0.5)
  } else if (x>53&x<=55){
    return(0.6)
  } else if (x>55&x<=60){
    return((x-49)/10)
  } else if (x>60&x<=62){
    return((x-48)/10)
  } else if (x>62&x<=64){
    return((x-47)/10)
  } else if (x==65){
    return(1.9)
  } else if (x==66){
    return(2.2)
  } else if (x==67){
    return(2.4)
  } else if (x==68){
    return(2.7)
  } else if (x==69){
    return(3.0)
  } else if (x==70){
    return(3.4)
  } else if (x==71){
    return(3.8)
  } else if (x==72){
    return(4.2)
  } else if (x==73){
    return(4.7)
  } else if (x==74){
    return(5.2)
  } else if (x==75){
    return(5.8)
  } else if (x==76){
    return(6.5)
  } else if (x==77){
    return(7.3)
  } else if (x==78){
    return(8.1)
  } else if (x==79){
    return(9.0)
  } else if (x==80){
    return(10.0)
  } else {
    return("10+")
  }
}

#MI
#smoker 1:0 family 1:0 sex:1 11/9

riskMI<-function(x,age,mi){
    t1<-mim[mim$Age==age,]
    t2<-t1[t1$start<=x & t1$end>=x,]
    risk<-t2$risk
    if(nrow(t2)==0){
      risk<-"30+"
    }
  return(risk)
}

temp<-Fabian %>% select("sample","age","sex","sbp","HDL","LDL","smoking","dia","TC","familyHeart")
colnames(temp)<-c("sample","age","sex","SBP","hdl","ldl","smoke","dia","trig","fam")
temp$smoke[temp$smoke<=2]<-0
temp$smoke[temp$smoke>2]<-1
temp$ldl<-temp$ldl*38.67
temp$hdl<-temp$hdl*38.67
temp$trig<-temp$trig*88.496
temp<-temp[temp$age>=35 & temp$age<=65,]
temp<-temp[complete.cases(temp), ]
temp$procamisscore<-0

for (i in 1:nrow(temp)){
  #scoreIS<-6*(2-temp$sex[i])+1*ceiling(temp$age[i])+9*temp$smoke[i]+7*temp$dia[i]+sf(temp$SBP[i],140,180,5,9)
  scoreIS<-6*temp$sex[i]+1*ceiling(temp$age[i])+9*temp$smoke[i]+7*temp$dia[i]+sf(temp$SBP[i],140,180,5,9)
  temp$procamisscore[i]<-riskIS(scoreIS)
}

temp$procammiscore<-0
mim<-fread("Procam_men")
miwm<-fread("Procam_woman")
for (i in 1:nrow(temp)){
  if (temp$sex[i]==0){
    #male
    scoreMI<-temp$fam[i]*5+temp$smoke[i]*12+temp$dia[i]*9+sf(temp$SBP[i],110,180,10,8)+sf(temp$trig[i],100,200,50,4)+(11-sf(temp$hdl[i],35,55,2,11))+sf(temp$ldl[i],100,196,5,20)
    temp$procammiscore[i]<-riskMI(scoreMI,temp$age[i],mim)
  } else if (temp$sex[i]==1){
    #female
    scoreMI<-temp$fam[i]*5+temp$smoke[i]*12+temp$dia[i]*11+sf(temp$SBP[i],110,180,10,8)+sf(temp$trig[i],100,200,50,4)+(11-sf(temp$hdl[i],35,55,2,11))+sf(temp$ldl[i],100,196,5,20)
    temp$procammiscore[i]<-riskMI(scoreMI,temp$age[i],miwm)
  }
}



#7.Cuore
#1 - [S(t)] ^ {EXP [β1 × age + β2 × SBP + β3 × TC +β4 ×HDL+β5 (if smoker)+β6 (if DM) + β7 (if BP treatment) - G(μ)]} 
#variable: age, sex, chol,HDL,SBP, BP treatment, DM, smoke
cuoscore<-function(pheno,i){
  if (pheno$sex[i]==0){
    Lc<-1-0.953^(exp(0.076*pheno$age[i]+0.013*pheno$SBP[i]+0.006*pheno$chol[i]+(-0.013)*pheno$hdl[i]+0.508*pheno$smoke[i]+0.462*pheno$dia[i]+0.490*pheno$BPtreatment[i]-6.583))
  } else if (pheno$sex[i]==1){
    Lc<-1-0.989^(exp(0.079*pheno$age[i]+0.016*pheno$SBP[i]+0.003*pheno$chol[i]+(-0.015)*pheno$hdl[i]+0.773*pheno$smoke[i]+0.339*pheno$dia[i]+0.590*pheno$BPtreatment[i]-6.016))
  }
  return(Lc)
}

temp<-Fabian %>% select("sample","age","sex","sbp","chol","HDL","smoking","dia","bpm")
colnames(temp)<-c("sample","age","sex","SBP","chol","hdl","smoke","dia","BPtreatment")
temp$smoke[temp$smoke<=2]<-0
temp$smoke[temp$smoke>2]<-1
temp$chol<-temp$chol*38.67
temp$hdl<-temp$hdl*38.67
temp<-temp[temp$age>=25 & temp$age<=84,]
temp<-temp[complete.cases(temp), ]
temp$CuoS<-0

for (i in 1:nrow(temp)){
  temp$CuoS[i]<-cuoscore(temp,i)
}
#unsure

#summary
d<-Fabian %>%select(sample,SCORE)
colnames(d)[2]<-"CVR"
d$group<-"SCORE"
d1<-Fabian %>% select(sample,FRS) 
colnames(d1)[2]<-"CVR"
d1$group<-"FRS"
d<-rbind(d,d1)
d1<-Fabian %>% select(sample,PCE) 
colnames(d1)[2]<-"CVR"
d1$group<-"PCE"
d<-rbind(d,d1)
d1<-Fabian %>% select(sample,ASSIGN) 
colnames(d1)[2]<-"CVR"
d1$group<-"ASSIGN"
d<-rbind(d,d1)
d1<-Fabian %>% select(sample,ASSIGN_SIMD) 
colnames(d1)[2]<-"CVR"
d1$group<-"ASSIGN_SIMD"
d<-rbind(d,d1)
d1<-Fabian %>% select(sample,Cuore) 
colnames(d1)[2]<-"CVR"
d1$group<-"Cuore"
d<-rbind(d,d1)
d1<-Fabian %>% select(sample,qrisk3) 
colnames(d1)[2]<-"CVR"
d1$group<-"QRISK3"
d<-rbind(d,d1)
ggplot(data=d)+geom_histogram(aes(x=CVR,group=group,fill=group),binwidth=1,alpha=0.5)+theme_classic()+xlab("CVD risk score")+ylab("Count")

# NA value fill 
# average risk for specific age;
test<-Fabian
test<-test[is.na(test$age)==FALSE,]
# fill in Na values
meansc<-aggregate(test$SCORE,by=list(test$age),FUN=mean,na.rm=T)
for(i in 1:nrow(test)){
  if(is.na(test$SCORE[i])==TRUE){
    #print(i)
    test$SCORE[i]<-filter(meansc,Group.1==test$age[i])$x
  }
}

meansc<-aggregate(test$FRS,by=list(test$age),FUN=mean,na.rm=T)
for(i in 1:nrow(test)){
  if(is.na(test$FRS[i])==TRUE){
    #print(i)
    test$FRS[i]<-filter(meansc,Group.1==test$age[i])$x
  }
}

meansc<-aggregate(test$PCE,by=list(test$age),FUN=mean,na.rm=T)
for(i in 1:nrow(test)){
  if(is.na(test$PCE[i])==TRUE){
    #print(i)
    test$PCE[i]<-filter(meansc,Group.1==test$age[i])$x
  }
}

meansc<-aggregate(test$Cuore,by=list(test$age),FUN=mean,na.rm=T)
for(i in 1:nrow(test)){
  if(is.na(test$Cuore[i])==TRUE){
    #print(i)
    test$Cuore[i]<-filter(meansc,Group.1==test$age[i])$x
  }
}
