dat2l<-split(pcdat2,pcdat2$worder)
wadres<-prevstatn(dat2l,wadres)
prevstatn<-function(dat2l,adres){
  adres<-data.frame("case"=0,"all"=0,"nra"=0,"pheno"=0,"decile"=0)
  for (i in 1:10){
    t<-dat2l[[i]]
    tn<-t[is.na(t$fbr)==FALSE,]
    bmih<-tn[tn$fbr=="Yes",]
    bmil<-tn[tn$fbr=="No",]
    bmih0<-statsub(bmih,"FYes")
    bmil0<-statsub(bmil,"FNo") #26.32
    adres0<-rbind(bmih0,bmil0)
    adres0$decile<-i
    adres<-rbind(adres,adres0)
  }
  return(adres)
}

dat2l<-split(pcdat2,pcdat2$uorder)
uadres<-prevstatn(dat2l,uadres)

1,1001,1002,1003
4,4001,4002,4003
test<-ethic
test$eth<-NA
test$eth[test$ethic==1]<-"white"
test$eth[test$ethic==1001]<-"white"
test$eth[test$ethic==1002]<-"white"
test$eth[test$ethic==1003]<-"white"
test$eth[test$ethic==4]<-"black"
test$eth[test$ethic==4001]<-"black"
test$eth[test$ethic==4002]<-"black"
test$eth[test$ethic==4003]<-"black"
fwrite(test,"~/Desktop/Fabian/pheno/21000_ethic_bw.ukb")
t2<-test[match(pcdat2$sample,test$f.eid),]
pcdat2$race<-t2$eth

test<-moi2
test$mbr<-"No"
test$mbr[test$value==5]<-"Yes"
t2<-test[match(bcdat2$sample,test$f.eid),]
bcdat2$mbr<-t2$mbr

test<-sibi2
test$sbr<-"No"
test$sbr[test$value==5]<-"Yes"

test<-fai2
test$fbr<-"No"
test$fbr[test$value==13]<-"Yes"
t2<-test[match(pcdat2$sample,test$f.eid),]
pcdat2$fbr<-t2$fbr

d<-read.table("Figure5/CAD_MET_uGRS",header=TRUE,sep="\t")
pheno1<-"21<MET"
pheno2<-"MET≤7.5"
d1<-d[d$group==pheno1,]
d2<-d[d$group==pheno2,]
abnorm0<-c(rep(d1$numAllele[1],d1$n[1]),rep(d1$numAllele[2],d1$n[2]),rep(d1$numAllele[3],d1$n[3]),rep(d1$numAllele[4],d1$n[4]),rep(d1$numAllele[5],d1$n[5]),rep(d1$numAllele[6],d1$n[6]),rep(d1$numAllele[7],d1$n[7]),rep(d1$numAllele[8],d1$n[8]),rep(d1$numAllele[9],d1$n[9]),rep(d1$numAllele[10],d1$n[10]))
norm0<-c(rep(d1$numAllele[1],d1$N[1]-d1$n[1]),rep(d1$numAllele[2],d1$N[2]-d1$n[2]),rep(d1$numAllele[3],d1$N[3]-d1$n[3]),rep(d1$numAllele[4],d1$N[4]-d1$n[4]),rep(d1$numAllele[5],d1$N[5]-d1$n[5]),rep(d1$numAllele[6],d1$N[6]-d1$n[6]),rep(d1$numAllele[7],d1$N[7]-d1$n[7]),rep(d1$numAllele[8],d1$N[8]-d1$n[8]),rep(d1$numAllele[9],d1$N[9]-d1$n[9]),rep(d1$numAllele[10],d1$N[10]-d1$n[10]))
o = outer(abnorm0, norm0, "-")
mean((o>0) + .5*(o==0))

abnorm2<-c(rep(d2$numAllele[1],d2$n[1]),rep(d2$numAllele[2],d2$n[2]),rep(d2$numAllele[3],d2$n[3]),rep(d2$numAllele[4],d2$n[4]),rep(d2$numAllele[5],d2$n[5]),rep(d2$numAllele[6],d2$n[6]),rep(d2$numAllele[7],d2$n[7]),rep(d2$numAllele[8],d2$n[8]),rep(d2$numAllele[9],d2$n[9]),rep(d2$numAllele[10],d2$n[10]))
norm2<-c(rep(d2$numAllele[1],d2$N[1]-d2$n[1]),rep(d2$numAllele[2],d2$N[2]-d2$n[2]),rep(d2$numAllele[3],d2$N[3]-d2$n[3]),rep(d2$numAllele[4],d2$N[4]-d2$n[4]),rep(d2$numAllele[5],d2$N[5]-d2$n[5]),rep(d2$numAllele[6],d2$N[6]-d2$n[6]),rep(d2$numAllele[7],d2$N[7]-d2$n[7]),rep(d2$numAllele[8],d2$N[8]-d2$n[8]),rep(d2$numAllele[9],d2$N[9]-d2$n[9]),rep(d2$numAllele[10],d2$N[10]-d2$n[10]))
o2 = outer(abnorm2, norm2, "-")
mean((o2>0) + .5*(o2==0))

quant<-quantile(caddat$NRA,probs = c(0.2,0.4,0.6,0.8,1))
caddat$quintile<-5
caddat[caddat$NRA<=quant[1],22]<-1
caddat[quant[1]<caddat$NRA & caddat$NRA<=quant[2],22]<-2
caddat[quant[2]<caddat$NRA & caddat$NRA<=quant[3],22]<-3
caddat[quant[3]<caddat$NRA & caddat$NRA<=quant[4],22]<-4
caddat2<-caddat[caddat$quintile==5 | caddat$quintile==1,]
d<-caddat2[is.na(caddat2$sex)==FALSE,]
d1<-d[d$sex==1,]
d2<-d[d$sex==2,]
rest(d1,d2,d3)

d<-caddat2[is.na(caddat2$smoke)==FALSE,]
d1<-d[d$smoke==1,]
d2<-d[d$smoke==0,]

d<-caddat2[is.na(caddat2$Diabetes)==FALSE,]
d1<-d[d$Diabetes==1,]
d2<-d[d$Diabetes==0,]

d<-caddat2[is.na(caddat2$bmi)==FALSE,]
d1<-d[d$bmi<30,]
d2<-d[d$bmi>=30,]

d<-caddat2[is.na(caddat2$MET)==FALSE,]
d1<-d[d$MET<7.5,]
d2<-d[d$MET>=21,]
d3<-d[d$MET>=7.5,]
d3<-d3[d3$MET<21,]
rest(d1,d2,d3)

quant<-quantile(bcdat2$NRA,probs = c(0.2,0.4,0.6,0.8,1))
bcdat2$quintile<-5
bcdat2[bcdat2$NRA<=quant[1],17]<-1
bcdat2[quant[1]<bcdat2$NRA & bcdat2$NRA<=quant[2],17]<-2
bcdat2[quant[2]<bcdat2$NRA & bcdat2$NRA<=quant[3],17]<-3
bcdat2[quant[3]<bcdat2$NRA & bcdat2$NRA<=quant[4],17]<-4
caddat2<-bcdat2[bcdat2$quintile==5 | bcdat2$quintile==1,]
d<-caddat2[is.na(caddat2$bmi)==FALSE,]
d1<-d[d$bmi<30,]
d2<-d[d$bmi>=30,]

d<-caddat2[is.na(caddat2$alchol_freq_class)==FALSE,]
d1<-d[d$alchol_freq_class==0,]
d2<-d[d$alchol_freq_class==1,]
d3<-d[d$alchol_freq_class==2,]

quant<-quantile(pcdat2$NRA,probs = c(0.2,0.4,0.6,0.8,1))
pcdat2$quintile<-5
pcdat2[pcdat2$NRA<=quant[1],14]<-1
pcdat2[quant[1]<pcdat2$NRA & pcdat2$NRA<=quant[2],14]<-2
pcdat2[quant[2]<pcdat2$NRA & pcdat2$NRA<=quant[3],14]<-3
pcdat2[quant[3]<pcdat2$NRA & pcdat2$NRA<=quant[4],14]<-4
caddat2<-pcdat2[pcdat2$quintile==1 | pcdat2$quintile==5,]
d<-caddat2[is.na(caddat2$fbr)==FALSE,]
d1<-d[d$fbr=="No",]
d2<-d[d$fbr=="Yes",]

rest<-function(d1,d2,d3){
  #abnorm0<-d1[d1$CC==1,2]
  #norm0<-d1[d1$CC==0,2]
  #o = outer(abnorm0, norm0, "-")
  #mo1<-mean((o>0) + .5*(o==0))
  #abnorm2<-d2[d2$CC==1,2]
  #norm2<-d2[d2$CC==0,2]
  #o2 = outer(abnorm2, norm2, "-")
  #mo2<-mean((o2>0) + .5*(o2==0))
  #abnorm3<-d3[d3$CC==1,2]
  #norm3<-d3[d3$CC==0,2]
  #o3 = outer(abnorm3, norm3, "-")
  #mo3<-mean((o3>0) + .5*(o3==0))
  fit1 <- glm(CC~NRA, family=binomial(link="log"), data=d1)
  co1<-concordance(fit1)[[1]]
  fit2 <- glm(CC~NRA, family=binomial(link="log"), data=d2)
  co2<-concordance(fit2)[[1]]
  fit3 <- glm(CC~NRA, family=binomial(link="log"), data=d3)
  co3<-concordance(fit3)[[1]]
  q<-data.frame("co1"=co1,"co2"=co2,"co3"=co3)
  #q<-data.frame("mo1"=mo1,"mo2"=mo2,"mo3"=mo3,"co1"=co1,"co2"=co2,"co3"=co3)
  return(q)
}

abnorm0<-d1[d1$CC==1,2]
norm0<-d1[d1$CC==0,2]
o = outer(abnorm0, norm0, "-")
mean((o>0) + .5*(o==0))
abnorm2<-d2[d2$CC==1,2]
norm2<-d2[d2$CC==0,2]
o2 = outer(abnorm2, norm2, "-")
mean((o2>0) + .5*(o2==0))
abnorm3<-d3[d3$CC==1,2]
norm3<-d3[d3$CC==0,2]
o3 = outer(abnorm3, norm3, "-")
mean((o3>0) + .5*(o3==0))



fit1 <- glm(CC~NRA, family=binomial(link="log"), data=d1)
concordance(fit1)
fit2 <- glm(CC~NRA, family=binomial(link="log"), data=d2)
concordance(fit2)
fit3 <- glm(CC~NRA, family=binomial(link="log"), data=d3)
concordance(fit3)

gfit1 <- coxph(Surv(age,CC) ~ NRA, data = d1)
concordance(gfit1)
gfit2 <- coxph(Surv(age,CC) ~ NRA, data = d2)
concordance(gfit2)
gfit3 <- coxph(Surv(age,CC) ~ NRA, data = d3)
concordance(gfit3)