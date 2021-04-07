setwd("~/Desktop/ISAR_Moritz_20200910/phenotypes/")
library(readxl)
library(lubridate)
library(nycflights13)
library(tidyverse)
library(reshape2)
library(data.table)

####FUNCTION DEFINITION##########
# count time duration from the baseline time 
timecount<-function(mpv,orig){
# find the baseline data
t0<-aggregate(mpv,by=list(mpv$KEY_PATIENT),FUN=function(x) min(x))
mpv$time<-dseconds(0)
for (i in 1:nrow(t0)){
  n<-which(mpv$KEY_PATIENT==t0[i,1])
  for (j in 1:length(n)){
    mpv[n[j],ncol(mpv)]<-mpv[n[j],orig]-mpv[n[1],orig]
  }
}
return(mpv)
}

convertch<-function(mpv,val0,orig){
  #define time point tag
  q<-ncol(mpv)
mpv$tag<-0
mpv[mpv$time>dseconds(0) & mpv$time<dseconds(21600),(q+1)]<-1 #6h
mpv[mpv$time>=dseconds(21600)&mpv$time<dseconds(86400),(q+1)]<-2 #6h-24h
mpv[mpv$time>=dseconds(86400)&mpv$time<dseconds(172800),(q+1)]<-3 #24h-48h
mpv[mpv$time>=dseconds(172800),(q+1)]<-4 #>48h
# remove beyond range rows
mpv0<-mpv[mpv$tag!=4,]
#NR,KEY,val,crd,time,tag
#mpv0<-mpv0 %>% select(1,2,4,6,10,11)
mpv0<-mpv0 %>% select(1,2,val0,orig,q,(q+1))
# for each time point keep the largest one
mpv1<-mpv0[mpv0$tag==1,]
t1<-aggregate(mpv1,by=list(mpv1$KEY_PATIENT),FUN=function(x) max(x))
mpv2<-mpv0[mpv0$tag==2,]
t2<-aggregate(mpv2,by=list(mpv2$KEY_PATIENT),FUN=function(x) max(x))
mpv3<-mpv0[mpv0$tag==3,]
t3<-aggregate(mpv3,by=list(mpv3$KEY_PATIENT),FUN=function(x) max(x))
mpv00<-mpv0[mpv0$tag==0,]
t0<-aggregate(mpv00,by=list(mpv00$KEY_PATIENT),FUN=function(x) min(x))
ta<-rbind(t0,t1,t2,t3)

#reshape data frame
q1<-ta %>% select(3,4,7)
q1$tag2<-paste(q1$tag,"val",sep="_")
q2<-ta %>% select(3,5,7)
q2$tag2<-paste(q2$tag,"dttm",sep="_")
colnames(q2)[2]<-colnames(q1)[2]
#q1$val<-as.character(q1$val)
#q2$val<-as.character(q2$val)
qa<-rbind(q1,q2)
colnames(qa)[2]<-"val"
qt<-reshape2::dcast(qa,KEY_PATIENT~tag2,value.var="val")
qt$`0_dttm`<-as.POSIXct(qt$`0_dttm`, origin = "1970-01-01", tz = "UTC") 
qt$`1_dttm`<-as.POSIXct(qt$`1_dttm`, origin = "1970-01-01", tz = "UTC") 
qt$`2_dttm`<-as.POSIXct(qt$`2_dttm`, origin = "1970-01-01", tz = "UTC") 
qt$`3_dttm`<-as.POSIXct(qt$`3_dttm`, origin = "1970-01-01", tz = "UTC") 
colnames(qt)<-c("KEY_PATIENT","datebaseline","valuebaseline","date0h6h","value0h6h","date6h24h","value6h24h","date24h48h","value24h48h")
return(qt)
}
# R automatically use CEST instead of UTC, convert it manually
#######Processing##########

mpv<-read_excel("303_MPV.xlsx")
mpv$val<-as.numeric(mpv$val)
mpv<-timecount(mpv,6)
qt<-convertch(mpv,4,6)
fwrite(qt,"MPV.res")

thro<-read_excel("303_THRO.xlsx")
thro$val<-as.numeric(thro$val)
thro<-timecount(thro,6)
qt<-convertch(thro,4,6)
fwrite(qt,"THRO.res")

madp<-read_excel("303_MPADP.xlsx")
madp$LABVAL<-as.numeric(madp$LABVAL)
madp<-timecount(madp,9)
qt<-convertch(madp,5,9)
fwrite(qt,"MPADP.res")

mpaspi<-read_excel("303_MPASPI.xlsx")
mpaspi$LABVAL<-as.numeric(mpaspi$LABVAL)
mpaspi<-timecount(mpaspi,9)
qt<-convertch(mpaspi,5,9)
fwrite(qt,"MPASPI.res")

ipf<-read_excel("303_IPF.xlsx")
ipf$LABVAL<-as.numeric(ipf$LABVAL)
ipf<-timecount(ipf,9)
qt<-convertch(ipf,5,9)
fwrite(qt,"IPF.res")

ipfab<-read_excel("303_IPFAB.xlsx")
ipfab$LABVAL<-as.numeric(ipfab$LABVAL)
ipfab<-timecount(ipfab,9)
qt<-convertch(ipfab,5,9)
fwrite(qt,"IPFAB.res")
