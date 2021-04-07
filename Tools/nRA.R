setwd("/home/pang/data/SCP_involve/Liability_zeng_HS_FB_20190923/G7/JE_snp/PRSice")
library(data.table)
G5<-read.table("G5.JE.vcf",header = T,as.is=T)
G6<-read.table("G6.JE.vcf",header=T,as.is = T)
G7<-read.table("G7.JE.vcf",header = T,as.is = T)
NumberofRA<-function(ss){
  new<-ss[,-(1:9)]
  #setDT(G5_new_t)[,lapply(.SD, function(i) stringr::str_count(i,'1'))][]
  res<-data.frame(nchar(gsub("[^1]","",as.matrix(new))))
  sum<-colSums(res)
  res<-rbind(res,sum)
  info<-ss[,1:9]
  info[nrow(ss)+1,]<-rep(0,9)
  res<-cbind(info,res)
  return(res)
}
G5_res<-NumberofRA(G5)
G6_res<-NumberofRA(G6)
G7_res<-NumberofRA(G7)

                     


