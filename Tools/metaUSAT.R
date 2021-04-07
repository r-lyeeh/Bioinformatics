# library path correction
myPaths<-.libPaths()
newpath<-"/raid3/software/centos-7/R/R-4.0.1/library"
myPaths<-c(myPaths[1],myPaths[3])
.libPaths(myPaths)

setwd("/home/pang/Desktop/CAD_IS_PAD_20200121/data")
# read GWAS stat 
cadmeta<-read.table("CAD/CAD_meta.txt",header=TRUE,sep=" ",as.is=TRUE)
ismeta<-read.table("IS/IS_meta.txt",header=TRUE,sep=" ",as.is=TRUE)
padmeta<-read.table("PAD/PAD_meta.txt",header=TRUE,sep=" ",as.is=TRUE)
# Get beta
# library("QCGWAS")
cadmeta$beta<-log(cadmeta$OR)
ismeta$beta<-log(ismeta$OR)
padmeta$beta<-log(padmeta$OR)
# uniform reference allele
library(bigsnpr)
# using cad as template; "chr", "pos", "a0" and "a1".
info_snp<-cadmeta[,1:5]
colnames(info_snp)<-c("chr","pos","id","a0","a1")
# convert sumstats "chr", "pos", "a0", "a1" and "beta"
ist<-ismeta[,1:5]
ist<-cbind(ist,ismeta$beta)
colnames(ist)<-c("chr","pos","id","a0","a1","beta")
ist_flip<-snp_match(ist,info_snp,strand_flip = FALSE,join_by_pos = TRUE,match.min.prop = 0.5)
# strand flip=FALSE, a0-a1 different among different trait
istmeta_flip<-ismeta[match(ist_flip$id.ss,ismeta$SNP),]
istmeta_flip$beta<-ist_flip$beta
padt<-padmeta[,1:5]
padt<-cbind(padt,padmeta$beta)
colnames(padt)<-c("chr","pos","id","a0","a1","beta")
padt_flip<-snp_match(padt,info_snp,strand_flip = FALSE,join_by_pos = TRUE,match.min.prop = 0.5)
# strand flip=FALSE, a0-a1 different among different trait
padtmeta_flip<-padmeta[match(padt_flip$id.ss,padmeta$SNP),]
padtmeta_flip$beta<-padt_flip$beta
# multivariate meta-analysis
require(devtools)
source_url("https://github.com/RayDebashree/metaUSAT/blob/master/metaUSAT_v1.17.R?raw=TRUE")
# Z.matrix K x p K-trait;row-snp
# p.threshold default 10-5
library(tidyverse)
cadp<-cadmeta %>% select(SNP, P,beta)
isp<-istmeta_flip %>% select(SNP,P,beta)
padp<-padtmeta_flip %>% select(SNP,P,beta)
cadis<-merge(cadp,isp,by.x = "SNP",by.y = "SNP")
colnames(cadis)<-c("SNP","CADp","CADb","ISp","ISb")
cadispad<-merge(cadis,padp,by.x="SNP",by.y="SNP")
colnames(cadispad)[6:7]<-c("PADp","PADb")
cadispad2<-cadispad
cadispad<-cadispad2 %>% distinct()
# non-unique values when setting 'row.names': ‘rs12787329’, ‘rs4130054’, ‘rs4301800’, ‘rs6684987’
# manually check
# 1573342 rs12787329 0.8837 -0.001300846 0.9503 -0.0006001801x 0.1883 -0.01222497
# 3779425 rs4130054 0.9822 -0.000100005x 0.4542 -0.007326775 0.5851 -0.004589452
# 3825018 rs4301800 0.9782 0.000299955 0.05124  0.01882175 0.3002  0.008662373
# 4911613 rs6684987 0.008983 -0.01816397 0.8625 -0.002002003 0.9593 -0.000500125x
rownames(cadispad)<-cadispad[,1]
cadispad<-cadispad[,-1]
# P.matrix K x p K-trait;row-snp;cell-p-value
# adjusted p-values
library("qvalue")
cadispadq<-cadispad
CAD_qvalue<-qvalue(cadispad$CADp,fdr.level = 0.05,pfdr = FALSE)
cadispadq$CADp<-CAD_qvalue$qvalues
IS_qvalue<-qvalue(cadispad$ISp,fdr.level = 0.05,pfdr = FALSE)
cadispadq$ISp<-IS_qvalue$qvalues
PAD_qvalue<-qvalue(cadispad$PADp,fdr.level = 0.05,pfdr = FALSE)
cadispadq$PADp<-PAD_qvalue$qvalues
rownames(cadispadq)<-cadispadq[,1]
cadispadq<-cadispadq[,-1]
# construct Z matrix ~ beta effect needs flip 
Zmatrix<-cadispadq %>% select(CADb,ISb,PADb,SNP)
colnames(Zmatrix)<-c("CAD","IS","PAD","SNP")
# prepare R matrix; 
Qmatrix<-cadispadq %>% select(CADp,ISp,PADp,SNP)
colnames(Qmatrix)<-c("CAD","IS","PAD","SNP")
# metaUST method seems not working well
# Too many SNPs, pick out p<0.1
cad2<-Qmatrix[Qmatrix$CAD<0.1,]
is2<-Qmatrix[Qmatrix$IS<0.1,]
pad2<-Qmatrix[Qmatrix$PAD<0.1,]
Qmatrix2<-rbind(cad2,is2,pad2)
Zmatrix2<-Zmatrix[match(Qmatrix2$SNP,Zmatrix$SNP),]
Qmatrix2<-Qmatrix2[,-4]
Zmatrix2<-Zmatrix2[,-4]
# R<-cor.pearson(cadispadq,cadispadq,p.threshold=0.05)
R2<-cor.pearson(Zmatrix2,Qmatrix2,p.threshold=0.05)
res<-Zmatrix2
colnames(res)<-c("T","P","SNP")
for (i in 1:100){
  zi<-as.matrix(Zmatrix2[i,])
  out<-metausat(Z=zi,R=R2,weights=1)
  res[i,1]<-out$T.metausat
  res[i,2]<-out$p.metausat
  res[i,3]<-rownames(Zmatrix2)[i]
}



