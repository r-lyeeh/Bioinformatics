# library path correction
myPaths<-.libPaths()
newpath<-"/raid3/software/centos-7/R/R-4.0.1/library"
myPaths<-c(myPaths[1],myPaths[2],newpath)
.libPaths(myPaths)
            
setwd("/home/pang/Desktop/CAD_IS_PAD_20200121/data")
# read GWAS stat 
cadmeta<-read.table("CAD/CAD_meta.txt",header=TRUE,sep=" ",as.is=TRUE)
ismeta<-read.table("IS/IS_meta.txt",header=TRUE,sep=" ",as.is=TRUE)
padmeta<-read.table("PAD/PAD_meta.txt",header=TRUE,sep=" ",as.is=TRUE)
# Get beta
cadmeta$beta<-log(cadmeta$OR)
ismeta$beta<-log(ismeta$OR)
padmeta$beta<-log(padmeta$OR)
# adjusted p-values
library("qvalue")
CAD_qvalue<-qvalue(cadmeta$P,fdr.level = 0.05,pfdr = FALSE)
cadmeta$p_adjust<-CAD_qvalue$qvalues
IS_qvalue<-qvalue(ismeta$P,fdr.level = 0.05,pfdr = FALSE)
ismeta$p_adjust<-IS_qvalue$qvalues
PAD_qvalue<-qvalue(padmeta$P,fdr.level = 0.05,pfdr = FALSE)
padmeta$p_adjust<-PAD_qvalue$qvalues
# Set direction
cadmeta[cadmeta$beta<0,15]<-"-"
ismeta[ismeta$beta<0,15]<-"-"
padmeta[padmeta$beta<0,15]<-"-"
cadmeta[cadmeta$beta>0,15]<-"+"
ismeta[ismeta$beta>0,15]<-"+"
padmeta[padmeta$beta>0,15]<-"+"
colnames(cadmeta)[15]<-"Direction"
colnames(ismeta)[15]<-"Direction"
colnames(padmeta)[15]<-"Direction"
cadmeta$MAF<-0.5
ismeta$MAF<-0.5
padmeta$MAF<-0.5
cadmeta$N<-643087
ismeta$N<-879922
padmeta$N<-688564
colnames(cadmeta)[1:3]<-c("Chr","BP","Marker")
colnames(ismeta)[1:3]<-c("Chr","BP","Marker")
colnames(padmeta)[1:3]<-c("Chr","BP","Marker")
colnames(cadmeta)[7]<-c("pValue")
colnames(ismeta)[7]<-c("pValue")
colnames(padmeta)[7]<-c("pValue")
# Convert SE
cadmeta$SE<-abs(cadmeta$beta/qnorm(cadmeta$pValue/2))
ismeta$SE<-abs(ismeta$beta/qnorm(ismeta$pValue/2))
padmeta$SE<-abs(padmeta$beta/qnorm(padmeta$pValue/2))
# Using bmass
library("bmass")
# Load GWAS statistics
library("tidyverse")
# Format: Chr BP Marker MAF A1 A2 Direction pValue N
# Direction: all +
# N: Simulate - all
# MAF: Simulate - 0.5
bmass_cad<-cadmeta %>% select("Chr","BP","Marker","MAF","A1","A2","Direction","pValue","N")
# G5:3929 G6:2825 G7:6524 CARDIOGRAMplusC4D:60801+123504=184305 UKBB:445504
#bmass_is<-ismeta %>% select("Chr","BP","Marker","MAF","A1","A2","Direction","pValue","N")
# megastroke:434418 UKBB:445504
#bmass_pad<-padmeta %>% select("Chr","BP","Marker","MAF","A1","A2","Direction","pValue","N")
# MVP:243060 UKBB:445504
#write.table(bmass_cad,"bmass_cad",sep="\t",row.names=FALSE,quote=FALSE)
#write.table(bmass_is,"bmass_is",sep="\t",row.names=FALSE,quote=FALSE)
#write.table(bmass_pad,"bmass_pad",sep="\t",row.names=FALSE,quote=FALSE)
bmass_is<-read.table("bmass_is_match_cad",sep="\t",header=TRUE,as.is=TRUE)
bmass_pad<-read.table("bmass_pad_match_cad",sep="\t",header=TRUE,as.is=TRUE)
# Be careful of colname lines
bmass_cad<-bmass_cad %>% drop_na()
bmass_is<-bmass_is %>% drop_na()
bmass_pad<-bmass_pad %>% drop_na()
# Select SNPs
Phenotypes <- c("bmass_cad", "bmass_is","bmass_pad")
# Load sig SNPs Chr BP
bmass_sigSNPs<-read.table("SNP/cadispad_fdr_raw",header=TRUE,as.is=TRUE,sep=" ")
# Excute bmass
# too many SNPs, select 2000000 SNPs
# as.numeric
bmass_cad[, c(1,2,4,8,9)] <- sapply(bmass_cad[, c(1,2,4,8,9)], as.numeric)
bmass_is[, c(1,2,4,8,9)] <- sapply(bmass_is[, c(1,2,4,8,9)], as.numeric)
bmass_pad[, c(1,2,4,8,9)] <- sapply(bmass_pad[, c(1,2,4,8,9)], as.numeric)
library('dplyr')
bmass_cad2<- arrange(bmass_cad, pValue)
bmass_cad2<-bmass_cad2[1:1000000,]
# remove duplicate
# error? arrange is not correct
bmass_is2<-bmass_is
bmass_is2<-bmass_is2[!duplicated(bmass_is2$Marker), ]
bmass_is2<- arrange(bmass_is2, pValue)
bmass_is2<-bmass_is2[1:1000000,]
bmass_pad2<- arrange(bmass_pad, pValue)
bmass_pad2<-bmass_pad2[1:1000000,]
#bmass_cad2<-bmass_cad[bmass_cad$pValue<0.05,]
#bmass_is2<-bmass_is[bmass_is$pValue<0.05,]
#bmass_pad2<-bmass_pad[bmass_pad$pValue<0.05,]
Phenotypes <- c("bmass_cad2", "bmass_is2","bmass_pad2")
bmassResults <- bmass(Phenotypes, bmass_sigSNPs)
newsnp<-bmassResults$NewSNPs$SNPs
#write.table(newsnp,"multitrait_newSNP",row.names=FALSE,quote=FALSE,sep="\t")

# Add SE
cad2<-merge(bmass_cad2,cadmeta,by="Marker")
is2<-merge(bmass_is2,ismeta,by="Marker")
pad2<-merge(bmass_pad2,padmeta,by="Marker")
cad1<-cad2 %>% select("Chr.x","BP.x","Marker","A1.x","A2.x","MAF.x","Direction.x","pValue.x","N.x","p_adjust","SE","beta")
colnames(cad1)<-c("Chr","BP","Marker","A1","A2","MAF","Direction","pValue","N","p_adjust","SE","BETA")
cad1$betac<-cad1$BETA
is1<-is2 %>% select("Chr.x","BP.x","Marker","A1.x","A2.x","MAF.x","Direction.x","pValue.x","N.x","p_adjust","SE","beta")
colnames(is1)<-c("Chr","BP","Marker","A1","A2","MAF","Direction","pValue","N","p_adjust","SE","BETA")
is1.1<-is1[is1$Direction=="+",]
is1.1$betac<-abs(is1.1$BETA)
is1.2<-is1[is1$Direction=="-",]
is1.2$betac<--abs(is1.2$BETA)
is1<-rbind(is1.1,is1.2)
pad1<-pad2 %>% select("Chr.x","BP.x","Marker","A1.x","A2.x","MAF.x","Direction.x","pValue.x","N.x","p_adjust","SE","beta")
colnames(pad1)<-c("Chr","BP","Marker","A1","A2","MAF","Direction","pValue","N","p_adjust","SE","BETA")
pad1.1<-pad1[pad1$Direction=="+",]
pad1.1$betac<-abs(pad1.1$BETA)
pad1.2<-pad1[pad1$Direction=="-",]
pad1.2$betac<--abs(pad1.2$BETA)
pad1<-rbind(pad1.1,pad1.2)
write.table(cad1,"cadmeta",row.names=FALSE,quote=FALSE,sep="\t")
write.table(is1,"ismeta",row.names=FALSE,quote=FALSE,sep="\t")
write.table(pad1,"padmeta",row.names=FALSE,quote=FALSE,sep="\t")

library(CPBayes)
# Using UK biobank
traitNames<-c("CAD","IS","PAD") # trait names
BetaHat<-c(0.0577970987262168,0.0602809308849299,0.0593061058974075) # vector of betas for each phenotype
SE<-c(0.0100516579209148,0.009842704,0.008809187) # vector of standard errors
SNP1<-"rs11066301" # lead SNP name
result <- cpbayes_uncor(BetaHat, SE, Phenotypes = traitNames, Variant = SNP1)
PleioSumm <- post_summaries(result, level = 0.05)
x<-forest_cpbayes(result, level = 0.05)
for (i in 1:nrow(leadsnp)){
  traitNames<-c("CAD","IS","PAD") # trait names
  BetaHat<-betares[i,] # vector of betas for each phenotype
  SE<-seres[i,] # vector of standard errors
  SNP1<-snpres[i,] # lead SNP name
  result <- cpbayes_uncor(BetaHat, SE, Phenotypes = traitNames, Variant = SNP1)
  PleioSumm <- post_summaries(result, level = 0.05)
  x<-forest_cpbayes(result, level = 0.05)
}

#Moloc for all
#SNP  BETA  SE
library(moloc)
cadx1<-cad1 %>% select("Marker","betac","SE","N","MAF")
colnames(cadx1)<-c("SNP","BETA","SE","N","MAF")
isx1<-is1 %>% select("Marker","betac","SE","N","MAF")
colnames(isx1)<-c("SNP","BETA","SE","N","MAF")
padx1<-pad1 %>% select("Marker","betac","SE","N","MAF")
colnames(padx1)<-c("SNP","BETA","SE","N","MAF")
# isx1 problem
#isx1<-isx1[!duplicated(isx1$SNP), ]
Phenotypes<-list(padx1,isx1,padx1)
moloc <- moloc_test(Phenotypes)
y<-moloc[1]$priors_lkl_ppa
write.table(y,"moloc.ppa",quote=FALSE,sep="\t")

#metabf
library("metabf")
betas<-cbind(cadx1$BETA,isx1$BETA,padx1$BETA)
betas <- cotsapas[,grep("\\.beta", names(cotsapas))]
#Collect the summary data from the cotsapas data frame.
ses <- cotsapas[,grep("\\.SE", names(cotsapas))]
meta.abf(betas, ses, prior.sigma = 0.2, prior.cor = "indep", log10 = TRUE)


library("tidyverse")
setwd("/home/pang/Desktop/CAD_IS_PAD_20200121/data/SNP")
cad<-read.table("cadmeta",as.is=TRUE,header=TRUE)
is<-read.table("ismeta",as.is=TRUE,header=TRUE)
pad<-read.table("padmeta",as.is=TRUE,header=TRUE)
cads<-arrange(cad,pValue)
# select top 0.05% SNPs
cadx<-cads[1:ceiling(nrow(cads)/2000),]
iss<-arrange(is,pValue)
isx<-iss[1:ceiling(nrow(iss)/2000),]
pads<-arrange(pad,pValue)
padx<-pads[1:ceiling(nrow(pads)/2000),]

write.table(cadx,"5000/cadx",quote=FALSE,row.names = FALSE,sep="\t")
write.table(isx,"5000/isx",quote=FALSE,row.names = FALSE,sep="\t")
write.table(padx,"5000/padx",quote=FALSE,row.names = FALSE,sep="\t")

library(lattice)
library(qqman)
manhattan(cad, chr="Chr", bp="BP", snp="Marker", p="pValue" )
