#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=TRUE)
library(ggplot2)
library(data.table)
setwd("/home/pang/Desktop/CAD_IS_PAD_20200121/data/METAL/final")
# trace("ldRACER",edit=TRUE)
# change to ldlink.nci.nih.gov
region<-args[1]
specific<-args[2]

moloc<-fread(paste("moloc.ppa_snp",region,sep="_"))
p<-moloc[moloc$region==specific,]
p$group<-c("CAD","IS","CADIS","PAD","CADPAD","ISPAD","CADISPAD")
leadsnp0<-p[7,4][[1]]
moloc2<-fread(paste("moloc.ppa",region,sep="_"))
p2<-moloc2[moloc2$region==specific,]
p2$Group<-c("CAD","CAD,IS","CAD,PAD","CAD,ISPAD","CAD,IS,PAD","IS","IS,PAD","CADPAD,IS","PAD","CADIS,PAD","CADIS","CADPAD","ISPAD","CADISPAD","zero")

#p3<-p2 %>% select(group,PPA)
#row.names(p3)<-p3$group
#p3<-p3[,-1]
#p3<-as.matrix(p3)
#barplot(t(p3),col = "red",horiz = TRUE,las=1)

library(CPBayes)
library(tidyverse)
# Using UK biobank
meta<-fread("CADISPAD_fourier_region")
traitNames<-c("CAD","IS","PAD") # trait names
t1<-meta[meta$SNP==leadsnp0,]
BetaHat<-t1 %>% select("CAD_BETA","IS_BETA","PAD_BETA") # vector of betas for each phenotype
SE<-t1 %>% select("CAD_SE","IS_SE","PAD_SE") # vector of standard errors
SNP1<-leadsnp0 # lead SNP name
BetaHat<-as.numeric(abs(BetaHat))
SE<-as.numeric(SE)
result <- cpbayes_uncor(BetaHat, SE, Phenotypes = traitNames, Variant = SNP1)
PleioSumm <- post_summaries(result, level = 0.05)
x<-forest_cpbayes(result, level = 0.05)
pdf(paste(args[1],args[2],leadsnp0,".pdf",sep=""),onefile = F)
ggplot(data=p2, aes(x=Group, PPA)) + geom_bar(fill="red",stat="identity")+ coord_flip()+theme_classic()
dev.off()
