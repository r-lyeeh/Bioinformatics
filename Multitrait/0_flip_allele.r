setwd("/home/pang/Desktop/CAD_IS_PAD_20200121/data/METAL/final")
library(data.table)
# read GWAS stat 
cadmeta<-fread("CAD_fourier_region")
ismeta<-fread("IS_fourier_region")
padmeta<-fread("PAD_fourier_region")
cadmeta<-fread("CAD_GRCh37_500kb")
ismeta<-fread("IS_GRCh37_500kb")
padmeta<-fread("PAD_GRCh37_500kb")
colnames(cadmeta)<-c("chr","BP","SNP","A1","A2","BETA","SE","P","direction","r0")
colnames(ismeta)<-c("chr","BP","SNP","A1","A2","BETA","SE","P","direction","r0")
colnames(padmeta)<-c("chr","BP","SNP","A1","A2","BETA","SE","P","direction","r0")

colnames(cadmeta)<-c("chr","BP","SNP","A1","A2","N","P","P(R)","OR","OR(R)","Q","I","direction","r0")
colnames(ismeta)<-c("chr","BP","SNP","A1","A2","N","P","P(R)","OR","OR(R)","Q","I","direction","r0")
colnames(padmeta)<-c("chr","BP","SNP","A1","A2","N","P","P(R)","OR","OR(R)","Q","I","direction","r0")

cadmeta$MAF<-0.5
ismeta$MAF<-0.5
padmeta$MAF<-0.5
cadmeta$N<-643087
ismeta$N<-879922
padmeta$N<-688564
#PRS
#cadmeta$N<-197583
#cadmeta$N<-184305 cardioonly
#ismeta$N<-455952
#padmeta$N<-243060
# uniform reference allele
# library(bigsnpr)
# not necessary using METAL result
# adjusted pvalue
library("qvalue")
CADqvalue<-qvalue(cadmeta$P,fdr.level=0.05,pfdr=FALSE)
cadmeta$qvalue<-CADqvalue$qvalues
ISqvalue<-qvalue(ismeta$P,fdr.level=0.05,pfdr=FALSE)
ismeta$qvalue<-ISqvalue$qvalues
PADqvalue<-qvalue(padmeta$P,fdr.level=0.05,pfdr=FALSE)
padmeta$qvalue<-PADqvalue$qvalues

cadmeta$BETA<-log(cadmeta$OR)
ismeta$BETA<-log(ismeta$OR)
padmeta$BETA<-log(padmeta$OR)

# Generate RMS method
region<-read.table("GRCh37_500kb",header=TRUE,as.is=TRUE)
#region<-read.table("fourier_region",header=TRUE,as.is=TRUE)
res<-data.frame(matrix(nrow=nrow(region), ncol=5))
colnames(res)<-c("cadp","isp","padp","r0","fs")
for (i in 1:nrow(region)){
  # First stage
  # MP Tippet's minimum p-value
  # MP = min(1<k<K) (pk)
  # using the smallest p-value among three traits
  cadp<-cadmeta[cadmeta$r0==region$r0[i],]
  res$cadp[i]<-min(cadp$qvalue)
  isp<-ismeta[ismeta$r0==region$r0[i],]
  res$isp[i]<-min(isp$qvalue)
  padp<-padmeta[padmeta$r0==region$r0[i],]
  res$padp[i]<-min(padp$qvalue)
  res$r0[i]<-region$r0[i]
  # Second Stage
  # Fisher statistics in an intervel
  #FS = -2*sum(log(pi)
  res$fs[i]<--2*(log(res$cadp[i])+log(res$isp[i])+log(res$padp[i]))
}
# Remove invalid rows
res <- res[!is.infinite(res$fs),]
res<-na.omit(res)
write.table(res,"METAL_fourier_region.res",row.names = FALSE,quote=FALSE)
write.table(res,"METAL_GRCh37_500kb.res",row.names = FALSE,quote=FALSE)

# AWFisher
library("AWFisher")
res<-read.table("/home/pang/Desktop/CAD_IS_PAD_20200121/data/RSM/MP_fdr2.res",header=TRUE,sep=" ",as.is=TRUE)
K <- 3 ## combining K studies
G <- nrow(res)  ## simulate G genes
p.values = matrix(runif(K*G), ncol=K)
cadq<-res$cadp
isq<-res$isp
padq<-res$padp
p.values[,1]<-cadq
p.values[,2]<-isq
p.values[,3]<-padq
resawf<- AWFisher_pvalue(p.values)
weight<-resawf$weights
weight<-as.data.frame(weight)
weight$sum<-rowSums( weight[,1:3] )
r0<-res$r0
fs<-res$fs
weight$r0<-r0
weight$fs<-fs
qvalue <- p.adjust(resawf$pvalue, "BH")
weight$p.adjust<-qvalue
#write.table(weight,"AWF_fdr2",row.names = FALSE,quote=FALSE)
weight$p<-resawf$pvalues
colnames(weight)[1:3]<-c("CAD","IS","PAD")
write.table(weight,"METAL_fourier_region.awf",row.names = FALSE,quote=FALSE)
write.table(weight,"METAL_GRCh37_500kb.awf",row.names = FALSE,quote=FALSE)
library(tidyverse)
library(moloc)
y0<-data.frame("prior"=0,"sumbf"=0,"logBF_locus"=0,"PPA"=0,"region"=0)
z0<-data.frame("region"=0,"coloc_ppas"=0,"best.snp.coloc"=0)
for (i in 1:nrow(region)){
  tryCatch({
    q<-region$r0[i]
    cadt<-cadmeta[cadmeta$r0==q,]
    cadx1<-cadt %>% select(SNP,BETA,SE,N,MAF)
    ist<-ismeta[ismeta$r0==q,]
    isx1<-ist %>% select(SNP,BETA,SE,N,MAF)
    padt<-padmeta[padmeta$r0==q,]
    padx1<-padt %>% select(SNP,BETA,SE,N,MAF)
    Phenotypes<-list(cadx1,isx1,padx1)
    moloc <- moloc_test(Phenotypes)
    y<-moloc$priors_lkl_ppa
    y$region<-q
    z<- moloc$best_snp
    z$region<-q
    y0<-rbind(y0,y)
    z0<-rbind(z0,z)
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
# error for no common SNPs

y0<-y0[-1,]
z0<-z0[-1,]
write.table(y0,"moloc.ppa_fourier_region",quote=FALSE,sep="\t")
write.table(z0,"moloc.ppa_snp_fourier_region",quote=FALSE,sep="\t")

write.table(y0,"moloc.ppa_GRCh37_500kb",quote=FALSE,sep="\t")
write.table(z0,"moloc.ppa_snp_GRCh37_500kb",quote=FALSE,sep="\t")

#stat
stat<-data.frame("CAD","IS","PAD","CADIS","CADPAD","ISPAD","CADISPAD")
#column: all_fourier;all_500kb;sig_fourier;sig_500kb;all_ISsnp_fourier;all_ISsnp_500kb;sig_ISsnp_fourier;sig_ISsnp_500kb
nrow(weight[weight$sum==3,])
nrow(weight[weight$CAD==1 & weight$IS==0 & weight$PAD==0,])
nrow(weight[weight$CAD==0 & weight$IS==0 & weight$PAD==1,])
nrow(weight[weight$CAD==0 & weight$IS==1 & weight$PAD==0,])
nrow(weight[weight$CAD==1 & weight$IS==1 & weight$PAD==0,])
nrow(weight[weight$CAD==0 & weight$IS==1 & weight$PAD==1,])
nrow(weight[weight$CAD==1 & weight$IS==0 & weight$PAD==1,])
sig<-weight[weight$p.adjust<0.05,]
nrow(sig[sig$sum==3,])
nrow(sig[sig$CAD==1 & sig$IS==0 & sig$PAD==0,])
nrow(sig[sig$CAD==0 & sig$IS==0 & sig$PAD==1,])
nrow(sig[sig$CAD==0 & sig$IS==1 & sig$PAD==0,])
nrow(sig[sig$CAD==1 & sig$IS==1 & sig$PAD==0,])
nrow(sig[sig$CAD==1 & sig$IS==0 & sig$PAD==1,])
nrow(sig[sig$CAD==0 & sig$IS==1 & sig$PAD==1,])

# Regional plot
library(RACER)
# Input data 
#gwas<-read.table("meta_flip_bp.all",header=TRUE,as.is=TRUE)
# Region chr8_19500000_20000000
cad12<-read.table("chr8_19500000_20000000_CAD",header=FALSE,as.is=TRUE)
is12<-read.table("chr8_19500000_20000000_IS",header=FALSE,as.is=TRUE)
pad12<-read.table("chr8_19500000_20000000_PAD",header=FALSE,as.is=TRUE)
cad_gwas_f = RACER::formatRACER(assoc_data = cad12, chr_col = 1, pos_col = 2, p_col = 7)
is_gwas_f = RACER::formatRACER(assoc_data = is12, chr_col = 1, pos_col = 2, p_col = 7)
pad_gwas_f = RACER::formatRACER(assoc_data = pad12, chr_col = 1, pos_col = 2, p_col = 7)
cad12[which.min(cad12$V7),]$V3
cad_gwas_f_ld = RACER::ldRACER(assoc_data = cad_gwas_f, rs_col = 3, pops = "EUR", lead_snp = "rs2197089")
is12[which.min(is12$V7),]$V3
is_gwas_f_ld = RACER::ldRACER(assoc_data = is_gwas_f, rs_col = 3, pops = "EUR", lead_snp = "rs1994786")
pad12[which.min(pad12$V7),]$V3
pad_gwas_f_ld = RACER::ldRACER(assoc_data = pad_gwas_f, rs_col = 3, pops = "EUR", lead_snp = "rs322")
pdf("chr8_19500000_20000000")
RACER::singlePlotRACER(assoc_data = cad_gwas_f_ld, chr =8, build = "hg19", plotby = "coord", start_plot = 19500000, end_plot = 20000000)
RACER::singlePlotRACER(assoc_data = is_gwas_f_ld, chr = 8, build = "hg19", plotby = "coord", start_plot = 19500000, end_plot = 20000000)
RACER::singlePlotRACER(assoc_data = pad_gwas_f_ld, chr = 8, build = "hg19", plotby = "coord", start_plot = 19500000, end_plot = 20000000)
mirrorPlotRACER(assoc_data1 = cad_gwas_f_ld, assoc_data2 = pad_gwas_f_ld, chr = 8,  build = "hg19", plotby = "coord", start_plot = 19500000, end_plot = 20000000)
scatterPlotRACER(assoc_data1 = cad_gwas_f_ld, assoc_data2 = pad_gwas_f_ld, chr = 8, name1 = "GWAS p_values", name2 = "19500000-113000000", region_start = 19500000, region_end = 20000000, ld_df = 1)
dev.off()

