setwd("/home/pang/Desktop/CAD_IS_PAD_20200121/data/METAL/final")
#load library
library(data.table)
library("qvalue")
library("AWFisher")
library(RACER)
# read GWAS stat 
# METAL meta
# fourier region
cadmeta<-fread("CAD_fourier_region")
ismeta<-fread("IS_fourier_region")
padmeta<-fread("PAD_fourier_region")
colnames(cadmeta)<-c("chr","BP","SNP","A1","A2","BETA","SE","P","direction","r0")
colnames(ismeta)<-c("chr","BP","SNP","A1","A2","BETA","SE","P","direction","r0")
colnames(padmeta)<-c("chr","BP","SNP","A1","A2","BETA","SE","P","direction","r0")

# plink meta
# 500 kb
cadmeta<-fread("CAD_GRCh37_500kb")
ismeta<-fread("IS_GRCh37_500kb")
padmeta<-fread("PAD_GRCh37_500kb")
colnames(cadmeta)<-c("chr","BP","SNP","A1","A2","N","P","P(R)","OR","OR(R)","Q","I","direction","r0")
colnames(ismeta)<-c("chr","BP","SNP","A1","A2","N","P","P(R)","OR","OR(R)","Q","I","direction","r0")
colnames(padmeta)<-c("chr","BP","SNP","A1","A2","N","P","P(R)","OR","OR(R)","Q","I","direction","r0")

cadmeta$MAF<-0.5
ismeta$MAF<-0.5
padmeta$MAF<-0.5
cadmeta$N<-643087
ismeta$N<-879922
padmeta$N<-688564
cadmeta$Z<-cadmeta$BETA/cadmeta$SE
ismeta$Z<-ismeta$BETA/ismeta$SE
padmeta$Z<-padmeta$BETA/padmeta$SE


#PRS
#cadmeta$N<-197583
#cadmeta$N<-184305 cardioonly
#ismeta$N<-455952
#padmeta$N<-243060
# uniform reference allele
# library(bigsnpr)
# not necessary using METAL result
# adjusted pvalue

CADqvalue<-qvalue(cadmeta$P,fdr.level=0.05,pfdr=FALSE)
cadmeta$qvalue<-CADqvalue$qvalues
ISqvalue<-qvalue(ismeta$P,fdr.level=0.05,pfdr=FALSE)
ismeta$qvalue<-ISqvalue$qvalues
PADqvalue<-qvalue(padmeta$P,fdr.level=0.05,pfdr=FALSE)
padmeta$qvalue<-PADqvalue$qvalues

cadmeta$BETA<-log(cadmeta$OR)
ismeta$BETA<-log(ismeta$OR)
padmeta$BETA<-log(padmeta$OR)

#manhattanplot
cadm<-cadmeta[cadmeta$qvalue<0.05,]
ism<-ismeta[ismeta$qvalue<0.05,]
padm<-padmeta[padmeta$qvalue<0.05,]
fwrite(cadm,"CAD_FDR_GRCh37_500kb")
fwrite(ism,"IS_FDR_GRCh37_500kb")
fwrite(padm,"PAD_FDR_GRCh37_500kb")


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
  res$cadp[i]<-min(abs(cadp$Z))
  res$cadb[i]<-cadp[cadp$Z==res$cadp[i],6]
  isp<-ismeta[ismeta$r0==region$r0[i]]
  res$isp[i]<-min(abs(isp$Z))
  res$isb[i]<-isp[isp$Z==res$isp[i],6]
  padp<-padmeta[padmeta$r0==region$r0[i],]
  res$padp[i]<-min(abs(padp$Z))
  res$padb[i]<-padp[padp$Z==res$padp[i],6]
  res$r0[i]<-region$r0[i]
  # Second Stage
  # Fisher statistics in an intervel
  #FS = -2*sum(log(pi))
  #res$fs[i]<--2*(log(res$cadp[i])+log(res$isp[i])+log(res$padp[i]))
  res$fs[i]<-1-pnorm(sum(res$cadp[i]*res$cadb[i],res$isp[i]*res$isb[i],res$padp[i]*res$padp[i])/sqrt(sum(res$cadb[i]^2,res$isb[i]^2,res$padb[i]^2)))
}
# Remove invalid rows
res <- res[!is.infinite(res$fs),]
res<-na.omit(res)
write.table(res,"METAL_fourier_region.res",row.names = FALSE,quote=FALSE)
write.table(res,"METAL_GRCh37_500kb.res",row.names = FALSE,quote=FALSE)

# AWFisher

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

#flip all risk alleles
t[t$beta<0,c("A1","A2")]<-t[t$beta<0,c("A2","A1")]
t$beta<-abs(t$beta)
test<-mutate(test,max=pmax(cadbeta,isbeta,padbeta))

## Function plot prev; per prev -wGRS
wanalysis<-function(dat,pheno){
  
  modComp  <- summary( glm(CC~wgrs,family="binomial",data=dat) )$coefficients
  #dat$grp  <- round(dat$NRA) #
  dat$grp<-cut(dat$wgrs,breaks = seq(min(dat$wgrs),max(dat$wgrs),by=0.05),include.lowest=TRUE)
  prev_p   <- aggregate(dat$CC,by=list(dat$grp),FUN=mean,na.rm=T)
  prev_n   <- aggregate(dat$CC,by=list(dat$grp),FUN=sum,na.rm=T)
  prev_N   <- aggregate(dat$CC,by=list(dat$grp),FUN=function(x) length(!is.na(x)))
  #allele_c <- aggregate(dat$NRA,by=list(dat$grp),FUN=mean)
  allele_c <- aggregate(dat$wgrs,by=list(dat$grp),FUN=mean)
  wgrs<-aggregate(dat$wgrs,by=list(dat$grp),FUN=mean)
  d        <- cbind.data.frame(N=prev_N$x,n=prev_n$x,p=prev_p$x,numAllele=allele_c$x,wgrs=wgrs$x)
  d$p_se   <- sqrt(d$p*(1-d$p)/d$N)
  d        <- d[which(d$N>=200),]
  d$px     <- d$p / d$numAllele
  d$px_se  <- d$p_se / d$numAllele
  
  #modelLogit  <- glm(cbind(n,N)~wgrs,data=d,family = binomial(link="logit"))
  #modelProbit <- glm(cbind(n,N)~wgrs,data=d,family = binomial(link="probit"))
  #modelLog    <- glm(cbind(n,N)~wgrs,data=d,family = binomial(link="log"))
  #modelLinear <- lm(p~wgrs,data=d,weights = N)
  
  modelLogit  <- glm(p~wgrs,data=d,family = binomial(link="logit"))
  modelProbit <- glm(p~wgrs,data=d,family = binomial(link="probit"))
  modelLog    <- glm(p~wgrs,data=d,family = binomial(link="log"))
  modelLinear <- lm(p~wgrs,data=d,weights = N)
  
  predLogit   <- predict(modelLogit,type = "response")
  predProbit  <- predict(modelProbit,type = "response")
  predLog     <- predict(modelLog,type = "response")
  predLin     <- predict(modelLinear,type="response")
  
  Rlogit      <- round(unlist(cor.test(d$p,predLogit)[c("estimate","conf.int")]),2)
  Rprobit     <- round(unlist(cor.test(d$p,predProbit)[c("estimate","conf.int")]),2)
  Rlog        <- round(unlist(cor.test(d$p,predLog)[c("estimate","conf.int")]),2)
  Rlin        <- round(unlist(cor.test(d$p,predLin)[c("estimate","conf.int")]),2)
  
  Cols <- c("seagreen","coral1","goldenrod")
  l    <- which(d[,"p_se"]>0)
  
  #  png(paste0("plots/",ds,".prevalence.fit.png"),width=1500,height=1250,res=250)
  #pdf("test.pdf")
  # plot(d[l,"numAllele"],d[l,"p"],pch=19,axes=FALSE,ylim=c(0,max(d$p+d$p_se)),main=ds,xlab="Mean number of risk alleles",ylab="Prevalence");axis(1);axis(2)
  plot(d[l,"wgrs"],d[l,"p"],pch=19,axes=FALSE,ylim=c(0,max(d$p+d$p_se)),xlim=c(min(d$wgrs-sd(d$wgrs)),max(d$wgrs+sd(d$wgrs))),main=pheno,xlab="Mean number of weighted GRS",ylab="Prevalence");axis(1);axis(2)
  for(i in 1:nrow(d)){
    if(d[i,"p_se"]>0){
      arrows(d[i,"wgrs"],d[i,"p"]-d[i,"p_se"],d[i,"wgrs"],d[i,"p"]+d[i,"p_se"],col="grey",angle=90,len=0.05)
      arrows(d[i,"wgrs"],d[i,"p"]+d[i,"p_se"],d[i,"wgrs"],d[i,"p"]-d[i,"p_se"],col="grey",angle=90,len=0.05)
    }
  }
  #  matlines(d$numAllele,cbind(predLogit,predProbit,predLog),col=Cols,lwd=2)
  #  abline(a=modelLinear$coefficients[1],b=modelLinear$coefficients[2],col=4,lty=3)
  matlines(d$wgrs,cbind(predLogit,predLog),col=Cols,lwd=2)
  abline(a=modelLog$coefficients[1],b=modelLog$coefficients[2],col=4,lty=3)
  quant<-quantile(dat$wgrs,probs=c(0.1,0.9))
  polygon(x=c(quant[1],quant[1],quant[2],quant[2]),y=c(0,max(d$p),max(d$p),0),col=rgb(0.2,0.2,0.2,0.2),border = NA)
  legend("topleft",
         legend=c(paste0("Logit: R=",Rlogit[1]," (95%CI: [",Rlogit[2],"-",Rlogit[3],"])"),
                  #                paste0("Probit: R=",Rprobit[1]," (95%CI: [",Rprobit[2],"-",Rprobit[3],"])"),
                  paste0("Log: R=",Rlog[1]," (95%CI: [",Rlog[2],"-",Rlog[3],"])")),
         #                paste0("Linear: R=",Rlin[1]," (95%CI: [",Rlin[2],"-",Rlin[3],"])")),
         #       col=c(Cols,"blue"),cex=0.9,box.lty=0,lty=c(1,1,1,3),lwd=2,bg="transparent")
         col=Cols,cex=0.9,box.lty=0,lwd=2,bg="transparent")
  #dev.off()
}



t<-fread("~/Desktop/CAD_IS_PAD_20200121/data/Final/PRS/share_flip/UKBB_PAD.all.score")
t2<-t%>% select("FID","0.00050005")
t<-RCV[match(t2$FID,RCV$schunkert),]
t2$PAD<-t$PAD
t2<-t2[is.na(t2$PAD)==FALSE,]
colnames(t2)<-c("sample","wgrs","CC")
dat<-t2

cairo_pdf("~/Desktop/CAD_IS_PAD_20200121/data/Final/figures/UKBB_PAD.pdf",width=4.63,height=5.64,pointsize=12,family="Univers")
wanalysis(dat,"UKBB PAD")
dev.off()
