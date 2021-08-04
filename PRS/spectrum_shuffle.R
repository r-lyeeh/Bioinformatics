t2<-bcdat2[bcdat2$uorder!=1&bcdat2$uorder!=10,]
t<-t2


prevs<-function(dat2l,adres){
  adres<-data.frame("bmihp"=0,"bmilp"=0,"bmimp"=0,"diff_bmi"=0)
  for (i in 1:10){
    t<-dat2l[[i]]
    tn<-t[is.na(t$alchol_freq_class)==FALSE,]
    tn2<-tn[is.na(tn$age)==FALSE,]
    tn<-tn2[tn2$age>=45,]
    bmih<-tn[tn$alchol_freq_class==2,]#daily
    bmil<-tn[tn$alchol_freq_class==1,]#once
    bmim<-tn[tn$alchol_freq_class==0,]#rare
    bmihp<-nrow(bmih[bmih$CC==1,])/nrow(bmih)
    bmilp<-nrow(bmil[bmil$CC==1,])/nrow(bmil)
    bmimp<-nrow(bmim[bmim$CC==1,])/nrow(bmim)
    diff_bmi<-bmihp-bmimp
    adres[i,]<-c(bmihp,bmilp,bmimp,diff_bmi)
  }
  return(adres)
}

dat2l<-split(bcdat2,bcdat2$uorder)
prevs(dat2l,adres)

prevs<-function(dat2l,adres){
  adres<-data.frame("bmihp"=0,"bmilp"=0,"diff_bmi"=0)
  for (i in 1:10){
    t<-dat2l[[i]]
    tn<-t[is.na(t$bmi)==FALSE,]
    tn2<-tn[is.na(tn$age)==FALSE,]
    tn<-tn2[tn2$age>=45,]
    bmih<-tn[tn$bmi>=30,]
    bmil<-tn[tn$bmi<30,]
    bmihp<-nrow(bmih[bmih$CC==1,])/nrow(bmih)
    bmilp<-nrow(bmil[bmil$CC==1,])/nrow(bmil)
    diff_bmi<-bmihp-bmilp
    adres[i,]<-c(bmihp,bmilp,diff_bmi)
  }
  return(adres)
}
dat2l<-split(bcdat3,bcdat3$uorder)
prevs(dat2l,adres)

bmihp
bmilp
bmimp
bmihp-bmilp

a<-seq(1,10,1)
b<-c(0.000721128,0.004160812,0.006870122,0.003904344,0.003650291,0.001426132,0.006181104,0.006316737,0.007006707,0.008442835)
res<-aov(a~b)
summary(res)


# repear 10,000
a<-runif(8000,0.96,0.99)
b<-runif(2000,0.92,0.96)
c<-c(a,b)
d<-runif(8000,0.96,0.99)
e<-runif(2000,0.92,0.96)
f<-c(d,e)
t<-cbind("Logit",c)
q<-cbind("Log",f)
t<-as.data.frame(t)
q<-as.data.frame(q)
colnames(t)<-c("Function","Number")
colnames(q)<-c("Function","Number")
w<-rbind(t,q)
w$Number<-as.numeric(w$Number)
ggplot(w, aes(x = Function, y = Number, fill = "#F8766D")) +
  geom_boxplot(alpha = 0.80) +
  theme(text = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none")
dat2$uorder<-10
dat2$uorder[1:42730]<-1
dat2$uorder[42731:84740]<-2
dat2$uorder[84741:127110]<-3
dat2$uorder[127111:169480]<-4
dat2$uorder[169481:211850]<-5
dat2$uorder[211851:254220]<-6
dat2$uorder[254221:296590]<-7
dat2$uorder[296591:338960]<-8
dat2$uorder[338961:381330]<-9
PRS<-read.table("Lipids/UKBB_182_ex_LP_30.profile",header=TRUE,as.is=TRUE,sep=" ")
dat<-PRS %>% select(1,6,3)
dat$CC<-dat$PHENO-1
dat<-dat[,-3]
colnames(dat)<-c("sample","NRA","CC")
dat2<-arrange(dat,NRA)
t<-fread("Lipids/UKBB_182_ex_LP_30.all.score")
t2<-t[match(dat2$sample,t$FID),]
dat2$wgrs<-t2$`5e-08`
exp(sum(dat2$wgrs)/sum(dat2$NRA))
t2<-caddat[match(dat2$sample,caddat$sample),]
t2$NRA<-dat2$NRA
t2$wgrs<-dat2$wgrs
t2$uorder<-dat2$uorder
t2$worder<-dat2$worder

#nLP
t1<-dat0[dat0$uorder==10,]
#t1<-dat0[dat0$worder==10,]
t1<-dat1[dat1$uorder==10,]
#t1<-dat1[dat1$worder==10,]
t2<-arrange(t1,lipidsNRA)
t2<-arrange(t1,NRA)

quant<-quantile(t2$NRA,probs = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
t2$uorder<-10
t2[t2$NRA<=quant[1],5]<-1
t2[quant[1]<t2$NRA & t2$NRA<=quant[2],5]<-2
t2[quant[2]<t2$NRA & t2$NRA<=quant[3],5]<-3
t2[quant[3]<t2$NRA & t2$NRA<=quant[4],5]<-4
t2[quant[4]<t2$NRA & t2$NRA<=quant[5],5]<-5
t2[quant[5]<t2$NRA & t2$NRA<=quant[6],5]<-6
t2[quant[6]<t2$NRA & t2$NRA<=quant[7],5]<-7
t2[quant[7]<t2$NRA & t2$NRA<=quant[8],5]<-8
t2[quant[8]<t2$NRA & t2$NRA<=quant[9],5]<-9
dat2l<-split(t2,t2$uorder)
t<-t2[t2$uorder!=1 & t2$uorder!=10,]
ures<-data.frame("nra"=0,"sd"=0,"p"=0,"ppa"=0,"N"=0,"nc"=0,"wgrs"=0,"sdw"=0)
for (i in 1:10){
  t<-dat2l[[i]]
  nra<-mean(t$NRA) # number of risk alleles
  sd<-sd(t$NRA) # sd of nRA
  N<-nrow(t) # Sample N
  nc<-sum(t$CC)
  p<-nc/N # prevalence
  ppa<-p/nra # prevalence per allele
  wgrs<-mean(t$lipidsNRA) # weighted GRS
  sdw<-sd(t$lipidsNRA) # sd of wGRS
  ures[i,]<-c(nra,sd,p,ppa,N,nc,wgrs,sdw)
}

quant<-quantile(t2$lipidsNRA,probs = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
t2$uorder<-10
t2[t2$lipidsNRA<=quant[1],5]<-1
t2[quant[1]<t2$lipidsNRA & t2$lipidsNRA<=quant[2],5]<-2
t2[quant[2]<t2$lipidsNRA & t2$lipidsNRA<=quant[3],5]<-3
t2[quant[3]<t2$lipidsNRA & t2$lipidsNRA<=quant[4],5]<-4
t2[quant[4]<t2$lipidsNRA & t2$lipidsNRA<=quant[5],5]<-5
t2[quant[5]<t2$lipidsNRA & t2$lipidsNRA<=quant[6],5]<-6
t2[quant[6]<t2$lipidsNRA & t2$lipidsNRA<=quant[7],5]<-7
t2[quant[7]<t2$lipidsNRA & t2$lipidsNRA<=quant[8],5]<-8
t2[quant[8]<t2$lipidsNRA & t2$lipidsNRA<=quant[9],5]<-9

t2_1<-t2[1:42370,]
t2_2<-t2[42371:84740,]
mean(t2_1$lipidsNRA)
mean(t2_2$lipidsNRA)
sd(t2_1$lipidsNRA)
sd(t2_2$lipidsNRA)
sum(t2_1$CC)/nrow(t2_1)
sum(t2_2$CC)/nrow(t2_2)
t2_1$newo<-1
t2_2$newo<-2
t3<-rbind(t2_1,t2_2)
t3$newo<-as.character(t3$newo)
fit <- glm(CC ~ 0+newo,data=t3,family="binomial")
# results like you have
results = coefficients(summary(fit))
# rename second column, SE for plotting
colnames(results)[2] = "SE"
results<-as.data.frame(results)
results$Estimate.ad<-results$Estimate+abs(results$Estimate[1])
results$OR<-exp(results$Estimate.ad)
results$CI.L<-exp(results$Estimate.ad-1.96*results$SE)       
results$CI.U<-exp(results$Estimate.ad+1.96*results$SE)
results

library(ggplot2)
#60,90,120,150,180,210,240
library("tidyverse")
# read uGRS
dat<-PRS %>% select(1,6,3)
dat$CC<-dat$PHENO-1
dat<-dat[,-3]
colnames(dat)<-c("sample","NRA","CC")
modComp  <- summary( glm(CC~NRA,family="binomial",data=dat) )$coefficients
dat$grp  <- round(dat$NRA) #cut(dat$NRA,breaks = seq(min(dat$NRA),max(dat$NRA),by=2),include.lowest=TRUE)
prev_p   <- aggregate(dat$CC,by=list(dat$grp),FUN=mean,na.rm=T)
prev_n   <- aggregate(dat$CC,by=list(dat$grp),FUN=sum,na.rm=T)
prev_N   <- aggregate(dat$CC,by=list(dat$grp),FUN=function(x) length(!is.na(x)))
allele_c <- aggregate(dat$NRA,by=list(dat$grp),FUN=mean)
d        <- cbind.data.frame(N=prev_N$x,n=prev_n$x,p=prev_p$x,numAllele=allele_c$x)
d$p_se   <- sqrt(d$p*(1-d$p)/d$N)
d        <- d[which(d$N>=200),]
d$px     <- d$p / d$numAllele
d$px_se  <- d$p_se / d$numAllele
modelLog    <- glm(cbind(n,N)~numAllele,data=d,family = binomial(link="log"))
predLog     <- predict(modelLog,type = "response")

PRS<-read.table("UKBB_182_final.profile",header=TRUE,as.is=TRUE,sep=" ")
modelLog180<-modelLog
x<-seq(40,320,1)
y<-predict(modelLog180,list(numAllele=x),type="response")
SNP<-c(180)
test0<-as.data.frame(cbind(x,y,SNP))
test<-test0
ggplot(test)+geom_line(aes(x=x,y=y,color=SNP,group=SNP))+xlab("Mean number of risk alleles")+ylab("Prevalence (%)")+theme_classic()

PRS<-read.table("../random/UKBB_91B.profile",header=TRUE,as.is=TRUE,sep=" ")
modelLog90<-modelLog
y<-predict(modelLog90,list(numAllele=x),type="response")
SNP<-c(90)
test0<-as.data.frame(cbind(x,y,SNP))
test<-rbind(test,test0)

PRS<-read.table("Lipids/UKBB_182_ex_LP.profile",header=TRUE,as.is=TRUE,sep=" ")
modelLog150<-modelLog
y<-predict(modelLog150,list(numAllele=x),type="response")
SNP<-c(150)
test0<-as.data.frame(cbind(x,y,SNP))
test<-rbind(test,test0)

PRS<-read.table("Estimate/test1.1.profile",header=TRUE,as.is=TRUE,sep=" ")
modelLog60<-modelLog
y<-predict(modelLog60,list(numAllele=x),type="response")
SNP<-c(60)
test0<-as.data.frame(cbind(x,y,SNP))
test<-rbind(test,test0)

PRS<-read.table("Estimate/test5.1.profile",header=TRUE,as.is=TRUE,sep=" ")
modelLog120<-modelLog
y<-predict(modelLog120,list(numAllele=x),type="response")
SNP<-c(120)
test0<-as.data.frame(cbind(x,y,SNP))
test<-rbind(test,test0)

modelLog210<-modelLog
y<-predict(modelLog210,list(numAllele=x),type="response")
SNP<-c(210)
test0<-as.data.frame(cbind(x,y,SNP))
test<-rbind(test,test0)

modelLog240<-modelLog
y<-predict(modelLog240,list(numAllele=x),type="response")
SNP<-c(240)
test0<-as.data.frame(cbind(x,y,SNP))
test<-rbind(test,test0)
test$SNP<-as.character(test$SNP)
test$SNP<-factor(test$SNP,levels=c(60,90,120,150,180,210,240))
ggplot(test)+geom_line(aes(x=x,y=y,color=SNP,group=SNP),size=2)+xlab("Mean number of risk alleles")+ylab("Prevalence (%)")+theme_classic()+ylim(0,0.5)+scale_x_continuous(breaks = c(100,150,200,250))+ labs(color = "Number of risk SNPs")

test0<-as.data.frame(cbind(d$numAllele,predLog))
test0$SNP<-"real"
colnames(test0)<-c("x","y","SNP")
test<-rbind(test,test0)
#638x430
pcad<-ggplot(test)+geom_line(aes(x=x,y=y,color=SNP,group=SNP),size=2)+xlab("Mean number of risk alleles")+ylab("Prevalence (%)")+theme_classic()+ylim(0,1)+scale_x_continuous(breaks = c(0,50,100,150,200,250,300))+ggtitle("CAD")
pbc<-ggplot(test)+geom_line(aes(x=x,y=y,color=SNP,group=SNP),size=2)+xlab("Mean number of risk alleles")+ylab("Prevalence (%)")+theme_classic()+ylim(0,1)+scale_x_continuous(breaks = c(0,50,100,150,200,250,300))+ggtitle("BC")
ppc<-ggplot(test)+geom_line(aes(x=x,y=y,color=SNP,group=SNP),size=2)+xlab("Mean number of risk alleles")+ylab("Prevalence (%)")+theme_classic()+ylim(0,1)+scale_x_continuous(breaks = c(0,50,100,150,200,250,300))+ggtitle("PC")
ggplot(test1)+geom_line(aes(x=x,y=y,color=Trait,group=group,linetype=SNP),size=1.5)+xlab("Mean number of risk alleles")+ylab("Prevalence (%)")+theme_classic()+ylim(0,1)+scale_x_continuous(breaks = c(0,50,100,150,200,250,300))+
  scale_linetype_manual(values=c("dotted", "solid"))


utils::data(hg19)
chr<-cad12[1,1]
start<-16902069
end<-19481547
colnames(hg19) = c("GENE_ID", "CHR", "TRX_START", "TRX_END", "LENGTH", "GENE_NAME", "TYPE")
gene_sub = hg19[hg19$CHR == 7,]
gene_sub = subset(gene_sub, gene_sub$TRX_START > (16902069-50000))
gene_sub = subset(gene_sub, gene_sub$TRX_END < (19481547+50000))
gene_sub = gene_sub[gene_sub$TYPE == "protein_coding",]
gene_sub = gene_sub[,c(3,4,6)]
gene_sub = reshape2::melt(gene_sub,id.vars = "GENE_NAME")
gene_sub$y_value = as.numeric(as.factor(gene_sub$GENE_NAME))
plot_lab = subset(gene_sub, gene_sub$variable == "TRX_END")

c = ggplot2::ggplot(gene_sub, ggplot2::aes_string(x = "value", y = "y_value")) +
  ggplot2::geom_line(ggplot2::aes_string(group = "GENE_NAME"), size = 2) + ggplot2::bw() +
  ggplot2::geom_text(data = plot_lab, ggplot2::aes_string(x = "value", y = "y_value", label = "GENE_NAME"),
                     hjust = -0.1,vjust = 0.3, size = 2.5) +
  ggplot2::theme(axis.title.y = ggplot2::element_text(color = "white", size = 28),
                 axis.text.y = ggplot2::element_blank(),
                 axis.ticks.y = ggplot2::element_blank()) + ggplot2::xlab(paste0("Chromosome ", chr, " Position")) +
  ggplot2::coord_cartesian(xlim = c(start,end), ylim = c(0,(max(gene_sub$y_value)+1)))

b=ggplot(gwas,aes(x=POS,y=LOG10P,color=disease))+geom_point(alpha=0.5)+xlab(paste(chr,"Position",sep="_")) + ylab("-log10(p-value)")+theme_minimal()
ggpubr::ggarrange(b, c, heights = c(3,1), nrow = 2, ncol = 1,
                  common.legend = TRUE, legend = "right")

cad_gwas_f_ld = RACER::ldRACER(assoc_data = cad_gwas_f, rs_col = 3, pops = "EUR", lead_snp = "rs2107595")
is_gwas_f_ld = RACER::ldRACER(assoc_data = is_gwas_f, rs_col = 3, pops = "EUR", lead_snp = "rs2107595")
pad_gwas_f_ld = RACER::ldRACER(assoc_data = pad_gwas_f, rs_col = 3, pops = "EUR", lead_snp = "rs2107595")

