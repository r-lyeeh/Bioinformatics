quant<-quantile(caddat2$NRA,probs=c(0.5,0.6,0.7,0.8,0.9))

ggplot()+geom_density(data=test,aes(x=numAllele,y=..density..*427306),fill="grey",alpha=0.5)+xlab("Number of risk alleles")+ylab("Number of samples")+theme_classic()+geom_vline(xintercept=192.3440,linetype = "dashed")+geom_vline(xintercept=194.6810,linetype = "dashed")+geom_vline(xintercept=197.3250,linetype = "dashed")+geom_vline(xintercept=201.0845 ,linetype = "dashed")+geom_line(data=d,aes(x=numAllele,y=predLogit*150000))+scale_y_continuous(name="Number of samples",sec.axis=sec_axis(~.*2.848707/427306,name="Prevalence"),expand = c(0, 0))
d$predLog<-predLog

ggplot()+xlab("Number of risk alleles")+ylab("CAD prevalence (%)")+geom_line(data=d,aes(x=numAllele,y=predLogit),color="seagreen",linetype="dashed",size=1)+geom_line(data=d,aes(x=numAllele,y=predLog),color="coral1",linetype="dashed",size=1)+geom_density(data=test,aes(x=numAllele,y=..density..*2.848707),fill="grey",alpha=0.5,color="white")+theme_classic()+geom_vline(xintercept=179.3010,linetype = "dashed")+geom_vline(xintercept=201.0845 ,linetype = "dashed")+geom_errorbar(data=d,aes(x=numAllele,ymax = p + p_se, ymin = p - p_se),color="grey")+xlim(155,225)+geom_point(data=d,aes(x=numAllele,y=p))


t<-caddat2[caddat2$uorder!=1 & caddat2$uorder!=10,]
t1<-t[t$sex==1,]
t2<-t1[t1$smoke==1,]
t3<-t2[t2$Diabetes==1,]
t3<-t3[is.na(t3$NRA)==FALSE,]
sum(t3$CC)/nrow(t3)

t3<-t2[t2$Diabetes==0,]
t3<-t3[is.na(t3$NRA)==FALSE,]
sum(t3$CC)/nrow(t3)

t2<-t1[t1$smoke==0,]
t3<-t2[t2$Diabetes==0,]
t3<-t3[is.na(t3$NRA)==FALSE,]
sum(t3$CC)/nrow(t3)
