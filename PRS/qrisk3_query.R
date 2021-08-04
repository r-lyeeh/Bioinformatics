#backup
test1$qriskclass<-1
test1$qriskclass[test1$qrisk3>3 & test1$qrisk3<=6]<-2
test1$qriskclass[test1$qrisk3>6 & test1$qrisk3<=9]<-3
test1$qriskclass[test1$qrisk3>9 & test1$qrisk3<=12]<-4
test1$qriskclass[test1$qrisk3>12 & test1$qrisk3<=15]<-5
test1$qriskclass[test1$qrisk3>15]<-6
table(test1$qriskclass)
analysis0<-function(test1,pheno,qclass){
  q1<-test1[test1$qriskclass==qclass,]
  q1<-q1[is.na(q1$NRA)==FALSE,]
  quant<-quantile(q1$NRA,probs=c(0.2,0.8,0.95),na.rm = TRUE)
  q1.1<-q1[q1$NRA<=quant[1],]
  c1<-pqstat(q1.1)
  q1.2<-q1[q1$NRA>quant[1] & q1$NRA<=quant[2],]
  c2<-pqstat(q1.2)
  q1.3<-q1[q1$NRA>quant[2] & q1$NRA<=quant[3],]
  c3<-pqstat(q1.3)
  q1.4<-q1[q1$NRA>quant[3],]
  c4<-pqstat(q1.4)
  
  d <- as.data.frame(rbind(c1,c2,c3,c4))
  colnames(d)<-c("N","n","p","numAllele")
  d$p_se   <- sqrt(d$p*(1-d$p)/d$N)
  modelLog    <- glm(cbind(n,N)~numAllele,data=d,family = binomial(link="log"))
  predLog     <- predict(modelLog,type = "response")
  Rlog        <- round(unlist(cor.test(d$p,predLog)[c("estimate","conf.int")]),2)
  
  Cols <- c("seagreen")
  l    <- which(d[,"p_se"]>0)
  
  plot(d[l,"numAllele"],d[l,"p"],pch=19,axes=FALSE,ylim=c(0,max(d$p+d$p_se)),xlim=c(min(d$numAllele-sd(d$numAllele)),max(d$numAllele+sd(d$numAllele))),main=pheno,xlab="Mean number of risk alleles",ylab="Prevalence");axis(1);axis(2)
  for(i in 1:nrow(d)){
    if(d[i,"p_se"]>0){
      arrows(d[i,"numAllele"],d[i,"p"]-d[i,"p_se"],d[i,"numAllele"],d[i,"p"]+d[i,"p_se"],col="grey",angle=90,len=0.05)
      arrows(d[i,"numAllele"],d[i,"p"]+d[i,"p_se"],d[i,"numAllele"],d[i,"p"]-d[i,"p_se"],col="grey",angle=90,len=0.05)
    }
  }
  matlines(d$numAllele,cbind(predLog),col=Cols,lwd=2)
  quant<-quantile(q1$NRA,probs=c(0.1,0.9))
  polygon(x=c(quant[1],quant[1],quant[2],quant[2]),y=c(0,max(d$p),max(d$p),0),col=rgb(0.2,0.2,0.2,0.2),border = NA)
  legend(x=min(d$numAllele),y=max(d$p+2*d$p_se),
         legend=c(paste0("Log: R=",Rlog[1]," (95%CI: [",Rlog[2],"-",Rlog[3],"])")),
         col=Cols,cex=0.9,box.lty=0,lwd=2,bg="transparent")
  
}
#par(mfrow=c(2,3))
#analysis0(test1,"UKBB CAD Qrisk <3%",1)
#analysis0(test1,"UKBB CAD Qrisk 3%-6%",2)
#analysis0(test1,"UKBB CAD Qrisk 6%-9%",3)
#analysis0(test1,"UKBB CAD Qrisk 9%-12%",4)
#analysis0(test1,"UKBB CAD Qrisk 12%-15%",5)
#analysis0(test1,"UKBB CAD Qrisk >15%",6)


plotan0<-function(test1,score){
  analysis0<-function(test1,class){
    q1<-filter(test1,scoreclass==class)
    #a1<-filter(test1,test1[,colnames(test1) %in% score]==class)
    #q1<-test1 %>% filter_at(vars(score), any_vars(. == class))
    
    q1<-q1[is.na(q1$NRA)==FALSE,]
    quant<-quantile(q1$NRA,probs=c(0.2,0.8,0.95),na.rm = TRUE)
    q1.1<-q1[q1$NRA<=quant[1],]
    c1<-pqstat(q1.1)
    q1.2<-q1[q1$NRA>quant[1] & q1$NRA<=quant[2],]
    c2<-pqstat(q1.2)
    q1.3<-q1[q1$NRA>quant[2] & q1$NRA<=quant[3],]
    c3<-pqstat(q1.3)
    q1.4<-q1[q1$NRA>quant[3],]
    c4<-pqstat(q1.4)
    
    d <- as.data.frame(rbind(c1,c2,c3,c4))
    colnames(d)<-c("N","n","p","numAllele")
    d$p_se   <- sqrt(d$p*(1-d$p)/d$N)
    modelLog    <- glm(p~numAllele,data=d,family = binomial(link="log"))
    predLog     <- predict(modelLog,type = "response")
    Rlog        <- round(unlist(cor.test(d$p,predLog)[c("estimate","conf.int")]),2)
    d$predLog<-predLog
    d$Rlog<-Rlog[[1]]
    d$Rlogb<-Rlog[[2]]
    d$Rlogu<-Rlog[[3]]
    d$class<-class
    return(d)
  }
  d1<-analysis0(test1,1)
  d2<-analysis0(test1,2)
  d3<-analysis0(test1,3)
  d4<-analysis0(test1,4)
  d<-rbind(d1,d2,d3,d4)
  d$class<-as.character(d$class)
  quant<-quantile(test1$NRA,probs=c(0.1,0.9))
  pol<-data.frame(x=c(quant[1],quant[1],quant[2],quant[2]),y=c(min(d$p),max(d$p),max(d$p),min(d$p)))
  #spline.d <- as.data.frame(spline(d$numAllele, d$predLog))
  p<-ggplot()+geom_polygon(data=pol, mapping=aes(x=x, y=y),fill="grey",alpha=0.4)+geom_point(data=d,aes(x=numAllele,y=p,color=class),pch=19)+geom_smooth(data=d,aes(x=numAllele,y=p,group=class,color=class),method = "glm",formula=y~log10(x),se=FALSE)+geom_errorbar(data=d,aes(x=numAllele,ymin=p-p_se, ymax=p+p_se), size=0.3,width=1)+theme_classic()+ylab("Prevalence")+xlab("Mean number of risk alleles")+ scale_color_manual(values=c("#33ccff","#66cc33","#ffcc00","#ff3300"),labels=c(paste0(score," <5%: Log: R=",d$Rlog[1]," (95%CI: [",d$Rlogb[1],"-",d$Rlogu[1],"])"),paste0(score," 5-10%: Log: R=",d$Rlog[5]," (95%CI: [",d$Rlogb[5],"-",d$Rlogu[5],"])"),paste0(score," 10-15%: Log: R=",d$Rlog[9]," (95%CI: [",d$Rlogb[9],"-",d$Rlogu[9],"])"),paste0(score," >15%: Log: R=",d$Rlog[13]," (95%CI: [",d$Rlogb[13],"-",d$Rlogu[13],"])")))+ theme(legend.title =element_blank(),legend.position="top")+guides(color=guide_legend(nrow=4))
  return(p)
}


plotan1<-function(test1,score){
  analysis0<-function(test1,class){
    q1<-filter(test1,scoreclass==class)
    #a1<-filter(test1,test1[,colnames(test1) %in% score]==class)
    #q1<-test1 %>% filter_at(vars(score), any_vars(. == class))
    
    q1<-q1[is.na(q1$NRA)==FALSE,]
    quant<-quantile(q1$NRA,probs=c(0.2,0.8,0.95),na.rm = TRUE)
    q1.1<-q1[q1$NRA<=quant[1],]
    c1<-pqstat(q1.1)
    q1.2<-q1[q1$NRA>quant[1] & q1$NRA<=quant[2],]
    c2<-pqstat(q1.2)
    q1.3<-q1[q1$NRA>quant[2] & q1$NRA<=quant[3],]
    c3<-pqstat(q1.3)
    q1.4<-q1[q1$NRA>quant[3],]
    c4<-pqstat(q1.4)
    
    d <- as.data.frame(rbind(c1,c2,c3,c4))
    colnames(d)<-c("N","n","p","numAllele")
    d$p_se   <- sqrt(d$p*(1-d$p)/d$N)
    modelLog    <- glm(cbind(n,N)~numAllele,data=d,family = binomial(link="log"))
    predLog     <- predict(modelLog,type = "response")
    Rlog        <- round(unlist(cor.test(d$p,predLog)[c("estimate","conf.int")]),2)
    d$predLog<-predLog
    d$Rlog<-Rlog[[1]]
    d$Rlogb<-Rlog[[2]]
    d$Rlogu<-Rlog[[3]]
    d$class<-class
    return(d)
  }
  d1<-analysis0(test1,1)
  d2<-analysis0(test1,2)
  d3<-analysis0(test1,3)
  d4<-analysis0(test1,4)
  d<-rbind(d1,d2,d3,d4)
  d$class<-as.character(d$class)
  quant<-quantile(test1$NRA,probs=c(0.1,0.9))
  pol<-data.frame(x=c(quant[1],quant[1],quant[2],quant[2]),y=c(min(d$p),max(d$p),max(d$p),min(d$p)))
  #spline.d <- as.data.frame(spline(d$numAllele, d$predLog))
  p<-ggplot()+geom_polygon(data=pol, mapping=aes(x=x, y=y),fill="grey",alpha=0.4)+geom_point(data=d,aes(x=numAllele,y=p,color=class),pch=19)+geom_line(data=d,aes(x=numAllele,y=predLog,group=class,color=class))+geom_errorbar(data=d,aes(x=numAllele,ymin=p-p_se, ymax=p+p_se), size=0.3,width=1)+theme_classic()+ylab("Prevalence")+xlab("Mean number of risk alleles")+ scale_color_manual(values=c("#33ccff","#66cc33","#ffcc00","#ff3300"),labels=c(paste0(score," <5%: Log: R=",d$Rlog[1]," (95%CI: [",d$Rlogb[1],"-",d$Rlogu[1],"])"),paste0(score," 5-10%: Log: R=",d$Rlog[5]," (95%CI: [",d$Rlogb[5],"-",d$Rlogu[5],"])"),paste0(score," 10-15%: Log: R=",d$Rlog[9]," (95%CI: [",d$Rlogb[9],"-",d$Rlogu[9],"])"),paste0(score," >15%: Log: R=",d$Rlog[13]," (95%CI: [",d$Rlogb[13],"-",d$Rlogu[13],"])")))+ theme(legend.title =element_blank(),legend.position="top")+guides(color=guide_legend(nrow=4))
  
  return(p)
}

plotan2<-function(test1,score){
  analysis0<-function(test1,class){
    q1<-filter(test1,scoreclass==class)
    #a1<-filter(test1,test1[,colnames(test1) %in% score]==class)
    #q1<-test1 %>% filter_at(vars(score), any_vars(. == class))
    
    q1<-q1%>%mutate(quantile=ntile(NRA,10))
    d<-data.frame("N"=0,"n"=0,"p"=0,"numAllele"=0)
    for (i in 1:10){
      qi<-q1[q1$quantile==i,]
      d[i,]<-pqstat(qi)
    }
    
    d$p_se   <- sqrt(d$p*(1-d$p)/d$N)
    modelLog    <- glm(p~numAllele,data=d,family = binomial(link="log"))
    predLog     <- predict(modelLog,type = "response")
    Rlog        <- round(unlist(cor.test(d$p,predLog)[c("estimate","conf.int")]),2)
    d$predLog<-predLog
    d$Rlog<-Rlog[[1]]
    d$Rlogb<-Rlog[[2]]
    d$Rlogu<-Rlog[[3]]
    d$class<-class
    return(d)
  }
  d1<-analysis0(test1,1)
  d2<-analysis0(test1,2)
  d3<-analysis0(test1,3)
  d4<-analysis0(test1,4)
  d<-rbind(d1,d2,d3,d4)
  d$class<-as.character(d$class)
  quant<-quantile(test1$NRA,probs=c(0.1,0.9))
  pol<-data.frame(x=c(quant[1],quant[1],quant[2],quant[2]),y=c(min(d$p),max(d$p),max(d$p),min(d$p)))
  #spline.d <- as.data.frame(spline(d$numAllele, d$predLog))
  p<-ggplot()+geom_polygon(data=pol, mapping=aes(x=x, y=y),fill="grey",alpha=0.4)+geom_point(data=d,aes(x=numAllele,y=p,color=class),pch=19)+geom_smooth(data=d,aes(x=numAllele,y=p,group=class,color=class),method = "glm",formula=y~log10(x),se=FALSE)+geom_errorbar(data=d,aes(x=numAllele,ymin=p-p_se, ymax=p+p_se), size=0.3,width=1)+theme_classic()+ylab("Prevalence")+xlab("Mean number of risk alleles")+ scale_color_manual(values=c("#33ccff","#66cc33","#ffcc00","#ff3300"),labels=c(paste0(score," <5%: Log: R=",d$Rlog[1]," (95%CI: [",d$Rlogb[1],"-",d$Rlogu[1],"])"),paste0(score," 5-10%: Log: R=",d$Rlog[5]," (95%CI: [",d$Rlogb[5],"-",d$Rlogu[5],"])"),paste0(score," 10-15%: Log: R=",d$Rlog[9]," (95%CI: [",d$Rlogb[9],"-",d$Rlogu[9],"])"),paste0(score," >15%: Log: R=",d$Rlog[13]," (95%CI: [",d$Rlogb[13],"-",d$Rlogu[13],"])")))+ theme(legend.title =element_blank(),legend.position="top")+guides(color=guide_legend(nrow=4))
  return(p)
}

plotan3<-function(test1,score){
  analysis0<-function(test1,class){
    q1<-filter(test1,scoreclass==class)
    #a1<-filter(test1,test1[,colnames(test1) %in% score]==class)
    #q1<-test1 %>% filter_at(vars(score), any_vars(. == class))
    
    q1<-q1%>%mutate(quantile=ntile(NRA,10))
    d<-data.frame("N"=0,"n"=0,"p"=0,"numAllele"=0)
    for (i in 1:10){
      qi<-q1[q1$quantile==i,]
      d[i,]<-pqstat(qi)
    }
    
    d$p_se   <- sqrt(d$p*(1-d$p)/d$N)
    modelLog    <- glm(p~numAllele,data=d,family = binomial(link="log"))
    predLog     <- predict(modelLog,type = "response")
    Rlog        <- round(unlist(cor.test(d$p,predLog)[c("estimate","conf.int")]),2)
    d$predLog<-predLog
    d$Rlog<-Rlog[[1]]
    d$Rlogb<-Rlog[[2]]
    d$Rlogu<-Rlog[[3]]
    d$class<-class
    return(d)
  }
  d1<-analysis0(test1,1)
  d2<-analysis0(test1,2)
  d3<-analysis0(test1,3)
  d4<-analysis0(test1,4)
  d<-rbind(d1,d2,d3,d4)
  d$class<-as.character(d$class)
  quant<-quantile(test1$NRA,probs=c(0.1,0.9))
  pol<-data.frame(x=c(quant[1],quant[1],quant[2],quant[2]),y=c(min(d$p),max(d$p),max(d$p),min(d$p)))
  #spline.d <- as.data.frame(spline(d$numAllele, d$predLog))
  p<-ggplot()+geom_polygon(data=pol, mapping=aes(x=x, y=y),fill="grey",alpha=0.4)+geom_point(data=d,aes(x=numAllele,y=p,color=class),pch=19)+geom_line(data=d,aes(x=numAllele,y=predLog,group=class,color=class))+geom_errorbar(data=d,aes(x=numAllele,ymin=p-p_se, ymax=p+p_se), size=0.3,width=1)+theme_classic()+ylab("Prevalence")+xlab("Mean number of risk alleles")+ scale_color_manual(values=c("#33ccff","#66cc33","#ffcc00","#ff3300"),labels=c(paste0(score," <5%: Log: R=",d$Rlog[1]," (95%CI: [",d$Rlogb[1],"-",d$Rlogu[1],"])"),paste0(score," 5-10%: Log: R=",d$Rlog[5]," (95%CI: [",d$Rlogb[5],"-",d$Rlogu[5],"])"),paste0(score," 10-15%: Log: R=",d$Rlog[9]," (95%CI: [",d$Rlogb[9],"-",d$Rlogu[9],"])"),paste0(score," >15%: Log: R=",d$Rlog[13]," (95%CI: [",d$Rlogb[13],"-",d$Rlogu[13],"])")))+ theme(legend.title =element_blank(),legend.position="top")+guides(color=guide_legend(nrow=4))
  return(p)
}

plotan4<-function(test1,score){
  analysis0<-function(test1,class){
    q1<-filter(test1,scoreclass==class)
    #a1<-filter(test1,test1[,colnames(test1) %in% score]==class)
    #q1<-test1 %>% filter_at(vars(score), any_vars(. == class))
    q1$CC<-q1$CC-1
    q1$grp  <- round(q1$NRA) #cut(dat$NRA,breaks = seq(min(dat$NRA),max(dat$NRA),by=2),include.lowest=TRUE)
    prev_p   <- aggregate(q1$CC,by=list(q1$grp),FUN=mean,na.rm=T)
    prev_n   <- aggregate(q1$CC,by=list(q1$grp),FUN=sum,na.rm=T)
    prev_N   <- aggregate(q1$CC,by=list(q1$grp),FUN=function(x) length(!is.na(x)))
    allele_c <- aggregate(q1$NRA,by=list(q1$grp),FUN=mean)
    d        <- cbind.data.frame(N=prev_N$x,n=prev_n$x,p=prev_p$x,numAllele=allele_c$x)
    d        <- d[which(d$N>=200),]
    d$p_se   <- sqrt(d$p*(1-d$p)/d$N)
    modelLog    <- glm(p~numAllele,data=d,family = binomial(link="log"))
    predLog     <- predict(modelLog,type = "response")
    Rlog        <- round(unlist(cor.test(d$p,predLog)[c("estimate","conf.int")]),2)
    d$predLog<-predLog
    d$Rlog<-Rlog[[1]]
    d$Rlogb<-Rlog[[2]]
    d$Rlogu<-Rlog[[3]]
    d$class<-class
    return(d)
  }
  d1<-analysis0(test1,1)
  d2<-analysis0(test1,2)
  d3<-analysis0(test1,3)
  d4<-analysis0(test1,4)
  d<-rbind(d1,d2,d3,d4)
  d$class<-as.character(d$class)
  quant<-quantile(test1$NRA,probs=c(0.1,0.9))
  pol<-data.frame(x=c(quant[1],quant[1],quant[2],quant[2]),y=c(min(d$p),max(d$p),max(d$p),min(d$p)))
  #spline.d <- as.data.frame(spline(d$numAllele, d$predLog))
  p<-ggplot()+geom_polygon(data=pol, mapping=aes(x=x, y=y),fill="grey",alpha=0.4)+geom_point(data=d,aes(x=numAllele,y=p,color=class),pch=19)+geom_line(data=d,aes(x=numAllele,y=predLog,group=class,color=class))+geom_errorbar(data=d,aes(x=numAllele,ymin=p-p_se, ymax=p+p_se), size=0.3,width=1)+theme_classic()+ylab("Prevalence")+xlab("Mean number of risk alleles")+ scale_color_manual(values=c("#33ccff","#66cc33","#ffcc00","#ff3300"),labels=c(paste0(score," <5%: Log: R=",d$Rlog[1]," (95%CI: [",d$Rlogb[1],"-",d$Rlogu[1],"])"),paste0(score," 5-10%: Log: R=",d$Rlog[5]," (95%CI: [",d$Rlogb[5],"-",d$Rlogu[5],"])"),paste0(score," 10-15%: Log: R=",d$Rlog[9]," (95%CI: [",d$Rlogb[9],"-",d$Rlogu[9],"])"),paste0(score," >15%: Log: R=",d$Rlog[13]," (95%CI: [",d$Rlogb[13],"-",d$Rlogu[13],"])")))+ theme(legend.title =element_blank(),legend.position="top")+guides(color=guide_legend(nrow=4))
  return(p)
}

t2<-caddat2[match(Fabian$sample,caddat2$sample),]
test0<-Fabian
test0$NRA<-t2$NRA
test1<-test0[is.na(test0$NRA)==FALSE,]
test1<-test1[is.na(test1$qrisk3_adj)==FALSE,]
test1$qrisk3<-test1$qrisk3_adj
#test1<-test1[is.na(test1$qrisk3)==FALSE,]
#test1$qrisk3<-test1$qrisk32
test1<-test1[is.na(test1$ASSIGN_SIMD)==FALSE,]
test1<-test1[is.na(test1$FRS)==FALSE,]
score<-"qrisk3"
score<-"FRS"
score<-"ASSIGN_SIMD"
#1.) determine in UKBB 10-year risk based on the ASCVD risk calculator (please ask Fabian Starnecker for help) √
#2.) pool all people with an estimated risk of <3% — 3-6% — 6-9% — 9-12%  — 12-15% —  >15%  as indicated by the the ASCVD risk calculator (we may have to regroup at some point)
test1$scoreclass<-1
test1<-test1%>%mutate(scoreclass=replace(scoreclass,!!as.symbol(score)>5&!!as.symbol(score)<=10,2))
test1<-test1%>%mutate(scoreclass=replace(scoreclass,!!as.symbol(score)>10&!!as.symbol(score)<=15,3))
test1<-test1%>%mutate(scoreclass=replace(scoreclass,!!as.symbol(score)>15,4))
table(test1$scoreclass)

#test1<-test1%>%mutate(quantile=ntile(NRA,10))
#rownames(test1)[]<-"scoreclass"

#3.) plot the logarithmic curve of CAD risk (based on the number of risk alleles and percentile of the weighted score) for each of the groups
pqstat<-function(q1.1){
  r1<-sum(q1.1$CC)/nrow(q1.1)-1
  N1<-nrow(q1.1)
  n1<-sum(q1.1$CC)-nrow(q1.1)
  alle1<-mean(q1.1$NRA)
  c<-cbind(N1,n1,r1,alle1)
  return(c)
}

#4.) determine risk for the 0-20% — 20-80% — 80-95% — 95-100% percentiles based on the PRS in each of the ASCVD risk groups
qriskprscvd<-function(test1,class){
  q1<-test1[test1$scoreclass==class,]
  q1<-q1[is.na(q1$NRA)==FALSE,]
  quant<-quantile(q1$NRA,probs=c(0.2,0.8,0.95),na.rm = TRUE)
  q1.1<-q1[q1$NRA<=quant[1],]
  r1<-sum(q1.1$CC)/nrow(q1.1)-1
  q1.2<-q1[q1$NRA>quant[1] & q1$NRA<=quant[2],]
  r2<-sum(q1.2$CC)/nrow(q1.2)-1
  q1.3<-q1[q1$NRA>quant[2] & q1$NRA<=quant[3],]
  r3<-sum(q1.3$CC)/nrow(q1.3)-1
  q1.4<-q1[q1$NRA>quant[3],]
  r4<-sum(q1.4$CC)/nrow(q1.4)-1
  r0<-c(r1,r2,r3,r4,nrow(q1))
  return(r0)
}
qpcres<-data.frame("0-20"=0,"20-80"=0,"80-95"=0,"95-100"=0,"sample_size"=0)
for (i in 1:6){
  qpcres[i,]<-qriskprscvd(test1,i)
}

test1$smoking[test1$smoking<=2]<-0
test1$smoking[test1$smoking>2]<-1
clp<-data.frame("age"=0,"sex"=0,"BMI"=0,"SBP"=0,"HDL"=0,"Chol"=0,"Smoking"=0,"BPT"=0,"FH"=0,"townsend"=0,"ethic"=0,"Dia1"=0,"Dia2"=0,"SLE"=0,"SMI"=0,"CKD"=0,"RA"=0,"MIG"=0,"ED"=0,"AF"=0,"aty"=0,"corti"=0)
for (i in 1:4){
  q1<-test1[test1$scoreclass==i,]
  clp[i,]<-cbind(mean(q1$age),mean(q1$sex),mean(q1$BMI,na.rm=TRUE),mean(q1$sbp),mean(q1$HDL,na.rm = TRUE),mean(q1$chol),mean(q1$smoking,na.rm=TRUE),mean(q1$bmpa,na.rm=TRUE),mean(q1$familyHeart,na.rm = TRUE),mean(q1$Townsend,na.rm=TRUE),mean((q1$ethic-1)/5,na.rm=TRUE),mean(q1$dia,na.rm=TRUE),mean(q1$dia2,na.rm=TRUE),mean(q1$SLE,na.rm=TRUE),mean(q1$SMI,na.rm=TRUE),mean(q1$CKD,na.rm=TRUE),mean(q1$RA,na.rm=TRUE),mean(q1$MIG,na.rm=TRUE),mean(q1$ED,na.rm=TRUE),mean(q1$AF,na.rm = TRUE),mean(q1$atypical,na.rm=TRUE),mean(q1$corticosteroids,na.rm=TRUE))
}  

  


#5.) determine the factor by which risk differs from the middle part of the distribution curve (20-80%), e.g. 0.5 for the 0-20 percentile or 1.8 for the 80-95 percentile
qpcrest<-qpcres
qpcrest$X0.20<-qpcrest$X0.20/qpcrest$X20.80
qpcrest$X80.95<-qpcrest$X80.95/qpcrest$X20.80
qpcrest$X95.100<-qpcrest$X95.100/qpcrest$X20.80
qpcrest$X20.80<-1
#6.) the question is, whether the factor by which risk is affected based on the PRS is the same across all ASCVD risk groups. My prediction is that the higher the ASCVD risk the higher is the factor in the 80-95% — 95-100% percentile groups.
qpcresv<-qpcres
qpcresv[2,]<-qpcresv[2,]/qpcresv[1,]
qpcresv[3,]<-qpcresv[3,]/qpcresv[1,]
qpcresv[4,]<-qpcresv[4,]/qpcresv[1,]
qpcresv[1,]<-1

pdf(paste0(score,"_01.pdf"),width=6.50,height=7.80)
plotan0(test1,score)
dev.off()
pdf(paste(score,"_02.pdf"),width=6.50,height=7.80)
plotan1(test1,score)
dev.off()
plotan2(test1,score)
pdf(paste(score,"_03.pdf"),width=6.50,height=7.80)
plotan3(test1,score)
dev.off()
plotan4(test1,score)
