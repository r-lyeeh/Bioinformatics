G6_dis<-read.table("G6.PRS_sum.all.score",header = T,as.is=T)
G6_sample<-read.table("G6.cad",header = T,as.is = T)
G6_dis$cad<-G6_sample$CAD
m0<-format(round(mean(G6_dis[G6_dis$cad==0,4]), 2), nsmall = 2)
m1<-format(round(mean(G6_dis[G6_dis$cad==1,4]), 2), nsmall = 2)
pdf("G6_pl.pdf")
ggplot(G6_dis,aes(x=G6_dis$X5.005e.05, fill=as.factor(G6_dis$cad))) +geom_histogram(binwidth=.03, alpha=.5, position="identity")+geom_vline(xintercept = as.numeric(m0), colour="red",linetype="dashed",alpha=0.5)+geom_text(aes(as.numeric(m0)-0.3,50,label = m0), color="red",size=5,alpha=0.5)+geom_vline(xintercept = as.numeric(m1), colour="steelblue",linetype="dashed",alpha=0.5)+geom_text(aes(as.numeric(m1)+0.3,50,label = m1), color="steelblue",size=5,alpha=0.5)+theme_bw()+labs(title="G6_wPRS",x="PRS",fill= "CAD")+scale_fill_manual(labels=c("Control","Case"), values = c("red", "steelblue"))
dev.off()
