library(ggplot2)
library(reshape2)
library(lme4)
library(boot)
library(gridExtra)
library(ggthemes)

###############健康自评、年龄、居住方式、孤独感(SAH为纵坐标）/加入：风玫瑰图################################

###############导入数据库，更改SRH数据类型#######################

old<-read.table("D:\\tupian2.txt",header=T)
#fix(old)
old$SRH<-as.character(old$SRH)

###############健康自评、年龄################################

tmp1<-melt(old[,c("SRH","Age")],id.vars="SRH")
plot1<-ggplot(tmp1,aes(x=SRH,y=value))+geom_jitter(alpha = .1,col=2)+geom_violin(alpha=.35)+scale_y_sqrt()+labs(y="年龄(岁)",x="自评健康")
plot1<-plot1+annotate("rect", xmin=c(0.5,2.5), xmax=c(2.5,5.5), ymin=57, ymax=Inf, fill=c(3,"blue"), alpha=0.1)+annotate("text", x=c(1.5,4), y=58, label =c("较差", "良好"),  size=4)

###############健康自评、居住方式(SAH为纵坐标)################################

old1<-read.table("D:\\srh1.txt",header=T)
old1$SRH<-as.character(old1$SRH)

plot2<-ggplot(old1)+geom_bar(aes(x=live,fill=SRH, shape=SRH, colour=SRH))+coord_polar()+labs(x="居住方式",y="比例(%)")

###############健康自评、孤独感(SAH为纵坐标)################################

old2<-read.table("D:\\srh2.txt",header=T)
old2$SRH<-as.character(old2$SRH)
old2$loneliness<-as.character(old2$loneliness)

plot3<-ggplot(old2)+geom_bar(aes(x=loneliness,fill=SRH, shape=SRH, colour=SRH))+coord_polar()+labs(x="孤独感",y="比例(%)")

###############居住方式、孤独感(SAH为纵坐标)################################

old3<-read.table("D:\\srh3.txt",header=T)
old3$Loneliness<-as.character(old3$Loneliness)
old3$live<-as.character(old3$live)

tmp4<-melt(old3[,c("Loneliness","live")], id.vars="live")
ggplot(tmp4, aes(factor(live), y = value, fill=factor(live)))+geom_boxplot()+facet_wrap(~variable, scales="free_y")
plot4<-ggplot(old3,aes(x=live,y=Loneliness))+stat_sum(aes(size=..n..,group = 1))+scale_size_area(max_size=10)+labs(x="居住方式",y="孤独感")


grid.arrange(plot1, plot2, plot3, plot4, ncol=2)
