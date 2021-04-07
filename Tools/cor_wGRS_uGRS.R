#Liability model 
#Correlation test for uGRS and wGRS - Figure 1C
#Load packages
library("ggpubr")
#Read data
cor<-read.table("cor_uGRS_wGRS.txt",sep="\t",header = TRUE)
colnames(cor)[5]<-"uGRS"
#Correlation test: Pearson, Spearman, Kendall
#Result check:1) covariation linear 2)normal distribution - shapiro.test()/ggpubr::ggqqplot() 3)shapiro-wilk test
#shapiro.test(cor$wGRS) 3-5000 sample size
ggqqplot(cor$wGRS,ylab="wGRS")
ggqqplot(cor$uGRS,ylab="uGRS")
#Plot
p1<-ggscatter(cor, x = "wGRS", y = "uGRS", 
              add = "reg.line", conf.int = TRUE, 
              cor.coef = TRUE, cor.method = "pearson",
              xlab = "wGRS", ylab = "uGRS",alpha=0.25,add.params = list(color = "blue", fill = "lightgray"))
#R=0.89, p<2.2e-16
res<-cor.test(cor$wGRS,cor$uGRS,method="pearson")
#t = 343.73, df = 29790, p-value < 2.2e-16,cor=0.8936658, 0.95CI=0.8913558~0.8959294
#res<-cor.test(cor$wGRS,cor$uGRS,method="kendall")
#kendallL z = 181.89, p-value < 2.2e-16, tau=0.7028687
#qplot(cor$wGRS,cor$uGRS, data = cor, geom = c("point", "smooth"), method = "lm", alpha = I(1 / 5))
#Similarities test
library(vegan)
