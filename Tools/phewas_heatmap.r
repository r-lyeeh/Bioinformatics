setwd("/home/pang/Desktop/Scattered/PheWAS/plink/disease")
phewas<-read.table("ukb_disease_gene_final",sep=" ",header = TRUE,as.is = TRUE)
phewas$logor <- log(phewas$OR)
#phewas$logorl95 <- log(phewas$L95)
#phewas$diff <- phewas$logor - phewas$logorl95
#phewas$seor <- phewas$diff / 1.96
#phewas$z2<-phewas$logor/phewas$seor
phewas$z <- phewas$logor / phewas$SE
#phewas$p2<-(1 - pnorm(abs(phewas$z))) * 2
#phewas$z2<-qnorm(phewas$P)
phewas2<-phewas[phewas$OR!=0,]
phewas2<-phewas2[,-6]
library(tidyr)
phewas3<-unite(phewas2,SNP,c("SNPID","Gene"),remove=TRUE)
phewas3<-phewas3[phewas3$SNP!="NA_NA",]
ggplot(phewas3,aes(x=SNP,y=Disease,fill=z))+geom_tile(colour="gray80") +theme_bw(10) + xlab("SNP") + ylab("Disease") + scale_fill_gradient2(low = "darkred", mid = "white", high = "darkgreen",midpoint = 0, space = "Lab", na.value = "grey10", guide = "colourbar", limits=c(-10, 3))+theme(axis.text.x = element_text(angle = 90, hjust = 1))

library(reshape2)
t4<-dcast(phewas3,formula=Disease~SNP,value.var ="z",fill=0)
rownames(t4)<-t4[,1]
t4<-t4[,-1]
library(tidyverse)
#rowcolor<-unique( phewas3[ , 2:3 ] )
#rowcolor2<-rowcolor %>% arrange(Disease)
rowcolor2<-read.table("heatmap_color.txt",header = FALSE,as.is = TRUE)
rowcolor2<-rowcolor2[,-3]
colnames(rowcolor2)<-c("Disease_group","Disease","Dcolor")
HC0<-dcast(phewas3,formula=Disease~SNP,value.var ="P",fill=1)
highlightCell0<-as.data.frame(which(HC0<0.05,arr.ind=TRUE))
highlightCell0$color<-"black"
highlightCell0$lwd<-2
library(heatmap3)
pdf("test.pdf")
heatmap3(phema,Rowv=NA,RowSideColors = rowcolor2$Dcolor,
         RowSideLabs=c("Disease_group"),highlightCell = highlightCell0,showRowDendro=FALSE,showColDendro = FALSE)
#heatmap3(phema,Rowv=NA,RowSideColors = rowcolor2$Dcolor,
RowSideLabs=c("Disease_group"),highlightCell = highlightCell0,showRowDendro=FALSE,showColDendro = FALSE,legendfun=function() showLegend(legend=c("Cancer","Cardiovascular_diseases","Digestive_diseases","Endocrine_disorders","Eye_diseases","Genito-urinary_diseases","Musculoskeletal_disease","Neurological","Respiratory_diseases"),col=c("deepskyblue","darkolivegreen1","darkorange","firebrick1","darkorchid1","darksalmon","darkseagreen1","gold","darkslategray1"),cex=1))
dev.off()
