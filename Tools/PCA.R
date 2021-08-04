cosmic<-read.csv("/Users/magica/Desktop/AGEING/CTD_genes_diseases.csv",header = TRUE)
immune<-cosmic[grep(pattern="lymphoid",cosmic$cosmic70),]
write.csv(immune,"/Users/magica/Desktop/HumanGenome/annovar/immune.csv")

coad<-read.csv("/Users/magica/Desktop/HNSL/pancancer/failtrmp/rsem.coad.norm.tsv",header = TRUE,sep='\t')
lusc<-read.csv("/Users/magica/Desktop/HNSL/pancancer/failtrmp/rsem.lusc.norm.tsv",header = TRUE,sep='\t')
kirp<-read.csv("/Users/magica/Desktop/HNSL/pancancer/failtrmp/rsem.kirp.norm.tsv",header = TRUE,sep='\t')
stad<-read.csv("/Users/magica/Desktop/HNSL/pancancer/tempraw/rsem.stad.norm.tsv",header = TRUE,sep='\t')
brca<-read.csv("/Users/magica/Desktop/HNSL/pancancer/tempraw/rsem.brca.norm.tsv",header = TRUE,sep='\t')
kipan<-read.csv("/Users/magica/Desktop/HNSL/pancancer/tempraw/rsem.kipan.norm.tsv",header = TRUE,sep='\t')
pca=function(data,group=1,title = NULL,legend="right",text=FALSE,xlim=NULL,ylim=NULL,size.point=3,size.text=2){
  group=as.factor(group)
  if (length(group)==1) legend.position = "none" else legend.position = legend
  
  # calculation
  pr=prcomp(data,cor=T)  # correlation matrix
  Comp=pr$sdev
  Comp.var=Comp*Comp
  
  Comp.1=paste("Comp.1(",round(Comp.var[1]/sum(Comp.var)*100,2),"%)",sep="") # The first PC (x axis)
  Comp.2=paste("Comp.2(",round(Comp.var[2]/sum(Comp.var)*100,2),"%)",sep="") # The second PC (y axis)
  
  pcomp=pr$x # result matrix; rows with samples, columns with principles;
  
  #plot
  label=as.character(rownames(pcomp))
  require(ggplot2)
  p=ggplot(as.data.frame(pcomp),aes(PC1,PC2))
  coord = coord_fixed(ratio =Comp.var[2]/Comp.var[1], xlim = xlim, ylim = ylim)
  theme1=theme_bw() + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  geom = geom_point(size=size.point,shape=16,aes(colour=group))
  theme2=theme(legend.position=legend.position)+theme(plot.title = element_text(size = 20))+theme(axis.title = element_text(size = 12))
  lab=labs(x=Comp.1,y=Comp.2,title=title,size=10)
  if (text==FALSE) plot=p+coord+theme1+geom+theme2+lab
  else  {point.text= geom_text(label=label,size=size.text,hjust=1,vjust=1);
  plot=p+coord+theme1+geom+theme2+lab+point.text}
  
  plot		
}
apply(stad,1,su)
sum(log(stad)