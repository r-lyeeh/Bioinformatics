##plot using R
########R codes
mdsplot<-function(workdirectory,file){
  setwd(workdirectory)
  setwd("2_sample/population_stratification")
  d = read.table("plink.mds",h=T)
  fam = read.table(paste("../",file,"_merged_QC_no_relative.fam",sep=""),header = F)
  ind = match(fam$V2,d$IID)
  d = d[ind,]
  
  
  ##fam$v6,control=1,case=2
  ##fam$v5,male=1,female=2
  d$pop = as.factor(ifelse(fam$V6==2,"case","control"))
  d$sex = as.factor(ifelse(fam$V5==2,"female","male"))
  ##case&female=red
  d$col=NA
  d$col[intersect(which(fam$V6==2),which(fam$V5==2))]="red"
  ##case&male="orange"
  d$col[intersect(which(fam$V6==2),which(fam$V5==1))]="orange"
  ##control&femal="blue"
  d$col[intersect(which(fam$V6==1),which(fam$V5==2))]="blue"
  ##control&male="green"
  d$col[intersect(which(fam$V6==1),which(fam$V5==1))]="green"
  ##control&femal="yellow"
  d$col[intersect(which(fam$V6==-9),which(fam$V5==2))]="yellow"
  ##control&male="pink"
  d$col[intersect(which(fam$V6==-9),which(fam$V5==1))]="pink"
  
  table(d$col)
  ###plot 1-6 pcs
  for(i in 1:5){
    for(j in (i+1):6){
      pdf(paste(paste("PC",i,j,sep="-"),".pdf",sep = ""))
      plot(d[,i+3],d[,j+3],xlab=colnames(d)[i+3],ylab = colnames(d)[j+3],col=d$col,pch=19)
      dev.off()
    }}
# remove samples
# d[d$C2<(mean(d$C2)-3*sd(d$C2)),]
# d[d$C2>(mean(d$C2)+3*sd(d$C2)),]
  # Calculate Outliers To Be Removed based on the Largest Two Principle Component
  pca_res<-d
  c1<-as.numeric(pca_res$C1)
  c2<-as.numeric(pca_res$C2)
  w_c1<-c(which(c1>5*sd(c1)+mean(c1)),which(c1<mean(c1)-5*sd(c1)))
  w_c2<-c(which(c2>5*sd(c2)+mean(c2)),which(c2<mean(c2)-5*sd(c2)))
  w_rm<-unique(c(w_c1,w_c2))
  indi_to_be_rm_by_pca<-pca_res[w_rm,1:2]
  write.table(indi_to_be_rm_by_pca,file=paste(file,"_merged_QC_no_relative_prune.mds.sample",sep=""),sep="\t",quote=F,col.names=F,row.names=F)
}

args<-commandArgs(TRUE)
mdsplot(args[1],args[2])