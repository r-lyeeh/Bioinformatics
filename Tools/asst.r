setwd("/home/pang/data/SCP_involve/MAFF_Moritz_20190927")
case<-read.table("case_maff.frq",sep='\t',header=FALSE,stringsAsFactors = FALSE)
control<-read.table("control_maff.frq",sep='\t',header=FALSE,stringsAsFactors = FALSE)
a<-matrix(1:4,nrow=2,ncol=2)
pvalue<-case[,1:2]
colnames(pvalue)<-c("CHRM","POS","chisq","fisher")
pvalue<-pvalue[-1,]

for (i in 2:47){
  a[1,]<-c(as.numeric(as.character(case$V4[i]))*as.numeric(strsplit(as.character(case$V6[i]),":")[[1]][2]),as.numeric(as.character(control$V4[i]))*as.numeric(strsplit(as.character(control$V6[i]),":")[[1]][2]))
  a[2,]<-c(as.numeric(as.character(case$V4[i]))*as.numeric(strsplit(as.character(case$V5[i]),":")[[1]][2]),as.numeric(as.character(control$V4[i]))*as.numeric(strsplit(as.character(control$V5[i]),":")[[1]][2]))
  pvalue[i-1,3]<-chisq.test(a)$p.value
  pvalue[i-1,4]<-fisher.test(a)$p.value
}
write.table(pvalue,"maff_asso_test.txt",sep="\t",col.names = TRUE)
