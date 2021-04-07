#setwd("Desktop/GRN/")
setwd("/home/pang/Desktop/MP_Johan_20200221/All9")
hsq_full<-read.table("hsq_res_snpmod_full_20200102",header = TRUE,as.is = FALSE,sep=" ")
hsq_noncad<-read.table("hsq_res_snpmod_nonsnpcad304ld02_20200102",header=TRUE,as.is = FALSE,sep=" ")

getpvalue<-function(hsq_full){
  dfname<-deparse(substitute(hsq_full))
  hsq_full$h2<-hsq_full$Variance_VgVgL/0.4*100
  hsq_full$ld02<-hsq_full$h2/hsq_full$nSnp_Final_ld02
  hsq_noncad$h2<-hsq_noncad$Variance_VgVgL/0.4*100
  hsq_noncad$ld02<-hsq_noncad$h2/hsq_noncad$nSnp_Final_ld02
  
  stnor<-function(x){x/sum(x)}
  hsq_full_s<-hsq_full
  hsq_full_s$Variance_VgVgL<-stnor(hsq_full$Variance_VgVgL)
  hsq_full_s$Se_VgVgL<-stnor(hsq_full$Se_VgVgL)
  hsq_full_s$h2<-stnor(hsq_full$h2)
  hsq_full_s$ld02<-stnor(hsq_full$ld02)
  #14,15,17,19
  
  hsq_full_p<-as.data.frame(matrix(0,nrow(hsq_full),nrow(hsq_full)))
  rownames(hsq_full_p)<-hsq_full[,1]
  colnames(hsq_full_p)<-hsq_full[,1]
  hsq_full_pld02<-hsq_full_p
  for (i in (1:(nrow(hsq_full_s)-1))){
    for (j in ((i+1):nrow(hsq_full_s))){
      x<-c(hsq_full_s$Se_VgVgL[i],hsq_full_s$h2[i])
      y<-c(hsq_full_s$Se_VgVgL[j],hsq_full_s$h2[j])
      x2<-c(hsq_full_s$Se_VgVgL[i],hsq_full_s$ld02[i])
      y2<-c(hsq_full_s$Se_VgVgL[j],hsq_full_s$ld02[j])
      p<-t.test(x,y,alternative = "two.sided",var.equal =TRUE)$p.value
      p2<-t.test(x2,y2,alternative = "two.sided",var.equal =TRUE)$p.value
      #z<-c(hsq_full$h2[i],hsq_full$h2[j])
      #pe<-c(hsq_full$nSnp_Final_ld02[i],hsq_full$nSnp_Final_ld02[j])/sum(hsq_full$nSnp_Final_ld02[i]+hsq_full$nSnp_Final_ld02[j])
      hsq_full_p[i,j]<-p
      hsq_full_pld02[i,j]<-p2
    }
  }
  write.table(hsq_full_p,paste(dfname,"_p.txt",sep=""),row.names = TRUE, col.names = TRUE,quote = FALSE,sep="\t")
  write.table(hsq_full_pld02,paste(dfname,"_pld02.txt",sep=""),row.names = TRUE, col.names = TRUE,quote = FALSE,sep="\t")
}

getpvalue(hsq_full)
getpvalue(hsq_noncad)

ttest<-function(hsq_full){ dfname<-deparse(substitute(hsq_full))
hsq_full$h2<-hsq_full$Variance_VgVgL/0.4*100
hsq_full$ld02<-hsq_full$h2/hsq_full$nSnp_Final_ld02
hsq_noncad$h2<-hsq_noncad$Variance_VgVgL/0.4*100
hsq_noncad$ld02<-hsq_noncad$h2/hsq_noncad$nSnp_Final_ld02

stnor<-function(x){x/sum(x)}
hsq_full_s<-hsq_full
#hsq_full_s$Variance_VgVgL<-stnor(hsq_full$Variance_VgVgL)
#hsq_full_s$Se_VgVgL<-stnor(hsq_full$Se_VgVgL)
#hsq_full_s$h2<-stnor(hsq_full$h2)
#hsq_full_s$ld02<-stnor(hsq_full$ld02)
#14,15,17,19

hsq_full_p<-as.data.frame(matrix(0,nrow(hsq_full),nrow(hsq_full)))
rownames(hsq_full_p)<-hsq_full[,1]
colnames(hsq_full_p)<-hsq_full[,1]
hsq_full_pld02<-hsq_full_p
for (i in (1:(nrow(hsq_full_s)-1))){
  for (j in ((i+1):nrow(hsq_full_s))){
    x<-c(hsq_full_s$Se_VgVgL[i],hsq_full_s$h2[i])
    y<-c(hsq_full_s$Se_VgVgL[j],hsq_full_s$h2[j])
    x2<-c(hsq_full_s$Se_VgVgL[i],hsq_full_s$ld02[i])
    y2<-c(hsq_full_s$Se_VgVgL[j],hsq_full_s$ld02[j])
    p<-t.test(x,y,alternative = "two.sided",var.equal =TRUE)$p.value
    p2<-t.test(x2,y2,alternative = "two.sided",var.equal =TRUE)$p.value
    #z<-c(hsq_full$h2[i],hsq_full$h2[j])
    #pe<-c(hsq_full$nSnp_Final_ld02[i],hsq_full$nSnp_Final_ld02[j])/sum(hsq_full$nSnp_Final_ld02[i]+hsq_full$nSnp_Final_ld02[j])
    hsq_full_p[i,j]<-p
    hsq_full_pld02[i,j]<-p2
  }
}
write.table(hsq_full_p,paste(dfname,"_p.ttest",sep=""),row.names = TRUE, col.names = TRUE,quote = FALSE,sep="\t")
write.table(hsq_full_pld02,paste(dfname,"_pld02.ttest",sep=""),row.names = TRUE, col.names = TRUE,quote = FALSE,sep="\t")
}
#ttest(hsq_full)
#ttest(hsq_noncad)

