setwd("Desktop/GRN/")
hsq_full<-read.table("hsq_full.p",header = TRUE,as.is = FALSE,sep="\t")
hsq_noncad<-read.table("hsq_noncad.p",header=TRUE,as.is = FALSE,sep="\t")

stnor<-function(x){x/sum(x)}
hsq_full_s<-hsq_full
hsq_full_s[,14]<-stnor(hsq_full[,14])
hsq_full_s[,15]<-stnor(hsq_full[,15])
hsq_full_s[,17]<-stnor(hsq_full[,17])
hsq_full_s[,19]<-stnor(hsq_full[,19])

hsq_full_p<-as.data.frame(matrix(0,nrow(hsq_full),nrow(hsq_full)))
rownames(hsq_full_p)<-hsq_full[,1]
colnames(hsq_full_p)<-hsq_full[,1]
hsq_full_pld02<-hsq_full_p
for (i in (1:(nrow(hsq_full_s)-1))){
  for (j in ((i+1):nrow(hsq_full_s))){
    #x<-c(hsq_full_s[i,15],hsq_full_s[i,17])
    #y<-c(hsq_full_s[j,15],hsq_full_s[j,17])
    #x2<-c(hsq_full_s[i,15],hsq_full_s[i,19])
    #y2<-c(hsq_full_s[j,15],hsq_full_s[j,19])
    #p<-t.test(x,y,alternative = "two.sided",var.equal =TRUE)$p.value
    #p2<-t.test(x2,y2,alternative = "two.sided",var.equal =TRUE)$p.value
    z<-c(hsq_full[i,17],hsq_full[j,17])
    pe<-c(hsq_full[i,13],hsq_full[j,13])/sum(hsq_full[i,13]+hsq_full[j,13])
    hsq_full_p[i,j]<-p
    hsq_full_pld02[i,j]<-p2
  }
}
write.table(hsq_full_p,"hsq_full_p.txt",row.names = TRUE, col.names = TRUE,quote = FALSE,sep="\t")
write.table(hsq_full_pld02,"hsq_full_pld02.txt",row.names = TRUE, col.names = TRUE,quote = FALSE,sep="\t")

hsq_noncad_s<-hsq_noncad
hsq_noncad_s[,14]<-stnor(hsq_noncad[,14])
hsq_noncad_s[,15]<-stnor(hsq_noncad[,15])
hsq_noncad_s[,17]<-stnor(hsq_noncad[,17])
hsq_noncad_s[,19]<-stnor(hsq_noncad[,19])

hsq_noncad_p<-as.data.frame(matrix(0,nrow(hsq_noncad),nrow(hsq_noncad)))
rownames(hsq_noncad_p)<-hsq_noncad[,1]
colnames(hsq_noncad_p)<-hsq_noncad[,1]
hsq_noncad_pld02<-hsq_noncad_p
for (i in (1:(nrow(hsq_noncad_s)-1))){
  for (j in ((i+1):nrow(hsq_noncad_s))){
    x<-c(hsq_noncad_s[i,14],hsq_noncad_s[i,15],hsq_noncad_s[i,17])
    y<-c(hsq_noncad_s[j,14],hsq_noncad_s[j,15],hsq_noncad_s[j,17])
    x2<-c(hsq_noncad_s[i,14],hsq_noncad_s[i,15],hsq_noncad_s[i,19])
    y2<-c(hsq_noncad_s[j,14],hsq_noncad_s[j,15],hsq_noncad_s[j,19])
    p<-t.test(x,y,alternative = "two.sided",var.equal =TRUE)$p.value
    p2<-t.test(x2,y2,alternative = "two.sided",var.equal =TRUE)$p.value
    hsq_noncad_p[i,j]<-p
    hsq_noncad_pld02[i,j]<-p2
  }
}
write.table(hsq_noncad_p,"hsq_noncad_p.txt",row.names = TRUE, col.names = TRUE,quote = FALSE,sep="\t")
write.table(hsq_noncad_pld02,"hsq_noncad_pld02.txt",row.names = TRUE, col.names = TRUE,quote = FALSE,sep="\t")
