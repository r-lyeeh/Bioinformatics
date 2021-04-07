setwd("Desktop/MAFF")
wes<-read.table("WES_CAD_sample.txt", sep=" ")
pde5a<-read.table("PDE5A.txt",sep=" ")
#control:905 cad:79
maff<-read.table("temp3.txt",sep=" ")
maff<-maff[,-1]
colnames(maff)<-c("V1","V2")
#control:266 cad:21
control0<-wes[wes$V2=="CONTROL",]
cad0<-wes[wes$V2!="CONTROL",]
#sampling t=3500
t1<-3500
t2<-1500
control1<-control0[sample(nrow(control0),t1,replace = FALSE),]
cad1<-cad0[sample(nrow(cad0),t2,replace = FALSE),]
#compare maff
maff1_1<-intersect(maff, control1)  
maff1_2<-intersect(maff, cad1)  
chim<-matrix(1:4,nrow=2,ncol=2)
chim[1,]<-c(nrow(maff1_1),nrow(maff1_2))
n0<-(t1-nrow(maff1_1))
n1<-(t2-nrow(maff1_2))
chim[2,]<-c(n0,n1)
chisq.test(chim)
#compare pde5a
pde5a1_1<-intersect(pde5a, control1)  
pde5a1_2<-intersect(pde5a, cad1)  
chip<-matrix(1:4,nrow=2,ncol=2)
chip[1,]<-c(nrow(pde5a1_1),nrow(pde5a1_2))
m0<-(t1-nrow(pde5a1_1))
m1<-(t2-nrow(pde5a1_2))
chip[2,]<-c(m0,m1)
chisq.test(chip)

#control:266 cad:21
control0<-wes[wes$V2=="CONTROL",]
cad0<-wes[wes$V2=="HARD",]
#sampling t=3500
t1<-1500
t2<-1500
control1<-control0[sample(nrow(control0),t1,replace = FALSE),]
cad1<-cad0[sample(nrow(cad0),t2,replace = FALSE),]
#compare maff
maff1_1<-intersect(maff, control1)  
maff1_2<-intersect(maff, cad1)  
chim<-matrix(1:4,nrow=2,ncol=2)
chim[1,]<-c(nrow(maff1_1),nrow(maff1_2))
n0<-(t1-nrow(maff1_1))
n1<-(t2-nrow(maff1_2))
chim[2,]<-c(n0,n1)
chisq.test(chim)
#compare pde5a
pde5a1_1<-intersect(pde5a, control1)  
pde5a1_2<-intersect(pde5a, cad1)  
chip<-matrix(1:4,nrow=2,ncol=2)
chip[1,]<-c(nrow(pde5a1_1),nrow(pde5a1_2))
m0<-(t1-nrow(pde5a1_1))
m1<-(t2-nrow(pde5a1_2))
chip[2,]<-c(m0,m1)
chisq.test(chip)

