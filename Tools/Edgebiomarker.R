### Construct Edge markers from nodes with known information;
args<-commandArgs(TRUE)
file1<-args[1]
file2<-args[2]

db1<-read.table(file1,header=TRUE,row.names=1)
db2<-read.table(file2,header=TRUE,row.names=1)

db1_mean<-c()
db2_mean<-c()
db1_std<-c()
db2_std<-c()

db1<-as.matrix(db1)
db2<-as.matrix(db2)

m<-nrow(db1)
p<-m+1
n<-nrow(db2)
q<-n+1
k1<-sqrt((ncol(db1)-1)/ncol(db1))
k2<-sqrt((ncol(db2)-1)/ncol(db2))

for (k in 1:m){
  db1_mean[k]<-mean(db1[k,])
  db1_std[k]<-sd(db1[k,])
  k=k+1
}
for (k in 1:n){
  db2_mean[k]<-mean(db2[k,])
  db2_std[k]<-sd(db2[k,])
  k=k+1
}

db11<-as.data.frame(db1)
db21<-as.data.frame(db2)

list<-read.table("KinaseSub.screen.uniq.txt")
for (i in 1:nrow(list)){
    a<-which(rownames(db1)==list[i,1])
    b<-which(rownames(db1)==list[i,2])
    for (j in 1:ncol(db1)){
         db11[p,j]<-(db1[a,j]-db1_mean[a])/(k1*db1_std[a])*(db1[b,j]-db1_mean[b])/(k1*db1_std[b])
         db11[p+1,j]<-(db1[a,j]-db2_mean[a])/(k2*db2_std[a])*(db1[b,j]-db2_mean[b])/(k2*db2_std[b])
         rownames(db11)[p]<-paste(rownames(db1)[a],rownames(db1)[b],"pos",sep="_")
         rownames(db11)[p+1]<-paste(rownames(db1)[a],rownames(db1)[b],"neg",sep="_")
         j=j+1
    }
    p=p+2
    
    for (j in 1:ncol(db2)){
         db21[q,j]<-(db2[a,j]-db1_mean[a])/(k1*db1_std[a])*(db2[b,j]-db1_mean[b])/(k1*db1_std[b])
         db21[q+1,j]<-(db2[a,j]-db2_mean[a])/(k2*db2_std[a])*(db2[b,j]-db2_mean[b])/(k2*db2_std[b])
         rownames(db21)[q]<-paste(rownames(db1)[a],rownames(db1)[b],"pos",sep="_")
         rownames(db21)[q+1]<-paste(rownames(db1)[a],rownames(db1)[b],"neg",sep="_")
         j=j+1
    }
    q=q+2

}

write.table(db11,file=paste(file1,".edge.txt",sep=""),quote=FALSE,sep="\t")
write.table(db21,file=paste(file2,".edge.txt",sep=""),quote=FALSE,sep="\t")


dbnew<-cbind(db11,db21)
n1<-ncol(db11)
n2<-ncol(dbnew)
pvalue<-apply(dbnew,1,function(x){t.test(x[1:n1],x[(n1+1):n2])$p.value})
dbnew1<-dbnew[pvalue<=0.05,]

write.table(dbnew1,file=paste(file1,".feature.sig.txt",sep=""),quote=FALSE,sep="\t")


