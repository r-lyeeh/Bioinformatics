##plot using R
########R codes
expIBD<-function(workdirectory,file){
  setwd(workdirectory)
  setwd("2_sample")
  genome = read.table(paste(file,"_merged_QC.genome",sep=""),h=T)
  fam = read.table(paste("../1_SNP/",file,"_merged_QC.fam",sep=""),h=F)
  head(genome)
  library(igraph)
  
  g = graph.data.frame(genome[,c(2,4)],directed = F)
  
  g.dgr = degree(g)
  
  rm.node = c()
  while(max(g.dgr)>1){
    node = names(g.dgr)[which.max(g.dgr)]
    rm.node = c(rm.node,node)
    g = delete.vertices(g,v = node)
    g.dgr = degree(g)
  }
  
  g.sub = get.data.frame(g)
  dim(g.sub)
  
  rm.node = union(rm.node, g.sub[,2])
  
  rm.node = fam[match(rm.node,fam$V2),]
  write.table(rm.node,file=paste(file,"_samples_excluded_by_IBD.txt",sep=""),row.names = F,col.names = F,quote = F,sep = "\t")
}
args<-commandArgs(TRUE)
expIBD(args[1],args[2])