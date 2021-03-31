for (i in 1:22)
  {
  print(paste("Processing chromosome ",i,"...",sep=""))
  data <- read.table(paste("chr",i,".alignments.snp.strand",sep=""), header=T, row.names=NULL)
  if(length(data[which(duplicated(data$pos)),]$pos>0)) data <- data[-which(duplicated(data$pos)),]
  data <- data[which(data$type %in% "strand"),]
  out <- data[3]
  write.table(out, paste("flip_chr",i,".txt", sep=""), quote=F, row.names=F, col.names=F)  
  }
