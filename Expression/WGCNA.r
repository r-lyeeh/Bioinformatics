#library path
myPaths <- .libPaths()   # get the paths
myPaths <- c("/home/pang/program/miniconda3/lib/R/library", myPaths)  # switch them
myPaths <- c(myPaths[2],myPaths[3],myPaths[4])  # switch them
myPaths <- c(myPaths[2],myPaths[1],myPaths[4])  # switch them
.libPaths(myPaths)  # reassign them
# Expression test 
library(data.table)
library(tidyverse)
library(DESeq2)
library(purrr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggrepel)
library(ggeasy)
library(limma)
library(VennDiagram)
setwd("~/Desktop/Expression_20210620/data/GTExv8/expression/gene_counts/")
setwd("~/Desktop/Expression_20210620/data/GTExv8/expression/transcript_counts/")
# GTEx_v8 209 variants take into consideration
gtex<-fread("~/Desktop/Expression_20210620/data/GTExv8/genotype/241_common_nra_GTEX.profile")
starnet<-fread("~/Desktop/Expression/data/STARNET/241_common_starnet_nra.profile")
ukbb<-fread("~/Desktop/Expression/data/UKBB/241_common_nra_UKBB.profile")
gtexpheno<-fread("~/Desktop/Expression_20210620/data/GTExv8/phenotype/phs000424.v8.pht002742.v8.p2.c1.GTEx_Subject_Phenotypes.GRU.txt")
gtexp<-gtexpheno%>%select("phv00169062.v8.p2.c1", "phv00169063.v8.p2.c1","phv00169061.v8.p2.c1","phv00169163.v8.p2.c1","phv00169162.v8.p2.c1","phv00169164.v8.p2.c1")
colnames(gtexp)<-as.character(gtexp[1,])
gtexp<-gtexp[-1,]
t2<-gtexp[match(gtex$IID,gtexp$SUBJID),]
gtex<-cbind(gtex,t2[,c(1,2,5)])
#Phv00169162.v5.p1	MHHRTATT	Heart attack, acute myocardial infarction, acute coronary syndrome
#phv00169163.v5.p1	MHHRTDIS	Ischemic Heart Disease (coronary artery disease (CAD), coronary heart disease, ischemic cardiomyopathy)
#phv00169164.v5.p1	MHHRTDISB	Heart Disease
gtexall<-fread("~/Desktop/Expression_20210620/data/GTExv8/genotype/GTEx_all.all.score")
gtexallselect<-gtexall%>%select("FID","IID","5e-08","5.005e-05","0.00010005","0.00100005","0.01","0.0500001","0.1","0.2","0.5","1")
colnames(gtexallselect)<-c("FID","IID","PRS1","PRS2","PRS3","PRS4","PRS5","PRS6","PRS7","PRS8","PRS9","PRS10")
t2<-gtexallselect[match(gtex$IID,gtexallselect$IID),]
gtex<-cbind(gtex,t2[,3:12])
colnames(gtex)[6]<-"PRS0"

qures<-function(x){
        qures<-data.frame("0-20"=0,"20-80"=0,"80-95"=0,">95"=0,"0-90"=0,">90"=0)
        quant<-quantile(ukbb$SCORESUM,probs=c(0.2,0.8,0.9,0.95))
        t1<-x[x$SCORESUM<quant[1],]
        t2<-x[x$SCORESUM>=quant[1] & x$SCORESUM<quant[2],]
        t3<-x[x$SCORESUM>=quant[2] & x$SCORESUM<quant[4],]
        t4<-x[x$SCORESUM>=quant[4],]
        t5<-x[x$SCORESUM<quant[3],]
        t6<-x[x$SCORESUM>=quant[3],]
        qures[1,]<-c(mean(t1$SCORESUM),mean(t2$SCORESUM),mean(t3$SCORESUM),mean(t4$SCORESUM),mean(t5$SCORESUM),mean(t6$SCORESUM))
        qures[2,]<-c(nrow(t1),nrow(t2),nrow(t3),nrow(t4),nrow(t5),nrow(t6))
        return(qures)
}

ggplot(sumall,aes(x=SCORESUM)) +
        geom_histogram(aes(fill=data),binwidth = 1,color="white",size=0.2,position="identity") +
        facet_grid(data ~ ., scales = "free_y") +xlab("Number of risk alleles")+ylab("Number of samples")+labs(fill="Datasets")+theme_bw()+ theme(legend.position='none')+geom_vline(xintercept = 227.5690, linetype="longdash",color = "black",alpha=0.8)
pdf("~/Desktop/Expression_20210620/data/Summary.pdf")

#expression

topbot<-function(x){
        load(x)
        exp<-t(tss.gene.tpm)
        colnames(exp)<-exp[1,]
        exp<-exp[-1,]
        exp<-as.data.frame(exp)
        exp$sample<-rownames(exp)
        exp<-exp %>% separate(sample,sep="-",into=c("p1","p2","p3","p4","p5"))
        cols=c("p1","p2")
        exp$samplef<-do.call(paste,c(exp[cols],sep="-"))
        t2<-gtex[match(exp$samplef,gtex$IID),]
        qs<-cbind(nrow(t2[t2$quantile==1,]),nrow(t2[t2$quantile==10,]))
        return(qs)
}

#create a list of the files from your target directory
file_list <- list.files(path="~/Desktop/Expression/data/GTExv8/expression/gene_counts/",pattern="*.txt")
# remove 34
file_list<-file_list[-34]

#initiate a blank data frame, each iteration of the loop will append the data from the given file to this variable
res <- data.frame()

#had to specify columns to get rid of the total column
for (i in 1:length(file_list)){
        qs<-cbind(file_list[i],topbot(file_list[i]))
        res<-rbind(res,qs)
}


# Co-expression network
library('gplots')
library('ggplot2')
library('knitr')
library('limma')
library('reshape2')
library('RColorBrewer')
library('WGCNA')

# differentially expressed gene
# RPKM >= 0.5 and gene-level read counts >= 10
##### seperate read counts into sub-tissue #####
file_list <- list.files(path="~/Desktop/Expression/data/GTExv8/expression")
genecount<-fread("gene_counts/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct")
genecount<-fread("transcripts_counts/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_expected_count.gct")
for (i in 1:length(file_list)){
        load(file_list[i])
        file<-str_split(file_list[i],"_")[[1]][1]
        expt<-tss.gene.tpm
        testcol<-match(colnames(expt),colnames(genecount))
        exp2<-subset(genecount,select=testcol)
        fwrite(exp2,paste(file,"transcripts_counts.txt"),sep="\t",quote=FALSE)
}

#setwd("/raid3/Data_public/GTEx_phs000424.v8.p2/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm_by_tissue")
#samp.attr = fread("../Annotation/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
#rna.samp.attr = samp.attr %>% filter(SMAFRZE=='RNASEQ')
#dim(rna.samp.attr)

#mytissues = table(rna.samp.attr$SMTSD)
#alltss.tsp.tpm = fread("../GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_transcript_tpm.gct")
#alltss.tsp.tpm = as.data.frame(alltss.tsp.tpm)
#alltss.rnaseq.samp.summary = c()
#for(tss in names(mytissues)){
#  print(tss)
#  tss.samp = rna.samp.attr$SAMPID[rna.samp.attr$SMTSD == tss]
#  tss.tsp.tpm = alltss.tsp.tpm[,c(1,2, colnames(alltss.tsp.tpm) %in% tss.samp)]
#  nsamp = ncol(tss.tsp.tpm)-2
#  print(nsamp)
#  alltss.rnaseq.samp.summary = rbind.data.frame(alltss.rnaseq.samp.summary,
#                                                c(tissue=tss, nsamp = nsamp))
#  save(tss.tsp.tpm,file=paste(tss,"_transcript_tpm.rda",sep = ""))
#}

# extract decile 1 and decile 10
gtex1<-gtex%>%mutate(quantile=ntile(SCORESUM,10))
gtex0<-rbind(gtex1[gtex1$quantile==5,],gtex1[gtex1$quantile==10,])

#UKBB 0-90
gtex1<-gtex[gtex$SCORESUM<227.5690,]
gtex10<-gtex[gtex$SCORESUM>=227.5690,]
gtex1$quantile<-"0-90"
gtex10$quantile<-">=90"
gtex1$quantile<-1
gtex10$quantile<-10
gtex0<-rbind(gtex1,gtex10)


file_list <- list.files(path="~/Desktop/Expression_20210620/data/GTExv8/expression/gene_counts/",pattern="*.txt")
file_list <- list.files(path="~/Desktop/Expression_20210620/data/GTExv8/expression/transcripts_counts/")
#file_list<-file_list[-35]
siggene<-data.frame("file"=0,"siggene"=0)
for (i in 1:length(file_list)){
        exp<-fread(file_list[i])
        file<-str_split(file_list[i],"_")[[1]][1]
        g<-diffplt(exp,file)
        #g<-fread(paste0(file,"_diffexp_decile1_decile10.tab"))
        siggene[i,]<-cbind(file,nrow(g[g$DiffExpressed!="no",]))
}

diffplt<-function(exp,file){
        exp<-as.data.frame(exp)
        ann<-exp[,1:2]
        colnames(exp)[1]<-"Name"
        rownames(exp)<-exp$Name
        exp<-exp[,-c(1,2)]
        exp<-t(exp)
        exp<-as.data.frame(exp)
        exp$sample<-rownames(exp)
        exp<-exp %>% separate(sample,sep="-",into=c("p1","p2","p3","p4","p5"))
        cols=c("p1","p2")
        exp$samplef<-do.call(paste,c(exp[cols],sep="-"))
        t2<-gtex0[match(exp$samplef,gtex0$IID),]
        t2<-t2[is.na(t2$IID)==FALSE,]
        testcol<-match(t2$IID,exp$samplef)
        #exp2<-subset(exp,select=testcol)
        exp2<-exp[testcol,]
        testrow<-match(exp$samplef,gtex0$IID)
        samples<-gtex0[testrow,]
        samples<-samples[is.na(samples$IID)==FALSE,]
        samples<-samples%>%mutate(quantile=case_when(quantile==1~"decile1",quantile==10~"decile10"))
        samples$quantile<-as.factor(samples$quantile)
        samples<-as.data.frame(samples)
        #raw_counts<-exp2[,-c(56201:56206)]
        raw_counts<-exp2[,-c((ncol(exp2)-5):ncol(exp2))]
        raw_counts<-t(raw_counts)
        # 基因在每个样本中平均表达量为10就要被过滤
        low_count_mask <- rowSums(raw_counts) < 10*ncol(raw_counts)
        raw_counts <- raw_counts[!low_count_mask,]
        sprintf("Removing %d low-count genes (%d remaining).",           sum(low_count_mask), sum(!low_count_mask))
        # col sample; row gene
        ##筛选方差前25%的基因## trans --> left number extremely low thus no filter 
        ##Otherwise the ajusted p-value will be extremely large
        #m.vars=apply(raw_counts,1,var)
        #expro.upper=raw_counts[which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25),na.rm=TRUE)[4]),]
        #dim(expro.upper)
        #raw_counts<-raw_counts[match(rownames(expro.upper),rownames(raw_counts)),]
        
        #### 第一步，构建DESeq2的DESeq对象
        #exprSet<-raw_counts # gene
        exprSet<-as.data.frame(raw_counts) %>%mutate_all(as.integer) #transcripts
        group_list=samples$quantile
        colData <- data.frame(row.names=colnames(exprSet),group_list=group_list)
        dds <- DESeqDataSetFromMatrix(countData = exprSet,colData = colData,
                                      design = ~ group_list)
        #### 第二步，进行差异表达分析
        #f  <- 'data/Step03-DESeq2_dds2.Rdata'
        #if(!file.exists(f)){
        dds2 <- DESeq(dds)
        # 保存差异表达分析结果
        #  save(dds2, file = "data/Step03-DESeq2_dds2.Rdata")
        #}
        #### 第二步，提取差异分析结果
        #load(file = "data/Step03-DESeq2_dds2.Rdata")
        # 提取差异分析结果，trt组对untrt组的差异分析结果
        res <- results(dds2,contrast=c("group_list","decile1","decile10"))
        resOrdered <- res[order(res$padj),]
        head(resOrdered)
        normalized_counts <- counts(dds2, normalized=TRUE)
        
        DEG <- as.data.frame(resOrdered)
        # 去除差异分析结果中包含NA值的行
        DESeq2_DEG = na.omit(DEG)
        
        # 火山图 vocano
        g <- DESeq2_DEG
        g <- transform(g, padj = -1*log10(g$padj))
        g<-g%>%mutate(DiffExpressed=case_when(log2FoldChange<= -1 & padj>-log10(0.05)~"down",log2FoldChange >=1 & padj>-log10(0.05)~"up",TRUE~ "no"))
        t2<-ann[match(rownames(g),ann$Name),]
        g$label<-t2$Description
        g$label[g$DiffExpressed == "no"] <- NA
        fwrite(g,paste0(file,"_diffexp_ukb_decile1_9_10.tab"),sep="\t")
        # plot adding up all layers we have seen so far
        cairo_pdf(paste0(file,"_ukb_decile1_9_10.pdf"),heigh=7.08,width=5.90)
        ggp<-ggplot(data=g, aes(x=log2FoldChange, y=padj, col=DiffExpressed, label=label)) +
                geom_point() + 
                theme_minimal() +
                ggtitle(paste0(file," expression difference \nin decile 1-9 and decile 10"))+xlab("log2FoldChange")+ylab ("-log10(Qvalue)")+ggeasy::easy_center_title()+
                geom_text_repel() +
                scale_color_manual(values=c("blue", "black", "red")) +
                geom_vline(xintercept=c(-1,1), col="red") +
                geom_hline(yintercept=-log10(0.05), col="red")
        print(ggp)
        dev.off()
        return(g)
}

## Venn diagram

lista<-fread("gene_list_groupA")
listb<-fread("gene_list_groupB")
listc<-fread("gene_list_groupC")
venn.diagram(list(A=lista$V1,B=listb$V1,C=listc$V1),fill=c("red","yellow","blue"),cex=1.5,filename="venn.png")


#LIMMA

#samples<-samples%>%mutate(group=case_when(group==1~"0-90",group==2~">=90"))
#samples$group<-factor(samples$group,levels=c("0-90",">=90"))
#library('Homo.sapiens') #加载物种注释，这里以人为例
#keytypes(Homo.sapiens) #查看注释包的内容，里面包含了各种的ID以及基因位置信息
# 比如要匹配Ensembl的原始ID到gene ID，并且显示染色体信息以及基因位置
#gene_ids <- head(keys(Homo.sapiens, keytype='ENSEMBL'), 2)
#select(Homo.sapiens, keytype='ENSEMBL', keys=gene_ids, 
#       columns=c('ALIAS', 'TXCHROM', 'TXSTART', 'TXEND'))
# 就会得到类似下面的数据
## 'select()' returned 1:many mapping between keys and columns
# 用热图
# add a colorbar along the heatmap with sample condition

#{
#  num_conditions <- nlevels(samples$quantile)
#  pal <- colorRampPalette(brewer.pal(num_conditions, "Set1"))(num_conditions)
#  pal<-c("#E41A1C","#377EB8")
#  cond_colors <- pal[as.integer(samples$quantile)]
#md.pattern(raw_counts)
#  heatmap.2(cor(raw_counts), RowSideColors=cond_colors,
#            trace='none', main='Sample correlations (raw)')
#}

#用聚类
#{
#  sampleTree = hclust(dist(cor(raw_counts)), method = "average")
#  pdf(file = "pre-sampleClustering.pdf", width = 12, height = 9)
#  par(cex = 0.6)
#  par(mar = c(0,4,2,0))
#  plot(sampleTree, main = "Sample clustering to detect outliers", 
#       sub="", xlab="", cex.lab = 1.5,
#       cex.axis = 1.5, cex.main = 2)
#  dev.off()
#}

# 如果发现跑偏的样本，除掉它
#if(F){
#  abline(h = 1.0, col = "red") #画一条辅助线,h的值自定义
# 比如这里设置把高于20的切除
#  clust = cutreeStatic(sampleTree, cutHeight = 1, minSize = 10)
#  table(clust) # 0代表切除的，1代表保留的
#  keepSamples = (clust==1)
#  datExpr = datExpr0[keepSamples, ]
#}

# 基因在每个样本中平均表达量为1就要被过滤
#low_count_mask <- rowSums(raw_counts) < ncol(raw_counts)
#filt_raw_counts <- raw_counts[!low_count_mask,]
#sprintf("Removing %d low-count genes (%d remaining).",           sum(low_count_mask), sum(!low_count_mask))

# 注意log使用要将真数+1 （真数就是这里的filt_raw_counts）
#log_counts <- log2(filt_raw_counts + 1)

#然后画个图看看
#x = melt(as.matrix(log_counts))
#colnames(x) = c('gene_id', 'sample', 'value')
#ggplot(x, aes(x=value, color=sample)) + geom_density()+theme(legend.position = "none")

#再画个热图
#heatmap.2(cor(log_counts), RowSideColors=cond_colors,
#          trace='none', main='Sample correlations (log2-transformed)')
# 首先，去除方差为0的基因，因为这些基因没有表现出任何差别，还有可能对下面构建模型产生误导
#log_counts <- log_counts[apply(log_counts, 1, var) > 0,]

# 构建design矩阵为差异分析作准备
#mod <- model.matrix(~0+samples$quantile)
# 如果需要考虑其他分组信息，只需要加进来就好
# 例如：mod <- model.matrix(~0+samples$tissue+samples$cellline)
#colnames(mod) <- levels(samples$quantile)

#fit <- lmFit(log_counts, design=mod)
# 生成一系列的两两组合差异比较矩阵
#condition_pairs <- t(combn(levels(samples$quantile), 2))  
#comparisons <- list()                                                                                    
#for (i in 1:nrow(condition_pairs)) {                                                           
#  comparisons[[i]] <- as.character(condition_pairs[i,]) 
#}   

# 设置一个储存差异基因的向量
#sig_genes <- c()
# 为每对差异比较矩阵都生成差异基因，并存储到sig_genes空向量中
#for (conds in comparisons[[1]]) {
# 这里conds如果中间有空格，就需要用make.names变成标准名
#  contrast_formula <- paste(conds, collapse=' - ')
# 然后就是标准limma流程
#  contrast_mat <- makeContrasts(contrasts=contrast_formula, levels=mod)
#  contrast_fit <- contrasts.fit(fit, contrast_mat)
#  eb <- eBayes(contrast_fit)
#p值可以自定义，这里设置的比较严格
#  sig_genes <- union(sig_genes, 
#                     rownames(topTable(eb, number=Inf, p.value=0.005)))
#}
# 最后把差异不显著的基因去除，留下DEGs
#log_counts <- log_counts[rownames(log_counts) %in% sig_genes,]



#### WGCNA ####
#options(stringsAsFactors= FALSE)
# switch to TPM/FPKM log(x+1)
WGCNAfun<-function(exp,cutoff){
        # trait file 
        exp<-as.data.frame(exp)
        ann<-exp[,1:2]
        colnames(exp)[1]<-"Name"
        rownames(exp)<-exp$Name
        exp<-exp[,-c(1,2)]
        exp<-t(exp)
        exp<-as.data.frame(exp)
        exp$sample<-rownames(exp)
        exp<-exp %>% separate(sample,sep="-",into=c("p1","p2","p3","p4","p5"))
        cols=c("p1","p2")
        exp$samplef<-do.call(paste,c(exp[cols],sep="-"))
        t2<-gtex0[match(exp$samplef,gtex0$IID),]
        t2<-t2[is.na(t2$IID)==FALSE,]
        testcol<-match(t2$IID,exp$samplef)
        #exp2<-subset(exp,select=testcol)
        exp2<-exp[testcol,]
        testrow<-match(exp$samplef,gtex0$IID)
        samples<-gtex0[testrow,]
        samples<-samples[is.na(samples$IID)==FALSE,]
        #samples$quantile<-as.factor(samples$quantile)
        samples<-as.data.frame(samples)
        samples$sample<-rownames(exp2)
        raw_counts<-exp2[,-c((ncol(exp2)-5):ncol(exp2))]
        raw_counts<-t(raw_counts)
        # 基因在每个样本中平均表达量为1就要被过滤
        low_count_mask <- rowSums(raw_counts) < 10*ncol(raw_counts)
        raw_counts <- raw_counts[!low_count_mask,]
        sprintf("Removing %d low-count genes (%d remaining).",           sum(low_count_mask), sum(!low_count_mask))
        #normalized_counts <- counts(dds, normalized=TRUE)
        #提示输入数据不能为整数，我还不知道WGCNA应该用什么作为输入数据
        normalized_counts<-cpm(raw_counts, normalized.lib.sizes=TRUE)
        # select genes
        #expro<-data.frame(column_to_rownames(exp2, var = "Name"))
        #colnames(expro)<-colnames(exp2[,-1])
        #expro<-expro[,-1]
        expro<-normalized_counts
        #expro<-as.data.frame(t(exp2[,-c((ncol(exp2)-5):ncol(exp2))]))
        #row gene col sample
        ##筛选方差前25%的基因##
        m.vars=apply(expro,1,var)
        #expro.upper=expro[which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25),na.rm=TRUE)[4]),]
        expro.upper=expro[which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.05),na.rm=TRUE)[20]),]
        datExpr=as.data.frame(t(expro.upper));
        df2 <- data.frame(sapply(datExpr, function(x) as.numeric(as.character(x))))
        rownames(df2)<-rownames(datExpr)
        datExpr<-df2
        # gene col; sample row
        
        ##样本聚类检查离群值##
        gsg =goodSamplesGenes(datExpr)
        #gsg = goodSamplesGenes(datExpr, verbose = 3);
        if(!gsg$allOK)
        {
                #Optionally, print the gene and sample names that were removed:
                if(sum(!gsg$goodGenes)>0)
                        printFlush(paste("Removinggenes:", paste(names(datExpr)[!gsg$goodGenes], collapse = ",")));
                if(sum(!gsg$goodSamples)>0)
                        printFlush(paste("Removingsamples:", paste(rownames(datExpr)[!gsg$goodSamples], collapse = ",")));
                #Remove the offending genes and samples from the data:
                datExpr= datExpr[gsg$goodSamples, gsg$goodGenes]
        }
        sampleTree= hclust(dist(datExpr), method = "average");
        # Plotthe sample tree: Open a graphic output window of size 12 by 9 inches
        # Theuser should change the dimensions if the window is too large or too small.
        #sizeGrWindow(12,9)
        #pdf(file= "Plots/sampleClustering.pdf", width = 12, height = 9);
        #par(cex= 0.6);
        par(mar= c(0,4,2,0))
        #par(mfrow = c(2, 3)) # Create a 2 x 2 plotting matrix
        plot(sampleTree,main = "Sample clustering to detect outliers", sub="",xlab="", cex.lab = 1.5,
             cex.axis= 1.5, cex.main = 2)
        # Plota line to show the cut
        abline(h= cutoff, col = "red");
        #Determine cluster under the line
        clust= cutreeStatic(sampleTree, cutHeight = cutoff, minSize = 10)
        table(clust)
        #clust 1 contains the samples we want to keep.
        keepSamples= (clust==1)
        datExpr= datExpr[keepSamples, ]
        
        traitData<- samples
        #traitData<-traitData%>%mutate(quantile=case_when(quantile=="0-90"~1,quantile==">=90"~2))
        #traitData<-traitData%>%mutate(quantile=case_when(quantile=="decile1"~1,quantile=="decile10"~2))
        #remove columns that hold information we do not need.
        allTraits= traitData[, -c(3,4,5) ];
        # Forma data frame analogous to expression data that will hold the clinical traits.
        femaleSamples= rownames(datExpr);
        traitRows= match(femaleSamples, allTraits$sample);
        datTraits= allTraits[traitRows, ];
        rownames(datTraits)<-datTraits$sample
        #datTraits<-datTraits[,2,drop=FALSE]
        datTraits<-datTraits[,-c(1,2,18)]
        df2 <- data.frame(sapply(datTraits, function(x) as.numeric(as.character(x))))
        datTraits<-df2
        #Re-cluster samples
        sampleTree2= hclust(dist(datExpr), method = "average")
        #Convert traits to a color representation: white means low, red means high, greymeans missing entry
        traitColors= numbers2colors(as.numeric(unlist(datTraits)), signed = FALSE);
        # Plotthe sample dendrogram and the colors underneath.
        #plotDendroAndColors(sampleTree2,traitColors,
        #                    groupLabels= names(datTraits),
        #                    main ="Sample dendrogram and trait heatmap")
        #Choose a set of soft-thresholding powers
        powers= c(c(1:10), seq(from = 12, to=26, by=2))
        # Callthe network topology analysis function
        sft =pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
        # Plotthe results:
        sizeGrWindow(9,5)
        par(mfrow= c(1,2));
        cex1 =0.9;
        # Scale-freetopology fit index as a function of the soft-thresholding power
        plot(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
             xlab="SoftThreshold (power)",ylab="Scale Free Topology Model Fit,signedR^2",type="n",
             main =paste("Scale independence"));
        text(sft$fitIndices[,1],-sign(sft$fitIndices[,3])*sft$fitIndices[,2],
             labels=powers,cex=cex1,col="red");
        # thisline corresponds to using an R^2 cut-off of h
        abline(h=0.90,col="red")
        # Meanconnectivity as a function of the soft-thresholding power
        plot(sft$fitIndices[,1],sft$fitIndices[,5],
             xlab="SoftThreshold (power)",ylab="Mean Connectivity", type="n",
             main =paste("Mean connectivity"))
        text(sft$fitIndices[,1],sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
        sft$powerEstimate
        net =blockwiseModules(datExpr, power = sft$powerEstimate,maxBlockSize = 20000,
                              TOMType= "unsigned", minModuleSize = 30,
                              reassignThreshold= 0, mergeCutHeight = 0.25,
                              numericLabels= TRUE, pamRespectsDendro = FALSE,
                              saveTOMs= TRUE,
                              saveTOMFileBase= "femaleMouseTOM",
                              verbose= 3)
        #merge modules
        #merge = mergeCloseModules(datExpr, unmergedLabels, cutHeight = 0.25, verbose = 3)
        # opena graphics window
        sizeGrWindow(12,9)
        #Convert labels to colors for plotting
        mergedColors= labels2colors(net$colors)
        # Plotthe dendrogram and the module colors underneath
        plotDendroAndColors(net$dendrograms[[1]],mergedColors[net$blockGenes[[1]]],
                            "Modulecolors",
                            dendroLabels= FALSE, hang = 0.03,
                            addGuide= TRUE, guideHang = 0.05)
        #discard the unassigned genes and focus on the rest
        #restGenes= (mergedColors != "grey")
        #diss1=1-TOMsimilarityFromExpr(datExpr[,restGenes], power = sft$powerEstimate)
        #colnames(diss1) =rownames(diss1) =SubGeneNames[restGenes]
        #hier1=flashClust(as.dist(diss1), method="average" )
        #plotDendroAndColors(hier1, dynamicColors[restGenes], "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")
        
        
        moduleLabels= net$colors
        moduleColors= labels2colors(net$colors)
        MEs =net$MEs;
        geneTree= net$dendrograms[[1]]
        ##结果保存###
        #save(MEs, moduleLabels, moduleColors, geneTree,
        #     file = "AS-green-FPKM-02-networkConstruction-auto.RData")
        #Define numbers of genes and samples
        nGenes= ncol(datExpr);
        nSamples= nrow(datExpr);
        #Recalculate MEs with color labels
        MEs0 =moduleEigengenes(datExpr, moduleColors)$eigengenes
        MEs =orderMEs(MEs0)
        moduleTraitCor= cor(MEs, datTraits, use = "p");
        moduleTraitPvalue= corPvalueStudent(moduleTraitCor, nSamples);
        sizeGrWindow(10,6)
        # Willdisplay correlations and their p-values
        textMatrix= paste(signif(moduleTraitCor, 2), "\n(",
                          signif(moduleTraitPvalue,1), ")", sep = "");
        dim(textMatrix)= dim(moduleTraitCor)
        par(mar= c(6, 8.5, 3, 3));
        #Display the correlation values within a heatmap plot
        labeledHeatmap(Matrix= moduleTraitCor,
                       xLabels= names(datTraits),
                       yLabels= names(MEs),
                       ySymbols= names(MEs),
                       colorLabels= FALSE,
                       colors= blueWhiteRed(50),
                       textMatrix= textMatrix,
                       setStdMargins= FALSE,
                       cex.text= 0.5,
                       zlim =c(-1,1),
                       main =paste("Module-trait relationships"))
        #modules above a conventional cutoff = 0.15 were considered significant
        
        # seperate blocks - should be similar results
        #bwnet = blockwiseModules(datExpr, maxBlockSize = 2000,
        #                         power = 8, TOMType = "unsigned", minModuleSize = 30,
        #                         reassignThreshold = 0, mergeCutHeight = 0.25,
        #                         numericLbels = TRUE,
        #                         saveTOMs = TRUE,
        #                         saveTOMFileBase = "femaleMouseTOM-blockwise",
        #                         verbose = 3)
        #bwLabels= matchLabels(bwnet$colors, moduleLabels);
        #Convert labels to colors for plotting
        #bwModuleColors= labels2colors(bwLabels)
        # opena graphics window
        #sizeGrWindow(6,6)
        # Plotthe dendrogram and the module colors underneath for block 1
        #plotDendroAndColors(bwnet$dendrograms[[1]],bwModuleColors[bwnet$blockGenes[[1]]],
        #                    "Modulecolors", main = "Gene dendrogram and module colors in block 1",
        #                    dendroLabels= FALSE, hang = 0.03,
        #                    addGuide= TRUE, guideHang = 0.05)
        #sizeGrWindow(12,9)
        #plotDendroAndColors(geneTree,
        #                    cbind(moduleColors,bwModuleColors),
        #                    c("Singleblock", "Multi blocks"),
        #                    main ="Single block gene dendrogram and module colors",
        #                    dendroLabels= FALSE, hang = 0.03,
        #                    addGuide= TRUE, guideHang = 0.05)
        ####### END ###########
        
        #Define numbers of genes and samples
        #nGenes= ncol(datExpr);
        #nSamples= nrow(datExpr);
        #Recalculate MEs with color labels
        #MEs0 =moduleEigengenes(datExpr, bwModuleColors)$eigengenes
        #MEs =orderMEs(MEs0)
        #bwmoduleTraitCor= cor(MEs, datTraits, use = "p");
        #bwmoduleTraitPvalue= corPvalueStudent(bwmoduleTraitCor, nSamples);
        #sizeGrWindow(10,6)
        # Willdisplay correlations and their p-values
        #textMatrix= paste(signif(bwmoduleTraitCor, 2), "\n(",
        #                  signif(bwmoduleTraitPvalue,1), ")", sep = "");
        #dim(textMatrix)= dim(bwmoduleTraitCor)
        #par(mar= c(6, 8.5, 3, 3));
        #Display the correlation values within a heatmap plot
        #labeledHeatmap(Matrix= bwmoduleTraitCor,
        #               xLabels= names(datTraits),
        #               yLabels= names(MEs),
        #               ySymbols= names(MEs),
        #               colorLabels= FALSE,
        #               colors= blueWhiteRed(50),
        #               textMatrix= textMatrix,
        #               setStdMargins= FALSE,
        #               cex.text= 0.5,
        #               zlim =c(-1,1),
        #               main =paste("Module-trait relationships"))
        
        
        
        #Define variable weight containing the group column of datTrait
        group= as.data.frame(datTraits$PRS5);
        names(group)= "group"
        #names (colors) of the modules
        modNames= substring(names(MEs), 3)
        geneModuleMembership= as.data.frame(cor(datExpr, MEs, use = "p"));
        MMPvalue= as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
        
        names(geneModuleMembership)= paste("MM", modNames, sep="");
        names(MMPvalue)= paste("p.MM", modNames, sep="");
        geneTraitSignificance= as.data.frame(cor(datExpr, group, use = "p"));
        GSPvalue= as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
        names(geneTraitSignificance)= paste("GS.", names(group), sep="");
        names(GSPvalue)= paste("p.GS.", names(group), sep="");
        
        module= "red"
        column= match(module, modNames);
        moduleGenes= moduleColors==module;
        sizeGrWindow(7,7);
        par(mfrow= c(1,1));
        verboseScatterplot(abs(geneModuleMembership[moduleGenes,column]),abs(geneTraitSignificance[moduleGenes, 1]),
                           xlab =paste("Module Membership in", module, "module"),
                           ylab ="Gene significance for polygenic risk score",
                           main =paste("Module membership vs. gene significance\n"),col = module)
        
        # module significance
        GeneSignificance=abs(geneTraitSignificance[[1]])
        # Next module significance is defined as average gene significance.
        ModuleSignificance=tapply(GeneSignificance, moduleColors, mean, na.rm=T)
        sizeGrWindow(8,7)
        par(mfrow = c(1,1))
        plotModuleSignificance(GeneSignificance,moduleColors)
        
        ###导出网络到Cytoscape#### Recalculate topological overlap if needed
        TOM = TOMsimilarityFromExpr(datExpr, power =sft$powerEstimate);
        # Read in the annotation file# annot = read.csv(file = "GeneAnnotation.csv");
        # Select modules需要修改，选择需要导出的模块颜色
        modules = c("red");
        # Select module probes选择模块探测
        probes = names(datExpr)
        inModule = is.finite(match(moduleColors, modules));
        modProbes = probes[inModule];
        #modGenes = annot$gene_symbol[match(modProbes, annot$substanceBXH)];
        # Select the corresponding Topological Overlap
        modTOM = TOM[inModule, inModule];
        dimnames(modTOM) = list(modProbes, modProbes)
        # Export the network into edge and node list files Cytoscape can read
        cyt = exportNetworkToCytoscape(modTOM,
                                       edgeFile = paste(file,"_edges_",modules, ".txt", sep=""),
                                       nodeFile = paste(file,"_nodes_",modules, ".txt", sep=""),
                                       weighted = TRUE,
                                       threshold = 0.02,
                                       nodeNames = modProbes,                               
                                       #altNodeNames = modGenes,
                                       nodeAttr = moduleColors[inModule]);
        
        ## 可视化基因网络## 
        # Calculate topological overlap anew: this could be done more efficiently by saving the TOM
        # calculated during module detection, but let us do it again here.
        dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = sft$powerEstimate);
        # Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
        plotTOM = dissTOM^7;
        # Set diagonal to NA for a nicer plot
        diag(plotTOM) = NA;
        # Call the plot function#sizeGrWindow(9,9)
        TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
        ##TOMplot(plotTOM, geneTree, bwmoduleColors, main = "Network heatmap plot, all genes")
        #随便选取1000个基因来可视化
        nSelect = 1000
        # For reproducibility, we set the random seed
        set.seed(10);
        select = sample(nGenes, size = nSelect);
        selectTOM = dissTOM[select, select];
        # There's no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
        selectTree = hclust(as.dist(selectTOM), method = "average")
        selectColors = moduleColors[select];
        # Open a graphical window#sizeGrWindow(9,9)
        # Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing# the color palette; setting the diagonal to NA also improves the clarity of the plot
        plotDiss = selectTOM^7;
        diag(plotDiss) = NA;
        TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")
        #此处画的是根据基因间表达量进行聚类所得到的各模块间的相关性图
        MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
        MET = orderMEs(MEs)
        #sizeGrWindow(7, 6) 
        plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2), plotDendrograms = FALSE, xLabelsAngle = 90)
}

file_list <- list.files(path="~/Desktop/Expression_20210620/data/GTExv8/expression/gene_counts/UKBB_decile1-9_10_gene_counts/",pattern="*.tab")
file_list <- list.files(path="~/Desktop/Expression_20210620/data/GTExv8/expression/gene_counts/GTExv8_decile1_10_gene/",pattern="*.tab")
twas<-fread("~/Desktop/Expression_20210620/data/TWAS_gene")
gwasgene<-fread("~/Desktop/Expression_20210620/data/GWAS_gene")
expgene<-fread("~/Desktop/Expression_20210620/data/case_control_exp")
#TWAS, GWAS overlap
siggene<-data.frame("file"=0,"siggene"=0,"GWAS"=0,"TWAS"=0,"expdiff"=0)
for (i in 1:length(file_list)){
        #exp<-fread(file_list[i])
        file<-str_split(file_list[i],"_")[[1]][1]
        #g<-diffplt(exp,file)
        g<-fread(paste0(file,"_diffexp_decile1_decile10.tab"))
        #g<-fread(paste0(file,"_diffexp_ukb_decile1_9_10.tab"))
        testtwas<-match(twas$Gene,g$label)
        testgwas<-match(gwasgene$Gene,g$label)
        testexp<-match(expgene$Gene,g$label)
        siggene[i,]<-cbind(file,nrow(g[g$DiffExpressed!="no",]),sum(is.na(testgwas)==FALSE),sum(is.na(testtwas)==FALSE),sum(is.na(testexp)==FALSE))
}
