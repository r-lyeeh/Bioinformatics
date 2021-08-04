msigdf.human<-msigdbr(species = "Homo sapiens",category="C2")
c2<-msigdf.human %>% dplyr::select(gs_id, entrez_gene) %>% as.data.frame
c2 <- msigdbr(species = "Homo sapiens",gs_cat== "C2")  %>% select(gs_id, entrez_gene) %>% as.data.frame
data(geneList, package="DOSE")
id = bitr(names(geneList), "ENTREZID", "SYMBOL", "org.Hs.eg.db")
geneList = geneList[names(geneList) %in% id[,1]]
names(geneList) = id[match(names(geneList), id[,1]), 2]   
de <- names(geneList)[1:100]
x <- enricher(de, TERM2GENE = c2)
y <- GSEA(geneList, TERM2GENE = c2)

dotplot(KEGG,split=".sign")+facet_grid(~.sign)
library(enrichplot)
#特定通路作图
gseaplot2(x,1,color="red",pvalue_table = T) # 按第一个做二维码图，并显示p值

kegg<-fread("../kegg_sep/Path_backup/C 00010 Glycolysis.txt",header=FALSE)
kegg<-fread("kegg_pathway_summ.txt")
library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)

t<-read.gmt("../c2.cp.kegg.v7.4.entrez.gmt")
t<-read.gmt("kegg_pathway_summ.txt")
x <- c("AASDH","ABCB11","ADAM12","ADAMTS16","ADAMTS18")
test = bitr(x, #数据集
            fromType="SYMBOL", #输入为SYMBOL格式
            toType="ENTREZID",  # 转为ENTERZID格式
            OrgDb="org.Hs.eg.db") #人类 数据库


res<-enricher(test$ENTREZID,TERM2GENE = c2)
gseaplot2(res,1,color="red",pvalue_table = T) # 按第一个做二维码图，并显示p值
dotplot(res)
barplot(res,showCategory = 20)
