library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)
library(pathview)
library(enrichplot)

kegg<-read.gmt("../c2.cp.kegg.v7.4.entrez.gmt")
t<-read.gmt("kegg_pathway_summ.txt")
msigdf.human<-msigdbr(species = "Homo sapiens",category="C2")
c2<-msigdf.human %>% dplyr::select(gs_name, entrez_gene) %>% as.data.frame
msigdf.human<-msigdbr(species = "Homo sapiens",category="C5")
c5<-msigdf.human %>% dplyr::select(gs_name, entrez_gene) %>% as.data.frame
msigdf.human<-msigdbr(species = "Homo sapiens",category="H")
h<-msigdf.human %>% dplyr::select(gs_name, entrez_gene) %>% as.data.frame

x<-fread("../../0602/CAD_leadsnp_gene_500kb",header=FALSE)
x<-fread("../../0602/test_top100_v2",header=FALSE)
test = bitr(x$V1, #数据集
            fromType="SYMBOL", #输入为SYMBOL格式
            toType="ENTREZID",  # 转为ENTERZID格式
            OrgDb="org.Hs.eg.db") #人类 数据库


res<-enricher(test$ENTREZID,TERM2GENE = c2)
res<-enricher(test$ENTREZID,TERM2GENE = c5)
res<-enricher(test$ENTREZID,TERM2GENE = h)
res<-enricher(test$ENTREZID,TERM2GENE = kegg)
res<-enricher(test$ENTREZID,TERM2GENE = t)

rrt<-res@result
#y <- GSEA(geneList, TERM2GENE = c2)
#gseaplot2(y,1,color="red",pvalue_table = T) # 按第一个做二维码图，并显示p值
dotplot(res)
#dotplot(res,split=".sign")+facet_grid(~.sign)
barplot(res,showCategory = 40)

egoall <- enrichGO(test$SYMBOL, OrgDb=org.Hs.eg.db, ont='ALL',pAdjustMethod='BH', pvalueCutoff=0.05, qvalueCutoff=0.2, keyType='SYMBOL')
dotplot(egoall,title='Top5 GO terms of each sub-class',showCategory=5,split='ONTOLOGY')+ facet_grid(ONTOLOGY~.,scale="free")

barplot(res)
fwrite(res@result,"../../0602/PAD_leadsnp_gene_1000kb_c2.res")
fwrite(res@result,"../../0602/PAD_leadsnp_gene_1000kb_c5.res")
fwrite(res@result,"../../0602/PAD_leadsnp_gene_1000kb_h.res")
fwrite(res@result,"../../0602/PAD_leadsnp_gene_1000kb_kegg.res")
fwrite(res@result,"../../0602/PAD_leadsnp_gene_1000kb_self.res")

library(ggplot2)
library(RColorBrewer)
pathway<-fread("0616/Pathway_leadsnp_1000kb_unique")
colnames(pathway)[4]<-"P"
pathway<-pathway[order(pathway$Class)]
pathway$Class<-factor(pathway$Class,levels=unique(pathway$Class))
numColors <- length(levels(pathway$Class)) # How many colors you need
getColors <- function(pal) brewer.pal(brewer.pal.info[pal, "maxcolors"], pal)
myPalette <- getColors("Set1")[1:numColors]
# unique
myPalette<-c("#E41A1C","#000000","#000000","#000000","#377EB8","#000000","#000000","#000000","#4DAF4A","#984EA3","#FF7F00","#000000","#d8ca3c","#000000","#000000","#A65628","#F781BF")
names(myPalette) <- levels(pathway$Class) 
pathway$color<-myPalette[pathway$Class]
pathway$ID<-factor(pathway$ID,levels=unique(pathway$ID))
color<-unique(cbind(pathway$ID,pathway$color))

pdf("0616/Pathway_leadsnp_1000kb_unique.pdf")
#ggplot(pathway,aes(x=ID,y=Trait))+geom_point(aes(color=P,size=AssociatedGene))+theme_classic()+scale_color_gradient(low = "red",high = "yellow")+ coord_flip()+theme(axis.text.y = element_text(colour=myPalette[pathway$Class]))
ggplot(pathway,aes(x=ID,y=Trait))+geom_point(aes(color=P,size=AssociatedGene))+theme_classic()+scale_color_gradient(low = "red",high = "yellow")+ coord_flip()+theme(axis.text.y = element_text(colour=color[,2]))
dev.off()



pdf("../../0602/CAD_GO_test3.pdf")
rrt<-egoall@result
