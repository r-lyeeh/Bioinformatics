# Expression test 
# GTEx_v8 209 variants take into consideration
gtex<-fread("~/Desktop/Expression_20210620/data/GTExv8/genotype/241_common_nra_GTEX.profile")
starnet<-fread("~/Desktop/Expression_20210620/data/STARNET/241_common_starnet_nra.profile")
ukbb<-fread("~/Desktop/Expression_20210620/data/UKBB/241_common_nra_UKBB.profile")
gtexpheno<-fread("/raid3/Data_public/GTEx_phs000424.v8.p2/81713/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v8.p2.c1.GRU/PhenotypeFiles/phs000424.v8.pht002742.v8.p2.c1.GTEx_Subject_Phenotypes.GRU.txt")
gtexp<-gtexpheno%>%select("phv00169061.v8.p2.c1","phv00169163.v8.p2.c1","phv00169162.v8.p2.c1","phv00169164.v8.p2.c1")
colnames(gtexp)<-as.character(gtexp[1,])
gtexp<-gtexp[-1,]
t2<-gtexp[match(gtex$IID,gtexp$SUBJID),]
gtex$MHHRTDIS<-t2$MHHRTDIS

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
library(DESeq2)
library("purrr")
library(dplyr)
library(tidyr)
library(stringr)
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
file_list <- list.files(path="~/Desktop/Expression_20210620/data/GTExv8/expression")

#initiate a blank data frame, each iteration of the loop will append the data from the given file to this variable
res <- data.frame()

#had to specify columns to get rid of the total column
for (i in 1:length(file_list)){
  qs<-cbind(file_list[i],topbot(file_list[i]))
  res<-rbind(res,qs)
}


#GTEX 
gtex1<-gtex[gtex$SCORESUM<227.5690,]
gtex10<-gtex[gtex$SCORESUM>=227.5690,]
gtex1$group<-"0-90"
gtex10$group<-">=90"
gtex0<-rbind(gtex1,gtex10)
exp<-load(file_list[5])
t2<-gtex0[match(exp$samplef,gtex0$IID),]
exp$group<-t2$group
exp$SCORESUM<-t2$SCORESUM

# Co-expression network
library('gplots')
library('ggplot2')
library('knitr')
library('limma')
library('reshape2')
library('RColorBrewer')
library('WGCNA')

options(stringsAsFactors= FALSE)
femData= read.csv("~/Desktop/Expression/FemaleLiver-Data/LiverFemale3600.csv")
expro= as.data.frame(t(femData[, -c(1:8)]))##转置矩阵
names(expro)= femData$substanceBXH;
rownames(expro)= names(femData)[-c(1:8)];
expro<-t(expro)
##筛选方差前25%的基因##
m.vars=apply(expro,1,var)
expro.upper=expro[which(m.vars>quantile(m.vars, probs = seq(0, 1, 0.25),na.rm=TRUE)[4]),]
dim(expro.upper)
datExpr=as.data.frame(t(expro.upper));
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

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
sizeGrWindow(12,9)
#pdf(file= "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex= 0.6);
par(mar= c(0,4,2,0))
plot(sampleTree,main = "Sample clustering to detect outliers", sub="",xlab="", cex.lab = 1.5,
     cex.axis= 1.5, cex.main = 2)
# Plota line to show the cut
abline(h= 15, col = "red");
#Determine cluster under the line
clust= cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
table(clust)
#clust 1 contains the samples we want to keep.
keepSamples= (clust==1)
datExpr= datExpr[keepSamples, ]
nGenes= ncol(datExpr)
nSamples= nrow(datExpr)
traitData= read.csv("~/Desktop/Expression/FemaleLiver-Data/ClinicalTraits.csv");
dim(traitData)
names(traitData)
#remove columns that hold information we do not need.
allTraits= traitData[, -c(31, 16)];
allTraits= allTraits[, c(2, 11:36) ];
dim(allTraits)
names(allTraits)
# Forma data frame analogous to expression data that will hold the clinical traits.
#femaleSamples= rownames(datExpr);
#traitRows= match(femaleSamples, allTraits$Mice);
#datTraits= allTraits[traitRows, -1];
#rownames(datTraits)= allTraits[traitRows, 1]
#Re-cluster samples
#sampleTree2= hclust(dist(datExpr), method = "average")
#Convert traits to a color representation: white means low, red means high, greymeans missing entry
#traitColors= numbers2colors(datTraits, signed = FALSE);
# Plotthe sample dendrogram and the colors underneath.
#plotDendroAndColors(sampleTree2,traitColors,
#                    groupLabels= names(datTraits),
#                    main ="Sample dendrogram and trait heatmap")
#Choose a set of soft-thresholding powers
powers= c(c(1:10), seq(from = 12, to=20, by=2))
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
net =blockwiseModules(datExpr, power = 6,
                      TOMType= "unsigned", minModuleSize = 30,
                      reassignThreshold= 0, mergeCutHeight = 0.25,
                      numericLabels= TRUE, pamRespectsDendro = FALSE,
                      saveTOMs= TRUE,
                      saveTOMFileBase= "femaleMouseTOM",
                      verbose= 3)
# opena graphics window
sizeGrWindow(12,9)
#Convert labels to colors for plotting
mergedColors= labels2colors(net$colors)
# Plotthe dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]],mergedColors[net$blockGenes[[1]]],
                    "Modulecolors",
                    dendroLabels= FALSE, hang = 0.03,
                    addGuide= TRUE, guideHang = 0.05)
moduleLabels= net$colors
moduleColors= labels2colors(net$colors)
MEs =net$MEs;
geneTree= net$dendrograms[[1]]
##结果保存###
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "AS-green-FPKM-02-networkConstruction-auto.RData")

# seperate blocks
#bwnet= blockwiseModules(datExpr, maxBlockSize = 2000,
#                        power= 6, TOMType = "unsigned", minModuleSize = 30,
#                        reassignThreshold= 0, mergeCutHeight = 0.25,
#                        numericLabels= TRUE,
#                        saveTOMs= TRUE,
#                        saveTOMFileBase= "femaleMouseTOM-blockwise",
#                        verbose= 3)
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
# Plotthe dendrogram and the module colors underneath for block 2
#plotDendroAndColors(bwnet$dendrograms[[2]],bwModuleColors[bwnet$blockGenes[[2]]],
#                    "Modulecolors", main = "Gene dendrogram and module colors in block 2",
#                    dendroLabels= FALSE, hang = 0.03,
#                    addGuide= TRUE, guideHang = 0.05)
#sizeGrWindow(12,9)
#plotDendroAndColors(geneTree,
#                    cbind(moduleColors,bwModuleColors),
#                    c("Singleblock", "2 blocks"),
#                    main ="Single block gene dendrogram and module colors",
#                    dendroLabels= FALSE, hang = 0.03,
#                    addGuide= TRUE, guideHang = 0.05)
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
               colors= greenWhiteRed(50),
               textMatrix= textMatrix,
               setStdMargins= FALSE,
               cex.text= 0.5,
               zlim =c(-1,1),
               main =paste("Module-trait relationships"))



#Define variable weight containing the weight column of datTrait
weight= as.data.frame(datTraits$weight_g);
names(weight)= "weight"
#names (colors) of the modules
modNames= substring(names(MEs), 3)
geneModuleMembership= as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue= as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership)= paste("MM", modNames, sep="");
names(MMPvalue)= paste("p.MM", modNames, sep="");
geneTraitSignificance= as.data.frame(cor(datExpr, weight, use = "p"));
GSPvalue= as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance)= paste("GS.", names(weight), sep="");
names(GSPvalue)= paste("p.GS.", names(weight), sep="");

module= "brown"
column= match(module, modNames);
moduleGenes= moduleColors==module;
sizeGrWindow(7,7);
par(mfrow= c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes,column]),abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab =paste("Module Membership in", module, "module"),
                   ylab ="Gene significance for body weight",
                   main =paste("Module membership vs. gene significance\n"),
                   cex.main= 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

###导出网络到Cytoscape#### Recalculate topological overlap if needed
TOM = TOMsimilarityFromExpr(datExpr, power =6);
# Read in the annotation file# annot = read.csv(file = "GeneAnnotation.csv");
# Select modules需要修改，选择需要导出的模块颜色
modules = c("pink");
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
                               edgeFile = paste("AS-green-FPKM-One-step-CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("AS-green-FPKM-One-step-CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,                               
                               #altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule]);

## 可视化基因网络## 
# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 14);
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7;
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function#sizeGrWindow(9,9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
#随便选取1000个基因来可视化
nSelect = 1000
nSelect = 600
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
sizeGrWindow(7, 6) 
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2), plotDendrograms = FALSE, xLabelsAngle = 90)


