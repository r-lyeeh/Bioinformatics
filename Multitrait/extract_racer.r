#!/usr/bin/env Rscript
args=commandArgs(trailingOnly=TRUE)
library(RACER)
library(ggplot2)
library(data.table)
setwd("/home/pang/Desktop/CAD_IS_PAD_20200121/data/PRS/regionalplot")
# trace("ldRACER",edit=TRUE)
# change to ldlink.nci.nih.gov

regionplot <- function(region,specific){
  region0<-fread(region)
  q<-region0[region0$r0==specific,]
  chr<-q[1,1][[1]]
  start<-q[1,2][[1]]
  end<-q[1,3][[1]]
  moloc<-fread(paste("moloc.ppa_snp",region,sep="_"))
  p<-moloc[moloc$region==specific,]
  leadsnp0<-p[7,4][[1]]
  # gene plot
  utils::data(hg19)
  colnames(hg19) = c("GENE_ID", "CHR", "TRX_START", "TRX_END", "LENGTH", "GENE_NAME", "TYPE")
  gene_sub = hg19[hg19$CHR == chr,]
  gene_sub = subset(gene_sub, gene_sub$TRX_START > (start-50000))
  gene_sub = subset(gene_sub, gene_sub$TRX_END < (end+50000))
  gene_sub = gene_sub[gene_sub$TYPE == "protein_coding",]
  gene_sub = gene_sub[,c(3,4,6)]
  gene_sub = reshape2::melt(gene_sub,id.vars = "GENE_NAME")
  gene_sub$y_value = as.numeric(as.factor(gene_sub$GENE_NAME))
  plot_lab = subset(gene_sub, gene_sub$variable == "TRX_END")
  c = ggplot2::ggplot(gene_sub, ggplot2::aes_string(x = "value", y = "y_value")) +
    ggplot2::geom_line(ggplot2::aes_string(group = "GENE_NAME"), size = 2) + ggplot2::theme_bw() +
    ggplot2::geom_text(data = plot_lab, ggplot2::aes_string(x = "value", y = "y_value", label = "GENE_NAME"),
                       hjust = -0.1,vjust = 0.3, size = 2.5) +
    ggplot2::theme(axis.title.y = ggplot2::element_text(color = "white", size = 28),
                   axis.text.y = ggplot2::element_blank(),
                   axis.ticks.y = ggplot2::element_blank()) + ggplot2::xlab(paste0("Chromosome ", chr, " Position")) +
    ggplot2::coord_cartesian(xlim = c(start,end), ylim = c(0,(max(gene_sub$y_value)+1)))
 
  #GWAS plot
  cad12<-fread(paste("CAD",region,specific,sep="_"))
  is12<-fread(paste("IS",region,specific,sep="_"))
  pad12<-fread(paste("PAD",region,specific,sep="_"))
  cad_gwas_f = RACER::formatRACER(assoc_data = cad12, chr_col = 1, pos_col = 2, p_col = 8)
  is_gwas_f = RACER::formatRACER(assoc_data = is12, chr_col = 1, pos_col = 2, p_col = 8)
  pad_gwas_f = RACER::formatRACER(assoc_data = pad12, chr_col = 1, pos_col = 2, p_col = 8)
  cad_gwas_f_ld = RACER::ldRACER(assoc_data = cad_gwas_f, rs_col = 3, pops = "EUR", lead_snp = leadsnp0)
  is_gwas_f_ld = RACER::ldRACER(assoc_data = is_gwas_f, rs_col = 3, pops = "EUR", lead_snp = leadsnp0)
  pad_gwas_f_ld = RACER::ldRACER(assoc_data = pad_gwas_f, rs_col = 3, pops = "EUR", lead_snp = leadsnp0)
  cad_gwas_f_ld$Disease<-"CAD"
  is_gwas_f_ld$Disease<-"IS"
  pad_gwas_f_ld$Disease<-"PAD"
  gwas<-rbind(cad_gwas_f_ld,is_gwas_f_ld,pad_gwas_f_ld)
  b=ggplot(gwas,aes(x=POS,y=LOG10P,color=Disease,shape=LD_BIN))+geom_point(alpha=0.8)+xlab(paste("Chromosome ",chr,"Position",sep=" ")) + ylab("-log10(p-value)")+theme_bw()
  ggpubr::ggarrange(b, c, heights = c(3,1), nrow = 2, ncol = 1,
                    common.legend = TRUE, legend = "right")
}
pdf(paste(args[1],args[2],".pdf",sep=""),onefile = F)
regionplot(args[1],args[2])
dev.off()

