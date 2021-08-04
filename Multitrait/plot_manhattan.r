library(qqman)
library(lattice)
setwd("/home/pang/Desktop/CAD_IS_PAD_20200121/data/SNP")
cad<-read.table("cadmeta",as.is=TRUE,header=TRUE)
is<-read.table("ismeta",as.is=TRUE,header=TRUE)
pad<-read.table("padmeta",as.is=TRUE,header=TRUE)
png("cadmeta.png")
manhattan(cad, chr="Chr", bp="BP", snp="Marker", p="pValue" )
dev.off()
png("ismeta.png")
manhattan(is, chr="Chr", bp="BP", snp="Marker", p="pValue" )
dev.off()
png("padmeta.png")
manhattan(pad, chr="Chr", bp="BP", snp="Marker", p="pValue" )
dev.off()

library(ggrepel)
cad01<-fread("manhplot/PAD.gwas")
#cad01<-ism%>%select("chr","BP","SNP","P")
colnames(cad01)<-c("CHR","BP","SNP","P")
cadsig<-fread("manhplot/PAD.sig3")
colnames(cadsig)<-c("CHR","BP","SNP","group","P","Gene","r0")
cad02<-cadsig %>% select(1,2,3,5)
cad03<-rbind(cad01,cad02)
deduped.data <- unique( cad03[ , 1:4 ] )
cad03<-deduped.data

# Prepare the dataset
don <- cad03 %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(cad03, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot) %>%
  
  # Add highlight and annotation information
  mutate( is_annotate=ifelse(SNP %in% cadsig$SNP, "yes", "no")) %>%
  mutate( is_highlight=ifelse(SNP %in% cadsig$SNP, "yes", "no"))

# Prepare X axis
#axisdf <- don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
axisdf <- don %>% group_by(CHR) %>% summarise(center=( max(BPcum) + min(BPcum) ) / 2 )

t2<-don[match(cadsig$SNP,don$SNP),]
cadsig$tot<-t2$tot
cadsig$BPcum<-t2$BPcum

#t2<-cadsig[match(don$SNO,cadsig$SNP),]
#don$Gene<-t2$Gene

cadsig$color<-"red"
cadsig$color[cadsig$group=="CADPAD"]<-"green"
cadsig$color[cadsig$group=="ISPAD"]<-"purple"
cadsig$color[cadsig$group=="CADISPAD"]<-"orange"

ggplot(don, aes(x=BPcum, y=-log10(P))) +
  
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 22 )) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  #scale_y_continuous(expand = c(0, 0),limits=c(0,23) ) +     # remove space between plot area and x axis
  
  # Add highlighted points
  geom_point(data=cadsig, color=cadsig$color, size=2) +
  geom_text(data=cadsig,aes(label=Gene),hjust=0, vjust=0,size=3) +
  # Add label using ggrepel to avoid overlapping
  #geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=SNP), size=2) +
  geom_hline(yintercept=6,color="red")+
  
  # Custom the theme:
  theme_bw() +ggtitle("PAD")+xlab("CHR")+
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
#1109 x627
cairo_pdf(filename='Figure03.pdf', width=9.25, height=5.38,pointsize = 12,family="Univers")


## summary
cad01<-fread("manhplot/summ.gwas")
#cad01<-ism%>%select("chr","BP","SNP","P")
colnames(cad01)<-c("CHR","BP","SNP","P")
cadsig<-fread("manhplot/summ.sig")
colnames(cadsig)<-c("CHR","BP","SNP","P","Gene")
cad02<-cadsig %>% select(1,2,3,4)
cad03<-rbind(cad01,cad02)
deduped.data <- unique( cad03[ , 1:4 ] )
cad03<-deduped.data

# Prepare the dataset
don <- cad03 %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(cad03, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot) %>%
  
  # Add highlight and annotation information
  mutate( is_annotate=ifelse(SNP %in% cadsig$SNP, "yes", "no")) %>%
  mutate( is_highlight=ifelse(SNP %in% cadsig$SNP, "yes", "no"))

# Prepare X axis
#axisdf <- don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
axisdf <- don %>% group_by(CHR) %>% summarise(center=( max(BPcum) + min(BPcum) ) / 2 )

t2<-don[match(cadsig$SNP,don$SNP),]
cadsig$tot<-t2$tot
cadsig$BPcum<-t2$BPcum

don$P<-factor(don$P,level=c("a","b","c","ab","ac","bc","abc"))
ggplot(don, aes(x=BPcum, y=P,color=P)) +
  geom_point()+
  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=Gene), size=2)+
  # Custom the theme:
  theme_bw() +ggtitle("Significant regions")+xlab("CHR")+ylab("Region")+
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
#1109 x627
cairo_pdf(filename='~/Desktop/CAD_IS_PAD_20200121/data/Final/figures/manha.pdf', width=9.25, height=5.38,pointsize = 12,family="Univers")
