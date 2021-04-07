#!/raid2/zeng/bin/Rscript
findit<-function(plink_ld_file,this_pwd) {
setwd(this_pwd)

#plink_ld_file<-"snp_NotAvail_chr14_200kbld09.ld"
#this_pwd<-"/raid2/zeng/Working/Heritability_20160120/res_snplist1_cadsnp59"

ld<-read.table(plink_ld_file,header=T,as.is=T)
this_chr<-ld$CHR_A[1]

bimAvail<-read.table(paste("/raid2/zeng/work_zeng_summary/work_zeng_as_main_author/20150624-20190221_h2Johan/Heritability_Johan_20150624/prep_20151023/pre_h2_9study_20151023/plinkfile_all9study_Imput/all9study_Imput_chr",this_chr,".bim",sep=""),header=F,as.is=T)
#bimAvail<-read.table(paste("/raid2/zeng/Working/Heritability/Heritability_Johan_20150624/h2_20151023/pre_h2_9study_20151023/plinkfile_all9study_Imput/all9study_Imput_chr",this_chr,".bim",sep=""),header=F,as.is=T)
snpAvail<-bimAvail$V2

snp_queried<-unique(ld$SNP_A)
proxy_out<-c()
for(i in 1:length(snp_queried)) {
	w_ld<-which(ld$SNP_A==snp_queried[i] & ld$SNP_A!=ld$SNP_B & ld$SNP_B%in%snpAvail)
	if(length(w_ld)>0) {
	this_ld<-ld[w_ld,]
	w_proxy<-which.max(this_ld$R2)
	proxy_out<-rbind(proxy_out, this_ld[w_proxy,])
	}
}
write.table(proxy_out,file=paste("snp_NotAvail_BestProxy_chr",this_chr,sep=""),sep="\t",quote=F,row.names=F,col.names=F)
#CHR_A	BP_A	SNP_A	CHR_B	BP_B	SNP_B	R2
}


args<-commandArgs(TRUE)
findit(args[1],args[2])


