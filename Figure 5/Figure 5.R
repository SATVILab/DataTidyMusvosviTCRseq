rm(list=ls())
graphics.off()
library(rio)


cluster_freqs <- import("Frequencies_of_Mtb_gliph_clusters.csv")


signficant_cluster_order_in_figure5 <- c("DRB1*15 SVAL","DRB1*15 VALL","DRB3*01 SRDN%P","DRB5*01 ALFG","DQA1*01 SLQG%GYE","DRB1*15 SVAL%GNT","DRB1*15 GEAK","DPA1*01 SLG%PNTE","DRB1*13 R%GDYG","DRB1*15 SVALLG%T","DQB1*06 GEAK","DRB1*03 %PGWGMNTE","DQA1*01 GEAK","DRB3*01 %PGWGMNTE","DQB1*06 R%SGGEAKNI","DPA1*02 S%RQGAGYG","DQA1*01 GGKG%QP","DRB5*01 ST%NTE","DQA1*01 SS%RGTE","DRB4*01 SP%RSE","DRB1*04 S%LAAGQET","DRB3*02 SRDKG%NQP","DPB1*105 SLG%WET","DRB1*03 S%EDRGNTE","DRB1*03 R%TGPNE","DRB3*01 S%EDRGNTE","DPB1*105 RVLG%NE","DRB3*01 R%TGPNE","DRB1*03 RVLG%NE","DRB1*03 RR%GPNE")

par(mar=c(3,6,3,3))
par(mfrow=c(6,5))
for(i in signficant_cluster_order_in_figure5){
  temp <- subset(cluster_freqs,cluster_freqs$`HLA allele:GLIPH2 combination`==i)
  boxplot(subset(temp,temp$Group=="Controller")$`Frequency (%)`,subset(temp,temp$Group=="Progressor")$`Frequency (%)`,cex=0,range = 0,col="white",border=c("blue","red"),ylab="Frequency (%)",las=1)
  points(rep(1,times=length(subset(temp,temp$Group=="Controller")$`Frequency (%)`)),subset(temp,temp$Group=="Controller")$`Frequency (%)`,pch=19,cex=1.5,col="blue")
  
  points(rep(2,times=length(subset(temp,temp$Group=="Progressor")$`Frequency (%)`)),subset(temp,temp$Group=="Progressor")$`Frequency (%)`,pch=19,cex=1.5,col="red")
  mtext(side = 3,line = 0.5,i)
  mtext(side = 3,line = 1.5,paste("p = ",round(wilcox.test(subset(temp,temp$Group=="Controller")$`Frequency (%)`,subset(temp,temp$Group=="Progressor")$`Frequency (%)`)$p.value,digits = 4),sep = ""))
  mtext(side = 1,at = 1,line = 1.5,length(subset(temp,temp$Group=="Controller")$`Frequency (%)`))
  mtext(side = 1,at = 2,line = 1.5,length(subset(temp,temp$Group=="Progressor")$`Frequency (%)`))
}


cmv_clusters <- read.csv("Figure 5d controller_vs_progressor_cmv_clusters.csv")
ebv_clusters <- read.csv("Figure 5d controller_vs_progressor_ebv_clusters.csv")
flu_clusters <- read.csv("Figure 5d controller_vs_progressor_flu_clusters.csv")
mtb_clusters <- read.csv("Figure 5d controller_vs_progressor_mtb_clusters.csv")

sign_cmv_clusters <- dim(subset(cmv_clusters,cmv_clusters$Mann.Wh_pvalue<=0.05&cmv_clusters$fdr<0.2))[1]/dim(cmv_clusters)[1]*100

sign_ebv_clusters <- dim(subset(ebv_clusters,ebv_clusters$Mann.Wh_pvalue<=0.05&ebv_clusters$fdr<0.2))[1]/dim(ebv_clusters)[1]*100

sign_flu_clusters <- dim(subset(flu_clusters,flu_clusters$Mann.Wh_pvalue<=0.05&flu_clusters$fdr<0.2))[1]/dim(flu_clusters)[1]*100

sign_mtb_clusters <- dim(subset(mtb_clusters,mtb_clusters$Mann.Wh_pvalue<=0.05&mtb_clusters$fdr<0.2))[1]/dim(mtb_clusters)[1]*100


par(mfrow=c(2,2))
barplot(c(sign_cmv_clusters,sign_ebv_clusters,sign_flu_clusters,sign_mtb_clusters),ylim=c(0,20),las=1,names.arg = c("cmv","ebv","flu","mtb"),ylab="TCR specificity group:HLA combinations (%)",main="figure 5d")


randomized_disease_outcome <- read.csv("Figure 5e.csv")
hist(randomized_disease_outcome$number_of_sig_clusters,breaks=30,border = "white",xlim=c(0,40),freq = F,ylab="Proportion",las=1,main = "figure 5e")
abline(v=dim(subset(mtb_clusters,mtb_clusters$Mann.Wh_pvalue<=0.05&mtb_clusters$fdr<0.2))[1])



cmv_clusters_for_cmv_analyis <- read.csv("Figure 5f cmv+ vs cmv- cmv clusters.csv")
ebv_clusters_for_cmv_analyis <- read.csv("Figure 5f cmv+ vs cmv- ebv clusters.csv")
flu_clusters_for_cmv_analyis <- read.csv("Figure 5f cmv+ vs cmv- Flu clusters.csv")
mtb_clusters_for_cmv_analyis <- read.csv("Figure 5f cmv+ vs cmv- Mtb clusters.csv")

sign_cmv_clusters_for_cmv_analyis <- dim(subset(cmv_clusters_for_cmv_analyis,cmv_clusters_for_cmv_analyis$Mann.Wh_pvalue<=0.05&cmv_clusters_for_cmv_analyis$fdr<0.2))[1]/dim(cmv_clusters_for_cmv_analyis)[1]*100

sign_ebv_clusters_for_cmv_analyis <- dim(subset(ebv_clusters_for_cmv_analyis,ebv_clusters_for_cmv_analyis$Mann.Wh_pvalue<=0.05&ebv_clusters_for_cmv_analyis$fdr<0.2))[1]/dim(ebv_clusters_for_cmv_analyis)[1]*100

sign_flu_clusters_for_cmv_analyis <- dim(subset(flu_clusters_for_cmv_analyis,flu_clusters_for_cmv_analyis$Mann.Wh_pvalue<=0.05&flu_clusters_for_cmv_analyis$fdr<0.2))[1]/dim(flu_clusters_for_cmv_analyis)[1]*100

sign_mtb_clusters_for_cmv_analyis <- dim(subset(mtb_clusters_for_cmv_analyis,mtb_clusters_for_cmv_analyis$Mann.Wh_pvalue<=0.05&mtb_clusters_for_cmv_analyis$fdr<0.2))[1]/dim(mtb_clusters_for_cmv_analyis)[1]*100



barplot(c(sign_cmv_clusters_for_cmv_analyis,sign_ebv_clusters_for_cmv_analyis,sign_flu_clusters_for_cmv_analyis,sign_mtb_clusters_for_cmv_analyis),ylim=c(0,30),las=1,names.arg = c("cmv","ebv","flu","mtb"),ylab="TCR specificity group:HLA combinations (%)",main="figure 5f")


