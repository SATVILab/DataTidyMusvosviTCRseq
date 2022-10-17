rm(list=ls())
graphics.off()


overlap_gliph2_tcrdist_correct_status_labels <- read.csv("Extended_Data_Figure_8a.csv")

par(mfrow=c(2,1))
pie(c(sum(!is.na(overlap_gliph2_tcrdist_correct_status_labels$Number.of.GLIPH2.TCRs)&overlap_gliph2_tcrdist_correct_status_labels$Number.of.overlap>0),
sum(!is.na(overlap_gliph2_tcrdist_correct_status_labels$Number.of.GLIPH2.TCRs)&overlap_gliph2_tcrdist_correct_status_labels$Number.of.overlap==0),
sum(is.na(overlap_gliph2_tcrdist_correct_status_labels$Number.of.GLIPH2.TCRs)&overlap_gliph2_tcrdist_correct_status_labels$Number.of.overlap==0)))

overlap_gliph2_tcrdist_randomised_status_labels <- read.csv("Extended_Data_Figure_8b.csv")
hist(overlap_gliph2_tcrdist_randomised_status_labels$prop_of_overlap,breaks = 40,ylab = "Count",xlab="Percent of differentially abundant clusters identified by GLIPH2 and TCRdist3 with CDR3β sequences (%)",main = "")
abline(v=sum(!is.na(overlap_gliph2_tcrdist_correct_status_labels$Number.of.GLIPH2.TCRs)&overlap_gliph2_tcrdist_correct_status_labels$Number.of.overlap>0)/dim(overlap_gliph2_tcrdist_correct_status_labels)[1]*100,lty=2,lwd=2)
  
  