rm(list=ls())
graphics.off()

significant_clusters_after_randomised_status <- read.csv("extended_data_figure_4.csv")
hist(significant_clusters_after_randomised_status$number_of_clusers_with_p_less_than_0.05,breaks =40,border = "white",freq = F,las=1,ylab="Relative frequency (%)",xlab="number of clusters (p<0.05)",ylim=c(0,0.08),main = "")
abline(v=33,lty=2,lwd=2)
