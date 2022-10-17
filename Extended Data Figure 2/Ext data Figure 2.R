
rm(list=ls())
graphics.off()

study_cluster_overlap <- read.csv("extended data figure 2 raw data.csv")
par(mar=c(6,6,6,6))
pie(table(study_cluster_overlap$var1))
