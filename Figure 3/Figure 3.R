rm(list=ls())
graphics.off()
library(tidyverse)
##Figure 3 plots

paired_blood_lung_tcr_data <- read.csv("Figure 3.csv")

par(mfrow=c(2,2))
for(i in c("Mtb","CMV","EBV","InfluenzaA")){
 temp <- subset(paired_blood_lung_tcr_data,paired_blood_lung_tcr_data$Pathogen==i)
  
 plot(NA,NA,ylim=c(0,max(temp$Blood,temp$Lung,na.rm = T)),xlim=c(0,1),las=1,xaxt="n",xlab="",ylab=paste("Freq. of ",i," CDR3b sequences (%)",sep = ""))
 mtext(side = 1,at = 0,"Blood",line = 0.5)
 mtext(side = 1,at = 1,"Lung",line = 0.5)
 
 for(i_1 in 1:dim(temp)[1]){
   points(c(0,1),c(temp[i_1,]$Blood,temp[i_1,]$Lung),col=c("#ee3224","#91278f"),type="o",pch=19)
 }
 mtext(side = 3,paste("p = ",wilcox.test(temp$Blood,temp$Lung,paired = T)$p.value))
}