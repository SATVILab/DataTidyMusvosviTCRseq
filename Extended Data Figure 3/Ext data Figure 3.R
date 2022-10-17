rm(list=ls())
graphics.off()

hla_association_data <- read.csv("Extended Data Figure 3.csv")

acs_hla_controller_vs_progressor_data <- subset(hla_association_data,hla_association_data$cohort=="ACS")
gc6_hla_controller_vs_progressor_data <- subset(hla_association_data,hla_association_data$cohort=="GC6")

acs_hla_controller_vs_progressor_data <- acs_hla_controller_vs_progressor_data[order(acs_hla_controller_vs_progressor_data$fisher_test_pvalue,decreasing = F),]
gc6_hla_controller_vs_progressor_data <- gc6_hla_controller_vs_progressor_data[order(gc6_hla_controller_vs_progressor_data$fisher_test_pvalue,decreasing = F),]

acs_hla_allele <- acs_hla_controller_vs_progressor_data[,1]
acs_hla_controller_vs_progressor_data <- t(acs_hla_controller_vs_progressor_data[,2:6])
colnames(acs_hla_controller_vs_progressor_data) <- acs_hla_allele
        
gc6_hla_allele <- gc6_hla_controller_vs_progressor_data[,1]
gc6_hla_controller_vs_progressor_data <- t(gc6_hla_controller_vs_progressor_data[,2:6])
colnames(gc6_hla_controller_vs_progressor_data) <- gc6_hla_allele

par(mfrow=c(2,1))
par(mar=c(8,4,4,4))
barplot(acs_hla_controller_vs_progressor_data[1:4,],las=2,col = c("#3f59a8","#498cca","#ee3625","#f79421"),las=2,ylim=c(0,100),yaxt="n",main = "Extended data figure 3 ACS",names.arg = paste(colnames(acs_hla_controller_vs_progressor_data[1:4,]),round(acs_hla_controller_vs_progressor_data[5,],digits = 2),sep = " = "))
axis(side = 2,at = seq(from=0,to=84,by=12),las=2,cex.axis=2)

barplot(gc6_hla_controller_vs_progressor_data[1:4,],las=2,col = c("#3f59a8","#498cca","#ee3625","#f79421"),las=2,ylim=c(0,40),yaxt="n",main = "Extended data figure 3 GC6",names.arg = paste(colnames(gc6_hla_controller_vs_progressor_data[1:4,]),round(gc6_hla_controller_vs_progressor_data[5,],digits = 2),sep = " = "))
axis(side = 2,at = seq(from=0,to=40,by=4),las=2,cex.axis=2)




