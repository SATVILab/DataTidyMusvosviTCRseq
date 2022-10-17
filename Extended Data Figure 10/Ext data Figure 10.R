rm(list=ls())
graphics.off()


durt_cells_freq <- read.csv("SuppFigure10_DURT_figure.csv")

par(mfrow=c(2,2))
for(i in c("Freq.MAIT.MATCH", "Freq.GammaDelta.Adaptive.Definition", "Freq.iNKT.Adaptive.Definition","Freq.GEM.Adaptive.Definition")){
  temp_controller <- subset(durt_cells_freq,durt_cells_freq$Group=="Controller",select = i)
  temp_progressor <- subset(durt_cells_freq,durt_cells_freq$Group=="Progressor",select = i)
  boxplot(temp_controller[,1],temp_progressor[,1],range = 0,cex=0,col="white")
  points(rep(1,times=length(temp_controller[,1])),temp_controller[,1],pch=19,cex=1,col="blue")
  points(rep(2,times=length(temp_progressor[,1])),temp_progressor[,1],pch=19,cex=1,col="red")
  mtext(side = 3,paste("p = ",round(wilcox.test(temp_controller[,1],temp_progressor[,1],exact = F)$p.value,digits = 4)))
  
}