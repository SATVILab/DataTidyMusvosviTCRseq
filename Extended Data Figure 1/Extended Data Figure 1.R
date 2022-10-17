rm(list=ls())
graphics.off()

par(mfrow=c(3,3))
par(mar=c(3,3,5,4))
cells_subsets <-  import("Extended Data Figure 1b.csv")
pie(c(cells_subsets$`# of CD4 T cells`,cells_subsets$`# of CD8 T cells`,cells_subsets$`# of MAIT cells`),labels=c("CD4","CD8","MAIT"),main="Ext data Figure 1b")

mfi_cd26 <- import("Extended Data Figure 1c.csv")
boxplot(mfi_cd26$cd26_cd4,mfi_cd26$cd26_cd8,mfi_cd26$cd26_mait,range=0,cex=0,col="white",ylim=c(0,24000),names=c("cd4","cd8","mait"),yaxt="n",main="Ext data Figure 1c")
axis(side = 2,at=seq(from=0,to=24000,by=4000),las=2)
points(rep(1,times=length(mfi_cd26$cd26_cd4)),mfi_cd26$cd26_cd4,pch=19,col="#a86faf")
points(rep(2,times=length(mfi_cd26$cd26_cd8)),mfi_cd26$cd26_cd8,pch=19,col="#784d23")
points(rep(3,times=length(mfi_cd26$cd26_mait)),mfi_cd26$cd26_mait,pch=19,col="#738838")

segments(x0=1,x1=2,y0=20000,y1=20000)
text(x = 1.5,y=20500,round(wilcox.test(mfi_cd26$cd26_cd4,mfi_cd26$cd26_cd8,exact = T)$p.value,digits = 4))

segments(x0=2,x1=3,y0=22000,y1=22000)
text(x = 2.5,y=22500,round(wilcox.test(mfi_cd26$cd26_cd8,mfi_cd26$cd26_mait,exact = T)$p.value,digits = 4))

segments(x0=1,x1=3,y0=24000,y1=24000)
text(x = 2,y=24500,round(wilcox.test(mfi_cd26$cd26_cd4,mfi_cd26$cd26_mait,exact = T)$p.value,digits = 4))


proportion_of_transcript_positive <- import("Extended Data Figure 1d.csv")

for(i in unique(proportion_of_transcript_positive$transcript)){
  temp <- subset(proportion_of_transcript_positive,proportion_of_transcript_positive$transcript==i)
  boxplot(temp$cd4,temp$cd8,temp$mait,range=0,cex=0,col="white",names=c("cd4","cd8","mait"),yaxt="n")
  points(rep(1,times=length(temp$cd4)),temp$cd4,pch=19,col="#a86faf")
  points(rep(2,times=length(temp$cd8)),temp$cd8,pch=19,col="#784d23")
  points(rep(3,times=length(temp$mait)),temp$mait,pch=19,col="#738838")
  mtext(side = 3,line = 0.5,paste("Ext data Figure 1d = ", i,sep = ""))
  mtext(side = 3,line = 1.5,paste("cd4 vs cd8 = ",round(wilcox.test(temp$cd4,temp$cd8,exact = T)$p.value,digits = 4),sep = ""))
  mtext(side = 3,line = 2.5,paste("cd8 vs mait = ",round(wilcox.test(temp$cd8,temp$mait,exact = T)$p.value,digits = 4),sep = ""))
  mtext(side = 3,line = 3.5,paste("cd4 vs mait = ",round(wilcox.test(temp$cd4,temp$mait,exact = T)$p.value,digits = 4),sep = ""))
}
