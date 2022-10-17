rm(list=ls())
graphics.off()
library(tidyverse)
##Figure 2 plots


acs_flow_data <- read.csv("Figure 2a and 2b.csv")
plot(NA,NA,ylim=c(0.001,10),xlim=c(0,3),las=1,xaxt="n",xlab="",ylab="Freq. of activated T cells (%)",log="y")
mtext(side = 1,at = 0,"PBS",line = 0.5)
mtext(side = 1,at = 1,"M.tb Lysate",line = 0.5)
mtext(side = 1,at = 2,"PBS",line = 0.5)
mtext(side = 1,at = 3,"M.tb Lysate",line = 0.5)

for(i in 1:dim(acs_flow_data)[1]){
  temp <- acs_flow_data[i,]
  if(temp$donorGrp=="Non Progressor"){
    points(c(0,1),c(temp$Unstim.Freq,temp$Mtb.Lysate.Freq),col="blue",type="o",pch=19)
  }else{
    points(c(2,3),c(temp$Unstim.Freq,temp$Mtb.Lysate.Freq),col="red",type="o",pch=19)
  }
}

mtext(side = 3,line = 0.5,paste("p = ",wilcox.test(subset(acs_flow_data,acs_flow_data$donorGrp=="Non Progressor")$Unstim.Freq,subset(acs_flow_data,acs_flow_data$donorGrp=="Non Progressor")$Mtb.Lysate.Freq)$p.value),at = 0.5)

mtext(side = 3,line = 0.5,paste("p = ",wilcox.test(subset(acs_flow_data,acs_flow_data$donorGrp=="Progressor")$Unstim.Freq,subset(acs_flow_data,acs_flow_data$donorGrp=="Progressor")$Mtb.Lysate.Freq)$p.value),at = 2.5)


boxplot(subset(acs_flow_data,acs_flow_data$donorGrp=="Non Progressor")$BCK.sub,subset(acs_flow_data,acs_flow_data$donorGrp=="Progressor")$BCK.sub,pch=0,col = "white",cex=0,ylab="Freq. of activated T cells (%)",ylim=c(0,2),las=2,range = 0)
points(rep(1,times=length(subset(acs_flow_data,acs_flow_data$donorGrp=="Non Progressor")$BCK.sub)),subset(acs_flow_data,acs_flow_data$donorGrp=="Non Progressor")$BCK.sub,pch=19,cex=2,col="blue")
points(rep(2,times=length(subset(acs_flow_data,acs_flow_data$donorGrp=="Progressor")$BCK.sub)),subset(acs_flow_data,acs_flow_data$donorGrp=="Progressor")$BCK.sub,pch=19,cex=2,col="red")
mtext(side = 1,at = 1,line = 0.5, "Controller")
mtext(side = 1,at = 2,line = 0.5, "Progressor")
mtext(side = 3,line = 0.5,paste("p = ",wilcox.test(subset(acs_flow_data,acs_flow_data$donorGrp=="Non Progressor")$BCK.sub,subset(acs_flow_data,acs_flow_data$donorGrp=="Progressor")$BCK.sub)$p.value))


acs_cdr3b_detected <- read.csv("Figure 2c.csv")
boxplot(subset(acs_cdr3b_detected,acs_cdr3b_detected$donorGrp=="Non Progressor")$Total.CD3b_templates.detected,subset(acs_cdr3b_detected,acs_cdr3b_detected$donorGrp=="Progressor")$Total.CD3b_templates.detected,pch=0,col = "white",cex=0,ylab="Freq. of activated T cells (%)",ylim=c(0,500),las=2,range = 0)


points(rep(1,times=length(subset(acs_cdr3b_detected,acs_cdr3b_detected$donorGrp=="Non Progressor")$Total.CD3b_templates.detected)),subset(acs_cdr3b_detected,acs_cdr3b_detected$donorGrp=="Non Progressor")$Total.CD3b_templates.detected,pch=19,cex=2,col="blue")
points(rep(2,times=length(subset(acs_cdr3b_detected,acs_cdr3b_detected$donorGrp=="Progressor")$Total.CD3b_templates.detected)),subset(acs_cdr3b_detected,acs_cdr3b_detected$donorGrp=="Progressor")$Total.CD3b_templates.detected,pch=19,cex=2,col="red")
mtext(side = 1,at = 1,line = 0.5, "Controller")
mtext(side = 1,at = 2,line = 0.5, "Progressor")
mtext(side = 3,line = 0.5,paste("p = ",wilcox.test(subset(acs_cdr3b_detected,acs_cdr3b_detected$donorGrp=="Non Progressor")$Total.CD3b_templates.detected,subset(acs_cdr3b_detected,acs_cdr3b_detected$donorGrp=="Progressor")$Total.CD3b_templates.detected)$p.value))


acs_flow_kinetics <- read.csv("Figure 2d.csv")


key.for.spline.dots_MtblysateReactive <-   unique.data.frame(subset(acs_flow_kinetics,select=c("PID","Number.of.Spline.dots","donorGrp")))

B1 <- subset(acs_flow_kinetics,!is.na(acs_flow_kinetics$Freq))

sp.frame <- subset(B1,B1$donorGrp=="Non Progressor",select=c("Slipe.Timepoint","Freq","PID"))
sp.frame$Slipe.Timepoint <- ifelse(sp.frame$Slipe.Timepoint>650,650,sp.frame$Slipe.Timepoint)
sp.resampler <- function() {
  n <- nrow(sp.frame)
  resample.rows <- sample(1:n,size=n,replace=TRUE)
  return(sp.frame[resample.rows,])
}

grid.300 <- seq(from=min(sp.frame$Slipe.Timepoint),to=max(sp.frame$Slipe.Timepoint),length.out=300)

sp.spline.estimator <- function(data,m=300) {
  # Fit spline to data, with cross-validation to pick lambda
  fit <- smooth.spline(x=data[,1],y=data[,2],cv=TRUE,df = 4)
  # Set up a grid of m evenly-spaced points on which to evaluate the spline
  eval.grid <- seq(from=min(sp.frame$Slipe.Timepoint),to=max(sp.frame$Slipe.Timepoint),length.out=m)
  # Slightly inefficient to re-define the same grid every time we call this,
  # but not a big overhead
  # Do the prediction and return the predicted values
  return(predict(fit,x=eval.grid)$y)  # We only want the predicted values
}


sp.spline.cis <- function(B,alpha,m=300) {
  spline.main <- sp.spline.estimator(sp.frame,m=m)
  # Draw B boottrap samples, fit the spline to each
  spline.boots <- replicate(B,sp.spline.estimator(sp.resampler(),m=m))
  # Result has m rows and B columns
  cis.lower <- 2*spline.main - apply(spline.boots,1,quantile,probs=1-alpha/2)
  cis.upper <- 2*spline.main - apply(spline.boots,1,quantile,probs=alpha/2)
  return(list(main.curve=spline.main,lower.ci=cis.lower,upper.ci=cis.upper,
              x=seq(from=min(sp.frame$Slipe.Timepoint),to=max(sp.frame$Slipe.Timepoint),length.out=m)))
}

#set.seed(100)
set.seed(100)
sp.cis <- sp.spline.cis(B=2000,alpha=0.05)

plot(sp.frame$Slipe.Timepoint,sp.frame$Freq,xlab="Days To TB",ylab="Freq",xlim=c(650,0),ylim=c(0,2),cex=0)



for(temp.donors in unique(sp.frame$PID)){
  temp.points <- subset(sp.frame,sp.frame$PID==temp.donors)
  points(temp.points$Slipe.Timepoint,temp.points$Freq,cex=3,col="#519cd6",pch=19)
  text(temp.points$Slipe.Timepoint,temp.points$Freq,labels=unique(subset(key.for.spline.dots_MtblysateReactive,key.for.spline.dots_MtblysateReactive$PID==temp.donors)$Number.of.Spline.dots),col="white",cex=0.7,font = 2)
  
}

lines(x=sp.cis$x,y=sp.cis$main.curve,col="blue",lwd=2)
mtext(side = 3,"Freq of Mtb reactive cells (BCK_sub)")

t_col <- function(color, percent = 20, name = NULL) {
  
  #	  color = color name
  #	percent = % transparency
  #	   name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100-percent)*255/100,
               names = name)
  
  ## Save the color
  invisible(t.col)
  
}

mycol <- t_col("blue", perc = 80, name = "lt.pink")

polygon(c(sp.cis$x,900,rev(sp.cis$x)),c(sp.cis$lower.ci,900,rev(sp.cis$upper.ci)),col =mycol,lwd=0.0001)


sp.frame <- subset(B1,B1$donorGrp=="Progressor",select=c("Slipe.Timepoint","Freq","PID"))
sp.frame$Slipe.Timepoint <- ifelse(sp.frame$Slipe.Timepoint>650,650,sp.frame$Slipe.Timepoint)
sp.resampler <- function() {
  n <- nrow(sp.frame)
  resample.rows <- sample(1:n,size=n,replace=TRUE)
  return(sp.frame[resample.rows,])
}

grid.300 <- seq(from=min(sp.frame$Slipe.Timepoint),to=max(sp.frame$Slipe.Timepoint),length.out=300)

sp.spline.estimator <- function(data,m=300) {
  # Fit spline to data, with cross-validation to pick lambda
  fit <- smooth.spline(x=data[,1],y=data[,2],cv=TRUE,df = 4)
  # Set up a grid of m evenly-spaced points on which to evaluate the spline
  eval.grid <- seq(from=min(sp.frame$Slipe.Timepoint),to=max(sp.frame$Slipe.Timepoint),length.out=m)
  # Slightly inefficient to re-define the same grid every time we call this,
  # but not a big overhead
  # Do the prediction and return the predicted values
  return(predict(fit,x=eval.grid)$y)  # We only want the predicted values
}


sp.spline.cis <- function(B,alpha,m=300) {
  spline.main <- sp.spline.estimator(sp.frame,m=m)
  # Draw B boottrap samples, fit the spline to each
  spline.boots <- replicate(B,sp.spline.estimator(sp.resampler(),m=m))
  # Result has m rows and B columns
  cis.lower <- 2*spline.main - apply(spline.boots,1,quantile,probs=1-alpha/2)
  cis.upper <- 2*spline.main - apply(spline.boots,1,quantile,probs=alpha/2)
  return(list(main.curve=spline.main,lower.ci=cis.lower,upper.ci=cis.upper,
              x=seq(from=min(sp.frame$Slipe.Timepoint),to=max(sp.frame$Slipe.Timepoint),length.out=m)))
}

#set.seed(100)
set.seed(100)

j=1
for(temp.donors in unique(sp.frame$PID)){
  temp.points <- subset(sp.frame,sp.frame$PID==temp.donors)
  points(temp.points$Slipe.Timepoint,temp.points$Freq,cex=3,col="red",pch=19)
  text(temp.points$Slipe.Timepoint,temp.points$Freq,labels=unique(subset(key.for.spline.dots_MtblysateReactive,key.for.spline.dots_MtblysateReactive$PID==temp.donors)$Number.of.Spline.dots),col="white",cex=0.7,font = 2)
  j=j+1
}


#points(sp.frame$Slipe.Timepoint,sp.frame$Freq,cex=1.5,col="red",pch=19)
sp.cis <- sp.spline.cis(B=2001,alpha=0.05)
lines(x=sp.cis$x,y=sp.cis$main.curve,col="red",lwd=2)
#lines(x=sp.cis$x,y=sp.cis$lower.ci,col="red",lwd=2)
#lines(x=sp.cis$x,y=sp.cis$upper.ci,col="red",lwd=2)

t_col <- function(color, percent = 20, name = NULL) {
  
  #	  color = color name
  #	percent = % transparency
  #	   name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               max = 255,
               alpha = (100-percent)*255/100,
               names = name)
  
  ## Save the color
  invisible(t.col)
  
}

mycol <- t_col("red", perc = 80, name = "lt.pink")

polygon(c(sp.cis$x,900,rev(sp.cis$x)),c(sp.cis$lower.ci,900,rev(sp.cis$upper.ci)),col =mycol,lwd=0.0001)  


TCR_b_visitcount <- read.csv("Figure 2e.csv")


TCR_beta_count_P <- TCR_b_visitcount[which(TCR_b_visitcount$donorGrp=="Progressor"), ]
TCR_beta_count_NP <- TCR_b_visitcount[which(TCR_b_visitcount$donorGrp=="NonProgressor"), ]

tcrcount <- TCR_beta_count_NP
ggplot(tcrcount, aes(donorId, log2(n), color=donorId)) +
  geom_jitter(size=log2(tcrcount$n), position=position_jitter(width=.2)) +
  scale_y_continuous(name="log2(Counts)", limits=c(-0.2, 5)) +
  facet_grid(visit ~ .) +
  theme_bw(base_size=16) + 
  theme(axis.text.x=element_text(angle=90), legend.position = "none")+
  ggtitle("Non Progressor")


tcrcount <- TCR_beta_count_P
ggplot(tcrcount, aes(donorId, log2(n), color=donorId)) +
  geom_jitter(size=log2(tcrcount$n), position=position_jitter(width=.2)) +
  scale_y_continuous(name="log2(Counts)", limits=c(-0.2, 5)) +
  facet_grid(visit ~ .) +
  theme_bw(base_size=16) + 
  theme(axis.text.x=element_text(angle=90), legend.position = "none")+
  ggtitle("Progressor")