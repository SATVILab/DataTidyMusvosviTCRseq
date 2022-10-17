rm(list=ls())
graphics.off()
library(effsize)

freq_of_clusters_in_controllers_and_progressors <- import("Frequencies_of_Mtb_gliph_clusters_in_controllers_and_progressors_with_spline_timepoint.csv")

signficant_cluster_order_in_figure5 <- c("DRB1*15 SVAL","DRB1*15 VALL","DRB3*01 SRDN%P","DRB5*01 ALFG","DQA1*01 SLQG%GYE","DRB1*15 SVAL%GNT","DRB1*15 GEAK","DPA1*01 SLG%PNTE","DRB1*13 R%GDYG","DRB1*15 SVALLG%T","DQB1*06 GEAK","DRB1*03 %PGWGMNTE","DQA1*01 GEAK","DRB3*01 %PGWGMNTE","DQB1*06 R%SGGEAKNI","DPA1*02 S%RQGAGYG","DQA1*01 GGKG%QP","DRB5*01 ST%NTE","DQA1*01 SS%RGTE","DRB4*01 SP%RSE","DRB1*04 S%LAAGQET","DRB3*02 SRDKG%NQP","DPB1*105 SLG%WET","DRB1*03 S%EDRGNTE","DRB1*03 R%TGPNE","DRB3*01 S%EDRGNTE","DPB1*105 RVLG%NE","DRB3*01 R%TGPNE","DRB1*03 RVLG%NE","DRB1*03 RR%GPNE")

par(mfrow=c(6,5))
significant_gliph2 <- subset(freq_of_clusters_in_controllers_and_progressors,freq_of_clusters_in_controllers_and_progressors$`HLA allele:GLIPH2 combination`%in%signficant_cluster_order_in_figure5)
for(i in signficant_cluster_order_in_figure5){
  


sp.frame <- subset(significant_gliph2,significant_gliph2$Group=="Controller"&significant_gliph2$`HLA allele:GLIPH2 combination`==i,select=c("Slipe.Timepoint","Frequency (%)","donor_id"))
sp.frame <- as.data.frame(sp.frame)
sp.frame$days_to_tb <- ifelse(sp.frame$Slipe.Timepoint>650,650,sp.frame$Slipe.Timepoint)
sp.resampler <- function() {
  n <- nrow(sp.frame)
  resample.rows <- sample(1:n,size=n,replace=TRUE)
  return(sp.frame[resample.rows,])
}

grid.300 <- seq(from=min(sp.frame$days_to_tb),to=max(sp.frame$days_to_tb),length.out=300)

sp.spline.estimator <- function(data,m=300) {
  # Fit spline to data, with cross-validation to pick lambda
  fit <- smooth.spline(x=data[,1],y=data[,2],cv=TRUE,df = 4)
  # Set up a grid of m evenly-spaced points on which to evaluate the spline
  eval.grid <- seq(from=min(sp.frame$days_to_tb),to=max(sp.frame$days_to_tb),length.out=m)
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
              x=seq(from=min(sp.frame$days_to_tb),to=max(sp.frame$days_to_tb),length.out=m)))
}

#set.seed(100)
set.seed(100)
sp.cis <- sp.spline.cis(B=2000,alpha=0.05)

plot(sp.frame$days_to_tb,sp.frame$`Frequency (%)`,xlab="Days To TB",ylab="Freq",xlim=c(650,0),ylim=c(0,max(subset(significant_gliph2,significant_gliph2$`HLA allele:GLIPH2 combination`==i)$`Frequency (%)`)),cex=0)



for(temp.donors in unique(sp.frame$donor_id)){
  temp.points <- subset(sp.frame,sp.frame$donor_id==temp.donors)
  points(temp.points$days_to_tb,temp.points$`Frequency (%)`,cex=3,col="#519cd6",pch=19)
  text(temp.points$days_to_tb,temp.points$`Frequency (%)`,labels=unique(subset(significant_gliph2,significant_gliph2$donor_id==temp.donors)$Number.of.Spline.dots),col="white",cex=0.7,font = 2)
  
}

lines(x=sp.cis$x,y=sp.cis$main.curve,col="blue",lwd=2)
mtext(side = 3,i)

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


sp.frame <- subset(significant_gliph2,significant_gliph2$Group=="Progressor"&significant_gliph2$`HLA allele:GLIPH2 combination`==i,select=c("Slipe.Timepoint","Frequency (%)","donor_id"))
sp.frame$days_to_tb <- ifelse(sp.frame$Slipe.Timepoint>650,650,sp.frame$Slipe.Timepoint)
sp.resampler <- function() {
  n <- nrow(sp.frame)
  resample.rows <- sample(1:n,size=n,replace=TRUE)
  return(sp.frame[resample.rows,])
}

grid.300 <- seq(from=min(sp.frame$days_to_tb),to=max(sp.frame$days_to_tb),length.out=300)

sp.spline.estimator <- function(data,m=300) {
  # Fit spline to data, with cross-validation to pick lambda
  fit <- smooth.spline(x=data[,1],y=data[,2],cv=TRUE,df = 4)
  # Set up a grid of m evenly-spaced points on which to evaluate the spline
  eval.grid <- seq(from=min(sp.frame$days_to_tb),to=max(sp.frame$days_to_tb),length.out=m)
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
              x=seq(from=min(sp.frame$days_to_tb),to=max(sp.frame$days_to_tb),length.out=m)))
}

#set.seed(100)
set.seed(100)

j=1
for(temp.donors in unique(sp.frame$donor_id)){
  temp.points <- subset(sp.frame,sp.frame$donor_id==temp.donors)
  points(temp.points$days_to_tb,temp.points$`Frequency (%)`,cex=3,col="red",pch=19)
  text(temp.points$days_to_tb,temp.points$`Frequency (%)`,labels=unique(subset(significant_gliph2,significant_gliph2$donor_id==temp.donors)$Number.of.Spline.dots),col="white",cex=0.7,font = 2)
  j=j+1
}


#points(sp.frame$days_to_tb,sp.frame$`Frequency (%)`,cex=1.5,col="red",pch=19)
sp.cis <- sp.spline.cis(B=200,alpha=0.05)
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
}




effect_size_df <- NULL
for(i in signficant_cluster_order_in_figure5){
  temp <- subset(significant_gliph2,significant_gliph2$`HLA allele:GLIPH2 combination`==i)
  temp_meta <- as.data.frame(i)
  colnames(temp_meta) <- "allele_cluster"
  temp_meta$all_cohorts_p_value <- round(wilcox.test(subset(temp,temp$Group=="Controller")$`Frequency (%)`,subset(temp,temp$Group=="Progressor")$`Frequency (%)`)$p.value,digits = 4)
  temp_meta$all_cohorts_cliff_delta <- round(cliff.delta(subset(temp,temp$Group=="Controller")$`Frequency (%)`,subset(temp,temp$Group=="Progressor")$`Frequency (%)`)$estimate,digits = 4)
  
  
  temp_cohort <- subset(temp,temp$cohort=="ACS")
  temp_meta$acs_p_value <- round(wilcox.test(subset(temp_cohort,temp_cohort$Group=="Controller")$`Frequency (%)`,subset(temp_cohort,temp_cohort$Group=="Progressor")$`Frequency (%)`)$p.value,digits = 4)
  temp_meta$acs_cliff_delta <- round(cliff.delta(subset(temp_cohort,temp_cohort$Group=="Controller")$`Frequency (%)`,subset(temp_cohort,temp_cohort$Group=="Progressor")$`Frequency (%)`)$estimate,digits = 4)
  
  
  
  temp_cohort <- subset(temp,temp$cohort=="GC6")
  temp_meta$gc6_p_value <- round(wilcox.test(subset(temp_cohort,temp_cohort$Group=="Controller")$`Frequency (%)`,subset(temp_cohort,temp_cohort$Group=="Progressor")$`Frequency (%)`)$p.value,digits = 4)
  temp_meta$gc6_cliff_delta <- round(cliff.delta(subset(temp_cohort,temp_cohort$Group=="Controller")$`Frequency (%)`,subset(temp_cohort,temp_cohort$Group=="Progressor")$`Frequency (%)`)$estimate,digits = 4)
  
  
  effect_size_df <- rbind.data.frame(effect_size_df,temp_meta)
}


par(mfrow=c(6,5))
par(mar=c(4,7,4,2))
for(i in signficant_cluster_order_in_figure5){
  temp <- subset(significant_gliph2,significant_gliph2$`HLA allele:GLIPH2 combination`==i)
  temp_acs<- subset(temp,temp$cohort=="ACS")
  temp_gc6<- subset(temp,temp$cohort=="GC6")
  
  
  boxplot(subset(temp,temp$Group=="Controller")$`Frequency (%)`,subset(temp,temp$Group=="Progressor")$`Frequency (%)`,NA,
          subset(temp_acs,temp_acs$Group=="Controller")$`Frequency (%)`,subset(temp_acs,temp_acs$Group=="Progressor")$`Frequency (%)`,NA,
          subset(temp_gc6,temp_gc6$Group=="Controller")$`Frequency (%)`,subset(temp_gc6,temp_gc6$Group=="Progressor")$`Frequency (%)`
          ,col=c("#7ecce6","#e57571",NA),names = c(length(subset(temp,temp$Group=="Controller")$`Frequency (%)`),length(subset(temp,temp$Group=="Progressor")$`Frequency (%)`),NA,length(subset(temp_acs,temp_acs$Group=="Controller")$`Frequency (%)`),length(subset(temp_acs,temp_acs$Group=="Progressor")$`Frequency (%)`),NA,length(subset(temp_gc6,temp_gc6$Group=="Controller")$`Frequency (%)`),length(subset(temp_gc6,temp_gc6$Group=="Progressor")$`Frequency (%)`)),cex.axis=1.5,las=1,cex=0)
  
  
  mtext(side = 3,line =1,subset(effect_size_df,effect_size_df$allele_cluster==i)$all_cohorts_p_value,at = 1.5,cex = 0.75)
  mtext(side = 3,line = 0,subset(effect_size_df,effect_size_df$allele_cluster==i)$all_cohorts_cliff_delta,at = 1.5,cex = 0.75)
  
  mtext(side = 3,line = 1,subset(effect_size_df,effect_size_df$allele_cluster==i)$acs_p_value,at = 4.5,cex = 0.75)
  mtext(side = 3,line = 0,subset(effect_size_df,effect_size_df$allele_cluster==i)$acs_cliff_delta,at = 4.5,cex = 0.75)
  
  mtext(side = 3,line = 1,subset(effect_size_df,effect_size_df$allele_cluster==i)$gc6_p_value,at = 7.5,cex = 0.75)
  mtext(side = 3,line = 0,subset(effect_size_df,effect_size_df$allele_cluster==i)$gc6_cliff_delta,at = 7.5,cex = 0.75)
  
  
  points(rep(1,times=length(subset(temp,temp$Group=="Controller")$`Frequency (%)`)),subset(temp,temp$Group=="Controller")$`Frequency (%)`,cex=1,pch=19,col="#364f92")
  points(rep(2,times=length(subset(temp,temp$Group=="Progressor")$`Frequency (%)`)),subset(temp,temp$Group=="Progressor")$`Frequency (%)`,cex=1,pch=19,col="#db3428")
  
  points(rep(4,times=length(subset(temp_acs,temp_acs$Group=="Controller")$`Frequency (%)`)),subset(temp_acs,temp_acs$Group=="Controller")$`Frequency (%)`,cex=1,pch=19,col="#364f92")
  points(rep(5,times=length(subset(temp_acs,temp_acs$Group=="Progressor")$`Frequency (%)`)),subset(temp_acs,temp_acs$Group=="Progressor")$`Frequency (%)`,cex=1,pch=19,col="#db3428")
  
  points(rep(7,times=length(subset(temp_gc6,temp_gc6$Group=="Controller")$`Frequency (%)`)),subset(temp_gc6,temp_gc6$Group=="Controller")$`Frequency (%)`,cex=1,pch=19,col="#364f92")
  points(rep(8,times=length(subset(temp_gc6,temp_gc6$Group=="Progressor")$`Frequency (%)`)),subset(temp_gc6,temp_gc6$Group=="Progressor")$`Frequency (%)`,cex=1,pch=19,col="#db3428")
  
  mtext(side = 2,line = 5,paste("Freq. of ",unlist(strsplit(i," ",fixed = T))[2]," (%)"),cex=0.8)
  mtext(side = 3,line = 2.5,i)
  
  mtext(side = 3,line = 1,at = -1,"p value =",cex=0.8)
  mtext(side = 3,line = 0,at = -1,"effect size =",cex=0.8)
  
}
