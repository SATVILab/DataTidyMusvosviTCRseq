rm(list=ls())
graphics.off()
library(tidyverse)
library(rio)
library(pheatmap)
##Figure 4 plots

cluster_matrix <- import("Figure 4a.csv")

subject_group <- cluster_matrix[, c("donorId", "Donor.Grouping")]
subject_group$Donor.Grouping <- as.factor(cluster_matrix$Donor.Grouping)
subject_group <- column_to_rownames(subject_group, "donorId")

cluster_matrix <- cluster_matrix[,3:dim(cluster_matrix)[2]]
rownames(cluster_matrix) <- rownames(subject_group)

pheatmap(cluster_matrix, cluster_rows = FALSE, cluster_cols = FALSE, scale = "column",
         annotation_row = subject_group,
         fontsize_row = 4, fontsize_col = 4, 
         cellwidth = 2.5, cellheight = 2.5,
         border_color = "grey70",
         color = colorRampPalette(c("white", "white", "navy"))(50),
         main = "All alleles")



number_of_np_donors_with_cluster <- matrix(NA,nrow = 1,ncol = 0)
number_of_p_donors_with_cluster <- matrix(NA,nrow = 1,ncol = 0)
for(i in colnames(cluster_matrix)){
  temp <- subset(cluster_matrix,rownames(cluster_matrix)%in%rownames(subset(subject_group,subject_group$Donor.Grouping=="NonProgressor")),select=i)
                 
  temp_meta <- as.data.frame(sum(temp[,1]>0))
  colnames(temp_meta)<-i
  number_of_np_donors_with_cluster <- cbind(number_of_np_donors_with_cluster,temp_meta)
  
  temp <- temp <- subset(cluster_matrix,rownames(cluster_matrix)%in%rownames(subset(subject_group,subject_group$Donor.Grouping=="Progressor")),select=i)
  temp_meta <- as.data.frame(sum(temp[,1]>0))
  colnames(temp_meta)<-i
  number_of_p_donors_with_cluster <- cbind(number_of_p_donors_with_cluster,temp_meta)
  
}

par(mfrow=c(2,1))
barplot(as.numeric(number_of_np_donors_with_cluster[1,]),col="#3a54a3",ylim=c(0,12),las=2,names.arg = colnames(number_of_np_donors_with_cluster),cex.names = 0.3)
barplot(as.numeric(number_of_p_donors_with_cluster[1,]),col="#eb2627",ylim=c(0,12),las=2,names.arg = colnames(number_of_p_donors_with_cluster),cex.names = 0.3)

cluster_ratio_filtered <- read.csv("Figure 4b.csv")
ggplot(cluster_ratio_filtered, aes(x = alleles, y = ratio, color = donorGrp)) +
  geom_boxplot(position = position_dodge(width = 0.8), outlier.shape = NA) +
  geom_point(alpha = 0.5, position = position_jitterdodge(0.2)) +
  coord_cartesian(ylim = c(0, 10))+
  scale_y_continuous(breaks=seq(0.0, 10, 2))+
  labs(y = "Number of Clusters per 100 CDR3s", x = "HLA allele") +
  theme_linedraw(base_size = 16, base_family = "Arial") +
  theme(axis.text.x= element_text(angle = 90, vjust = 0.5)) 


par(mfrow=c(1,1))
acs_hla_association <- import("Figure 4c.csv")
barplot(as.matrix(acs_hla_association[2:dim(acs_hla_association)[2]]),las=2,col = c("#3f59a8","#498cca","#ee3625","#f79421"),ylim=c(0,100),yaxt="n")
axis(side = 2,at = seq(from=0,to=84,by=12),las=2,cex.axis=2)




