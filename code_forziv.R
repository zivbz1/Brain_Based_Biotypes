# "Biotyping" analysis - Wards agglomerative hierarchical clust ---------------------------------------------------------
#data import and distribution checking
library(readxl)
roi_data <- read_excel("roi_data_ziv.xlsx")
roi_data <- na.omit(roi_data) #remove all cases with any missing data (n=40 total - 22 fz1, 18 fz2) - interpolation seems like a bad idea in this application
roi_data_raw <- roi_data #maintain a raw non-scaled dataset, since roi_data will be scaled with both freezes in next setp

#1- outlier correction and data characterization
library("PerformanceAnalytics")
chart.Correlation(roi_data_raw[1:69,4:12], histogram=TRUE, pch=19)
roi_datafz1<- roi_data_raw[1:69,]
roi_datafz1[,c(4:12)] <- scale(roi_datafz1[,c(4:12)]) #scale to z-distribution
roi_datafz1[,4:12][roi_datafz1[4:12]>=3] <- 3 #replace outliers by capping
roi_datafz1[,4:12][roi_datafz1[4:12]<=-3] <- -3 #replace outliers by capping

library("Hmisc")
mat1 <- rcorr(as.matrix(roi_datafz1[,4:12]))
mat1
mat1p <- mat1$P
mat2 <- cor(roi_datafz1[,4:12], use="complete.obs")
round(mat2, 2)

library(corrplot)
col <- colorRampPalette(c("#4477AA", "#77AADD", "#FFFFFF","#EE9988", "#BB4444"))
corrplot(mat2, method="color", col=col(200),  
         type="upper", order="hclust", 
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # Combine with significance
         p.mat = mat1p, sig.level = 0.05, insig = "blank", 
         # hide correlation coefficient on the principal diagonal
         diag=FALSE)


library(factoextra)
distancefz1 <- get_dist(roi_datafz1[,c(4:12)])
fviz_dist(distancefz1, show_labels = FALSE)

#2- cluster
library(cluster)
library(factoextra)
library(NbClust)
library(dendextend)
hc1 <- agnes(roi_datafz1[4:12], method = "ward")
row.names(hc1$data)<-c(1:69)
hc1$ac
pltree(hc1, cex = 0.6, hang = -1, main = "Clustering, Cohort 1") 

nb1 <- NbClust(roi_datafz1[,4:12], distance = "euclidean", min.nc = 2,
               max.nc = 6, method = "ward.D2", index="alllong")
fviz_nbclust(nb1)

fviz_nbclust(roi_datafz1[,4:12], hcut, method = "wss", k.max=6) +
  geom_vline(xintercept = 4, linetype = 2, color="steelblue")

fviz_nbclust(roi_datafz1[,4:12], hcut, method = "silhouette", k.max=6)

gap_stat <- clusGap(roi_datafz2[,4:12], FUN = hcut, K.max = 6, nboot=500)
fviz_gap_stat(gap_stat)

sub_grp<- cutree(as.hclust(hc1), k = 4)

cols_branches <- c("#E69F00", "#56B4E9", "#009E73", "#999999")
plot(col_dend1<-color_branches(as.dendrogram(hc1), k=4, col=cols_branches))


#4- Visualize clusters - PCA method
library(ggplot2)
roi_datafz1$sub_grp <- sub_grp
roi_datafz1$sub_grp <- as.factor(roi_datafz1$sub_grp)

p1 <- ggplot(roi_datafz1, aes(x=Dim.1, y=Dim.2, color=sub_grp)) #show clusters against PC1 & 2 
p1 + geom_point() + scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73", "#999999")) +
  xlim(-4.5, 4.5) + ylim(-4, 4.5) + 
  labs(x="PC1 - Threat", y="PC2 - Reward", color="Cluster")  +
  theme(axis.title.x = element_text(color="black", size=14, face="bold"), 
        axis.title.y = element_text(color="black", size=14, face="bold"),
        axis.text.x = element_text(color="black", size=12),
        axis.text.y = element_text(color="black", size=12))

#4B - Visualize clusters - Pattern across ROIs
library(ggplot2)
library("reshape2")
library(plyr)

datafz1_long <- melt(roi_datafz1, id.vars=c("PID", "site", "freeze", "sub_grp"), variable.name="roi")
datafz1_long <- datafz1_long[1:621,]
fz1summary <- ddply(datafz1_long, c("sub_grp", "roi"), summarise, 
                    N=length(value), mean=mean(value), sd=sd(value), se=sd/sqrt(N))
fz1summary$sub_grp<-as.factor(fz1summary$sub_grp)
fz1summary$subgrp_lab<-recode(fz1summary$sub_grp, '1'="Reactive/Disinhibited",
                              '2'="Low Reward/High Threat", '3'="High Reward", '4'="Inhibited")

fz1summary$task<- sub("threat.*", "threat", fz1summary$roi)
fz1summary$task<- sub("reward.*", "reward", fz1summary$task)
fz1summary$task<- sub("inhib.*", "inhib", fz1summary$task)

ggplot(data=fz1summary, aes(x=roi, y=mean, group=sub_grp, color=sub_grp)) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.1) + 
  geom_line() + facet_wrap(fz1summary$subgrp_lab, ncol=2) +
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73", "#999999")) +
  ylim(-2, 2) +
  labs(x=" ", y="Contrast Est., Z", color="Cluster") +
  theme(axis.title.y = element_text(color="black", size=14, face="bold"),
        axis.text.y = element_text(color="black", size=12),
        axis.text.x = element_text(color="black", size=4))

#5 - Bootstrap
library(fpc) 
cboot.hclust1<-clusterboot(roi_datafz1[,4:12], B=1000, bootmethod="boot", clustermethod=hclustCBI, method="ward.D2", k=4, scaling=FALSE, recover=0.6)
cboot.hclust1$bootmean 
cboot.hclust1$bootrecover 

###code for checking to make sure cboot does the same as agnes###
library(dendextend)
hc1.D2<-hclust(dist(roi_datafz1[,4:12]), method="ward.D2") 
plot(hc1.D2,  hang=-1)
dend_list <- dendlist(as.dendrogram(hc1), as.dendrogram(hc1.D2))
cor.dendlist(dend_list, method = "cophenetic")
dend_list %>% untangle(method="step1side") %>% tanglegram()

#6 - Training and testing approach to assess replication
#We might need to think about how to do this; probably need to meet/brainstorm after you have a dataset ready
library(class)
library(gmodels)
library(caret)

roifz1_noclust3 <- subset(roi_datafz1,sub_grp!="3") #optional - this was to drop our non-replicated cluster from fz1
roi_datafz2$sub_grp<- recode(roi_datafz2$sub_grp, "3"="4") #assign same number codes to qualitatively similar clusters

clus_pred_fz1train <- knn(train=roifz1_noclust3[,4:12], test=roi_datafz2[,4:12], cl=roifz1_noclust3$sub_grp, k=3) #forward test of replication
fz2train_vs_fz2clust <- CrossTable(roi_datafz2$sub_grp, clus_pred_fz1train, digits=2, chisq=TRUE)
confusionMatrix(roi_datafz2$sub_grp, clus_pred_fz1train) # get accuracy of the predicted labels


#7- MANOVA across all the outcome trajectories
library(tidyverse)
library(ggpubr)
library(rstatix)
detach("package:dplyr", unload=TRUE)
library(plyr)
library(dplyr)

outcomes_long2z <- outcomes_long2
outcomes_long2z[,9:15]<- scale(outcomes_long2z[,9:15])
outcomes_long2z[,18:19]<- scale(outcomes_long2z[,18:19])

summary(man1 <- manova(cbind(disab, pcl, dep, dissoc) ~ (sub_grp * time) + (sub_grp*time2) + freeze + Error(time/PID), data=outcomes_long2z)) 
#produced p=.03 for sub_grp but no interaction with time. Use man2 which includes subs outx

summary(man2 <- manova(cbind(disab, pcl, dep, dissoc,alcsqrt,thcsqrt) ~  (sub_grp * time) +(sub_grp*time2) + freeze + Error(time/PID), data=outcomes_long2z))
#F(12,1190)=1.91, p=.03 for sub_grp, F(6,594)=8.86, p=2.85*10^-9 for time^2.  


#MANOVA graph
sxs_emm = outcomes_long2z %>% 
  group_by(PID) %>% 
  summarise_at(c("disab","pcl","dep","dissoc", "alcsqrt", "thcsqrt"), mean) %>% 
  ungroup()
mergetemp <- select(alldata, "PID","sub_grp")
sxs_emm <- merge(sxs_emm, mergetemp, by="PID") #produce dataset that's collapsed over timepts

sxsemm_long<- na.omit(melt(sxs_emm, id.vars=c("PID", "sub_grp"), variable.name="outx")) #outcomes as rows
sxsemm_summary <- ddply(sxsemm_long, c("sub_grp", "outx"), summarise, 
                        N=length(value), mean=mean(value), sd=sd(value), se=sd/sqrt(N)) #get descriptives for graph
sxsemm_summary$sub_grp<-as.factor(sxsemm_summary$sub_grp)

ggplot(data=sxsemm_summary, aes(x=outx, y=mean, group=sub_grp, color=sub_grp)) + 
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=0.1) + 
  geom_line() +
  facet_wrap(~sub_grp)+
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#009E73", "#999999")) +
  labs(x=" ", y="Standardized score", color="Cluster") +
  theme(axis.title.x = element_text(color="black", size=14, face="bold"), 
        axis.title.y = element_text(color="black", size=14, face="bold"),
        axis.text.x = element_text(color="black", size=12),
        axis.text.y = element_text(color="black", size=12)) 

library(ggradar)
sxs_emm_xclust = read_xlsx("sxs_emm_xclust.xlsx")
ggradar(sxs_emm_xclust, group.colours=c("#E69F00", "#56B4E9", "#999999"), 
        values.radar = c("Z=-1", "0", "1"))  ###use this plot format###


#8 - Follow up analysis: Differences in each single outcome as a function of cluster
library(lme4)
library(lmerTest)
library(sjstats)

sink("./log_outcometests.txt")

disab_m1 <- lmer(disab ~ sub_grp*time + sub_grp*time2 + freeze + (1|PID), outcomes_long2)
print("Disability")
anova(disab_m1)
lsmeansLT(disab_m1)
difflsmeans(disab_m1)

pcl_m1 <- lmer(pcl ~ sub_grp*time + sub_grp*time2 + freeze + (time|PID), outcomes_long2)
print("PCL")
anova(pcl_m1)
lsmeansLT(pcl_m1)
difflsmeans(pcl_m1)

dep_m1 <- lmer(dep ~ sub_grp*time + sub_grp*time2 + freeze  + (time|PID), outcomes_long2)
print("Dep")
anova(dep_m1)
lsmeansLT(dep_m1)
difflsmeans(dep_m1)

dissoc_m1 <- lmer(dissoc ~ sub_grp*time + sub_grp*time2 + freeze + (1|PID), outcomes_long2)
print("Dissoc")
anova(dissoc_m1)
lsmeansLT(dissoc_m1)
difflsmeans(dissoc_m1)

alc_m1 <- lmer(alcsqrt ~ sub_grp*time + sub_grp*time2 + freeze + (1|PID), outcomes_long2)
print("Alcohol")
anova(alc_m1)
lsmeansLT(alc_m1)
difflsmeans(alc_m1) #Don't use this model - too skewed - loglinear model conducted in SPSS

thc_m1 <- lmer(thcsqrt ~ sub_grp*time + sub_grp*time2 + freeze + (time|PID), outcomes_long2)
print("Marijuana")
anova(thc_m1)
lsmeansLT(thc_m1)
difflsmeans(thc_m1) #Don't use this model - too skewed - loglinear model conducted in SPSS

sink()


timelabs=c("Pre", "2wk", "8wk", "3mo", "6mo") #make better labels for plots

ggplot(data=outcomes_long2, aes(x=time, y=dissoc, color=sub_grp)) +  #Dissoc graph
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#999999")) +
  geom_smooth(level=0.682, span=1)+
  facet_wrap(outcomes_long2$sub_grp, ncol=3)+
  labs(x="Time", y="Dissociation symptoms", color="Cluster") + 
  scale_x_continuous(labels=timelabs) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold"), 
        axis.title.y = element_text(color="black", size=14, face="bold"),
        axis.text.x = element_text(color="black", size=12, angle=40),
        axis.text.y = element_text(color="black", size=12))

ggplot(data=outcomes_long2, aes(x=time, y=pcl, color=sub_grp)) +  #PTSD graph
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#999999")) +
  geom_smooth(level=0.682, span=1)+
  facet_wrap(outcomes_long2$sub_grp, ncol=3)+
  labs(x="Time", y="PTSD symptoms", color="Cluster") + 
  scale_x_continuous(labels=timelabs) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold"), 
        axis.title.y = element_text(color="black", size=14, face="bold"),
        axis.text.x = element_text(color="black", size=12, angle=40),
        axis.text.y = element_text(color="black", size=12))

ggplot(data=outcomes_long2, aes(x=time, y=dep, color=sub_grp)) +  #depression graph
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#999999")) +
  geom_smooth(level=0.682, span=1)+
  facet_wrap(outcomes_long2$sub_grp, ncol=3)+
  labs(x="Time", y="Depression symptoms", color="Cluster") +  scale_x_continuous(labels=timelabs) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold"), 
        axis.title.y = element_text(color="black", size=14, face="bold"),
        axis.text.x = element_text(color="black", size=12, angle=40),
        axis.text.y = element_text(color="black", size=12))

ggplot(data=outcomes_long2, aes(x=time, y=anx, color=sub_grp)) +  #anx graph
  scale_color_manual(values=c("#E69F00","#56B4E9", "#999999")) +
  geom_smooth(level=0.682, span=1)+
  facet_wrap(outcomes_long2$sub_grp, ncol=4)+
  labs(x="Time", y="Anxiety symptoms", color="Cluster") + scale_x_continuous(labels=timelabs) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold"), 
        axis.title.y = element_text(color="black", size=14, face="bold"),
        axis.text.x = element_text(color="black", size=12, angle=40),
        axis.text.y = element_text(color="black", size=12))

ggplot(data=outcomes_long2, aes(x=time, y=disab, color=sub_grp)) +  #disab graph
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#999999")) +
  geom_smooth(level=0.682, span=1)+
  facet_wrap(outcomes_long2$sub_grp, ncol=3)+
  labs(x="Time", y="Disability", color="Cluster") + scale_x_continuous(labels=timelabs) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold"), 
        axis.title.y = element_text(color="black", size=14, face="bold"),
        axis.text.x = element_text(color="black", size=12, angle=40),
        axis.text.y = element_text(color="black", size=12))

ggplot(data=outcomes_long2, aes(x=time, y=alcsqrt, color=sub_grp)) +  #alc graph
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#999999")) +
  geom_smooth(level=0.682, span=1)+
  facet_wrap(outcomes_long2$sub_grp, ncol=4)+
  labs(x="Time", y="Alcohol (sqrt drinks/mo)", color="Cluster") + scale_x_continuous(labels=timelabs) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold"), 
        axis.title.y = element_text(color="black", size=14, face="bold"),
        axis.text.x = element_text(color="black", size=12, angle=40),
        axis.text.y = element_text(color="black", size=12))

ggplot(data=outcomes_long2, aes(x=time, y=thcsqrt, color=sub_grp)) +  #thc graph
  scale_color_manual(values=c("#E69F00", "#56B4E9", "#999999")) +
  geom_smooth(level=0.682, span=1)+
  facet_wrap(outcomes_long2$sub_grp, ncol=4)+
  labs(x="Time", y="Marijuana (sqrt days/month)", color="Cluster") + scale_x_continuous(labels=timelabs) +
  theme(axis.title.x = element_text(color="black", size=14, face="bold"), 
        axis.title.y = element_text(color="black", size=14, face="bold"),
        axis.text.x = element_text(color="black", size=12, angle=40),
        axis.text.y = element_text(color="black", size=12))



