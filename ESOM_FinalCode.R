library(tidyverse)
library(kohonen)
library(Umatrix)
library(randomForest)
library(cluster)
library(RColorBrewer)
library(vegan)
library(factoextra)
library(reshape2)
library(lemon)
library(gridExtra)
library(ggpubr)
library(caret)

#set working directories
setwd("/Users/keirajohnson/Box Sync/Keira_Johnson/GGM_Paper2Materials/Codes and Data")

#load data
#setwd(dname)
param <- c("Bicarbonate","Calcium", "Chloride","Potassium","Magnesium", "Sodium","Sulfate", "Silica")
Rtypes <- c("Sdep", "Conglomerate", "Carbonates", "Water", "Evaporites",
            "Mudstone", "Igmet", "Sandstone")
LCtypes <- c("Developed", "Cropland", "GrassShrub", "Forest", "Wetland", "IceSnow", "Barren")
attributes <- c("date", "siteno","River.x.x","joinkey","lat.x","long.x",
                "area","TDS","Flow", "pH", "temp",
                param, Rtypes)
id.vars <- c("date", "siteno", "lat.x", "long.x", "River.x.x","joinkey", "area",
             "TDS", "Flow", "pH", "temp")
dat <- readRDS("allriv_compwts04212022.RDS")

#center each variable to mean of 0 and unit variance
dat1 <- dat %>% select(Bicarbonate:Barren)
id.dat <- dat %>% select(all_of(id.vars))
raw.scaled <- as.data.frame(scale(dat1, center=TRUE,scale=TRUE),stringsAsFactors = FALSE)
raw.scaled <- as.data.frame(cbind(id.dat, raw.scaled))

dat_melt<-melt(dat, id.vars = id.vars)

chem.somdat <- raw.scaled[,param]

raw.scaled <- reshape2::melt(raw.scaled, id.vars=id.vars)

scale_chem <- raw.scaled[raw.scaled$variable %in% param,]

####Run ESOM####

#### Large SOMs (ESOM) with umatrix ####
#setwd(dname)
dat <- readRDS("allriv_compwts04212022.RDS")

#Normalizing data - needs to save mean and variance for later
z.means <- apply(dat[,param],2,mean, na.rm=TRUE)
z.vars <- apply(dat[,param],2, sd, na.rm=TRUE)
#write.csv (rbind(z.means,z.vars),"Data means and vars.csv",row.names=FALSE)

X1 <- scale(dat[,param], center=TRUE, scale = TRUE) #These are the normalized input data for the ESOM

# Make the ESOM
#res <-  esomTrain(X1, InitMethod = "norm_mean_2std", Epochs=100, Key = 1:nrow(X1))
#saveRDS(res, "ESOM_04242022.RDS") #also have "ESOM_03022022.RDS", "ESOM_04052022.RDS", "ESOM_04072022.RDS", "ESOM_04242022.RDS"

#ESOM Grace made - saved output for
res <- readRDS("ESOM_04242022.RDS")

#Make and plot the UMatrix
umatrix <-  umatrixForEsom(res$Weights, Lines=res$Lines, 
                           Columns=res$Columns, Toroid=res$Toroid)

pdf("Umatrix_plot_Dec2023.pdf", width = 8, height = 5.5, family = "Times")

plotMatrix(umatrix, Toroid=FALSE, BmSize = 2,
           DrawLegend = TRUE, Clean = TRUE)

dev.off()

#showMatrix3D(res$Umatrix)

#add unique ID for BMU
BestMatches <- as.data.frame(res$BestMatches)
#row of BM -1 times number of columns (82) plus column of BM (not sure why its like this..)
ids <- ((res$BestMatches[,2]-1)*length(unique(res$BestMatches[,3])))+res$BestMatches[,3]
#add unique ID to BMU
BestMatches<-cbind(BestMatches,ids)
names(BestMatches) <- c("obs", "row", "col", "id")

#Convert data to unnormalized results
#str(res)
Wts.raw <- res$Weights
Wts.unnorm <- (z.vars)*t(Wts.raw)+z.means
Wts.unnorm <- t(Wts.unnorm)

colnames(Wts.unnorm) <- colnames(dat[,param]) #move the column names over

#plot unnormalized compositional weights for each solute
pplots <- list()
x <- 1
for(i in 1:dim(Wts.unnorm)[2]){
  cname <- colnames(Wts.unnorm)[i]
  test <- matrix(Wts.unnorm[,i], nrow = 50, ncol=82, byrow=TRUE)
  pplots[[x]] <- plotMatrix(test, Toroid=FALSE, Clean = TRUE, Title=cname, DrawLegend=F, ColorStyle = "Pmatrix")
  x <- x+1
}

plot_legend<-plotMatrix(test, Toroid=FALSE, Clean = TRUE, Title=cname, DrawLegend=T, ColorStyle = "Pmatrix")+
  theme(legend.position = "bottom")+labs(fill="Compositional Solute Concentration")+
  guides(fill=guide_colorbar(title.position = "top"))

plot_legend

leg<-g_legend(plot_legend)

hlay <- rbind(c(1,1,1,2,2,2),
              c(1,1,1,2,2,2),
              c(4,4,4,5,5,5),
              c(4,4,4,5,5,5),
              c(6,6,6,7,7,7),
              c(6,6,6,7,7,7),
              c(8,8,8,9,9,9),
              c(8,8,8,9,9,9),
              c(3,3,3,10,10,10))

#plot
#library(gridExtra)
#setwd(fname)
pdf("SOM_parameterplots_March2024.pdf", family = "Times")
grid.arrange(pplots[[1]],pplots[[2]],leg, pplots[[3]],pplots[[4]],pplots[[5]],
             pplots[[6]],pplots[[7]],pplots[[8]], nrow=2, ncol=5, layout_matrix=hlay)
dev.off()

####Cluster ESOM####
weights_clust<-res$Weights

#kmeans using 6 clusters
set.seed(123)
kmeans_cluster<-kmeans(weights_clust, iter.max=50, nstart=50, centers = 6)

clusts<-kmeans_cluster$cluster

kmeans_clusts_mat<-matrix(clusts, nrow = 50, ncol=82, byrow=TRUE)

kmeans_clusts_mat_melt<-melt(kmeans_clusts_mat)

colnames(kmeans_clusts_mat_melt)<-c("row", "col", "clust")

BestMatches_cluster<-merge(BestMatches, kmeans_clusts_mat_melt, by=c("row", "col"))

table(clusts)

#make silhouette plots for clusters - order largest to smallest
cluster1<-c("2")
cluster2<-c("4")
cluster3<-c("5")
cluster4 <- c("3")
cluster5 <- c("6")
cluster6<- c("1")

#rename clusters so they are ordered largest to smallest
BestMatches_cluster$clust <- ifelse(BestMatches_cluster$clust %in% cluster1,1,
                                    ifelse(BestMatches_cluster$clust %in% cluster2, 2,
                                           ifelse(BestMatches_cluster$clust %in% cluster3,3,
                                                  ifelse(BestMatches_cluster$clust %in% cluster4,4,
                                                         ifelse(BestMatches_cluster$clust %in% cluster5, 5,
                                                                ifelse(BestMatches_cluster$clust %in% cluster6,6, BestMatches_cluster$clust))))))

#add clusters to data and plot composition of each cluster
dat$obs <- seq(1,nrow(dat),1)
cdat <- inner_join(dat,BestMatches_cluster[,c("obs", "clust")], by=c("obs"))

#write.csv(cdat, "Data_wClusters.csv")

colors <- brewer.pal(n=6,"Set1")

colors<-c("#4DAF4A", "#E41A1C", "#FF7F00", "#984EA3", "#FFFF33", "#377EB8")

pdf("Umatrix_Clusters_Jan2024.pdf", width = 8, height = 5.5)

#plot umatrix with clusters on top
plotMatrix(umatrix, as.matrix(BestMatches_cluster[,c("obs","row","col")]), Toroid=F, BmSize=4, DrawLegend=F, Clean=T, 
           Cls=BestMatches_cluster$clust, ClsColors = colors)

dev.off()

pdf("Sil_Width_Jan2024.pdf", width = 8, height = 4, family = "Times")

fviz_silhouette(sil)+
  labs(title="", fill="Cluster", col="Cluster", y="Silhouette width")+
  scale_color_manual(values=c("1" = colors[1], "2"= colors[2], "3" = colors[3],
                              "4" = colors[4], "5"=colors[5], "6"=colors[6]))+
  scale_fill_manual(values=c("1" = colors[1], "2"= colors[2], "3" = colors[3],
                             "4" = colors[4], "5"=colors[5], "6"=colors[6]))

dev.off()

id.vars <- c("date", "siteno", "lat.x", "long.x", "River.x.x","joinkey", "area",
             "TDS", "Flow", "pH", "temp","MAP","Alkalinity",
             "Aluminum","Bromide","Nitrate","Iron", "Phosphorus",Rtypes[-4], LCtypes, "clust", "obs")
cdat_m <- reshape2::melt(cdat, id.vars=id.vars)
labs <- c(paste("Cluster", seq(1,6,1)))

cdat_m$variable <- as.character(cdat_m$variable)
cdat_m$value<-as.numeric(cdat_m$value)
cdat_m$clust <- factor(cdat_m$clust)

pdf("Solute_Distribution_Jan2024.pdf", width = 8, height = 7, family = "Times")

ggplot(data=cdat_m, aes(x=variable, y=value, color=clust))+
  geom_boxplot(position=position_dodge(width=0.8))+
  labs(x="Solute", y="Proportion of concentration")+  
  scale_color_manual(values=c("1" = colors[1], "2"= colors[2], "3" = colors[3],
                              "4" = colors[4], "5"=colors[5], "6"=colors[6]))+
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle=45, hjust = 1),legend.position="none", text = element_text(size=20))+
  facet_wrap(~clust, nrow=3)

dev.off()

####Random Forest on ESOM Clusters####
dat <- cdat
dat$clust <- factor(dat$clust)
mod.attributes <- c(Rtypes[-4], LCtypes, "MAP")
dat <- dat[complete.cases(dat[,c(mod.attributes,"clust")]),c(mod.attributes, "clust")]

# #split into test and train
set.seed(123)
ind <- sample(2, nrow(dat), replace = TRUE, prob = c(0.7, 0.3))
train <- dat[ind==1,]
test <- dat[ind==2,]

#tune mtry
set.seed(123)
mtry <- tuneRF(y =train$clust,
               x=train[,mod.attributes],stepFactor=2,
               improve=0.01, ntree=500)


table(train$clust)

set.seed(123)
rf_test <- randomForest(clust~., mtry=2, data=train,
                        proximity=F, ntree=1000, sampsize=c(838,838,838,838,838,838), strata=train$clust,
                        replace=T)

pdf("RF_Ntree_Test.pdf", width = 8, height = 6)

plot(rf_test)

dev.off()

set.seed(123)
rf <- randomForest(clust~., mtry=2, data=train,
                   proximity=F, ntree=200, sampsize=c(838,838,838,838,838,838), strata=train$clust,
                   replace=T, importance=T)

rf

p1<-predict(rf, train)

confusionMatrix(p1, train$clust)

p2<-predict(rf, test)

test_clust_count<-data.frame(table(test$clust))
colnames(test_clust_count)<-c("Reference", "Sum")

conf_df<-confusionMatrix(p2, test$clust)
conf_df_table<-data.frame(conf_df$table)
conf_df_table<-merge(conf_df_table, test_clust_count, by="Reference")
conf_df_table$prop<-conf_df_table$Freq/conf_df_table$Sum
conf_df_table$same<-ifelse(conf_df_table$Reference==conf_df_table$Prediction, "yes","no")

#visualize matrix
CM_plot<-ggplot(conf_df_table, aes(Prediction, Reference))+geom_raster(aes(fill=same))+
  scale_fill_manual(values=c("yes"="forestgreen", "no"="salmon"))+
  geom_text(aes(label=round(prop, 2)), size=8, family="Times")+theme_classic()+labs(x="Predicted Cluster",y="Actual Cluster",fill="", tag="a")+
  theme(legend.position = "null", text = element_text(size = 25, family = "Times"))

importance_df<-data.frame(importance(rf))
importance_df$driver<-rownames(importance_df)

vars_order<-importance_df %>%
  dplyr::arrange(desc(MeanDecreaseAccuracy), driver) %>%
  dplyr::select(driver)

importance_df$driver<-factor(importance_df$driver, levels = vars_order$driver)

importance_melt<-melt(importance_df[,-8], id.vars = "driver")

importance_melt$driver<-factor(importance_melt$driver, levels = vars_order$driver)

overall_imp_accuracy<-ggplot(importance_df, aes(MeanDecreaseAccuracy, reorder(driver, MeanDecreaseAccuracy)))+
  geom_bar(stat = "identity", fill="black")+
  theme_classic()+labs(x="Mean Decrease Accuracy (%)", y="", tag="b")+
  theme(text = element_text(size=20, family = "Times"))

overall_imp_gini<-ggplot(importance_df, aes(MeanDecreaseGini, reorder(driver, MeanDecreaseGini)))+
  geom_bar(stat = "identity", fill="black")+
  theme_classic()+labs(x="Mean Decrease Gini Score", y="", tag="b")+
  theme(text = element_text(size=20, family = "Times"))

pdf("RF_performance_gini_accuracy.pdf", width = 14, height = 7)

ggarrange(overall_imp_accuracy, overall_imp_gini)

dev.off()


###plot most important variables across clusters
import_factors<-dat[,c("clust", "MAP", "Carbonates", "GrassShrub", "Wetland", "Forest", "Sandstone",
                       "Mudstone", "Conglomerate", "Sdep")]

import_factors_melt<-melt(import_factors, id.vars = c("clust"))

pdf("MostImportVars_March2024.pdf", width = 14, height = 14, family="Times")

ggplot(import_factors_melt, aes(clust, value))+geom_boxplot(aes(fill=clust), alpha=0.8)+
  facet_wrap(~variable, scales = "free")+theme_classic()+
  theme(text = element_text(size = 20))+labs(x="", y="Driver Value")+
  scale_fill_manual(values=colors)+scale_color_manual(values=colors)+
  theme(axis.text.x = element_blank())+labs(fill="Cluster", col="Cluster")

dev.off()

tiff("Importance_ConfusionMatrix_March2024.tiff", width = 14, height = 7, units = "in", res = 300)

ggarrange(CM_plot, overall_imp_accuracy, widths = c(0.55, 0.45))

dev.off()

