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
library(cetcolor)
library(lubridate)

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

#for four panel plot
p1<-plotMatrix(umatrix, Toroid=FALSE, BmSize = 2, Clean = TRUE, DrawLegend = F)+ggtitle("(a) ESOM Topography")

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

# wss<-list()
# 
# for (i in 1:20) {
#   print(i)
#   set.seed(123)
#   
#   km.out <- kmeans(weights_clust, centers = i, nstart=50)
#   # Save total within sum of squares to wss variable
#   wss[i] <- km.out$tot.withinss
# }
# 
# pdf("ScreePlot_Clusters.pdf", width = 8, height = 6)
# 
# plot(1:20, wss, type = "b", 
#      xlab = "Number of Clusters", 
#      ylab = "Within groups sum of squares")
# 
# dev.off()

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

write.csv(cdat, "CDAT_ESOM.csv")

cdat2<-inner_join(dat, BestMatches, by="obs")

cdat2<-cdat2 %>%
  group_by(id) %>%
  mutate(mean_year=mean(year(date)), mean_month=mean(month(date)))

cdat2 <- cdat2 %>%
  mutate(mean_year_group = round(mean_year/10)*10)

cdat2$year_cluster<-group_indices(cdat2, mean_year_group)

pdf("Clusters_by_Time.pdf", width = 10, height = 7)

ggplot(cdat, aes(x=month(date), fill=as.character(clust)))+geom_density(aes(group=as.character(clust)), alpha=0.5)+
  theme_classic()+scale_fill_manual(values = colors)+scale_x_continuous(labels=seq(1,12,1), breaks = seq(1,12,1))+
  labs(x="Month of Sample Collection", y="Density", fill="Cluster")+
  theme(text = element_text(size = 20))

ggplot(cdat, aes(x=month(date), fill=as.character(clust)))+geom_density(aes(group=as.character(clust)), alpha=0.5)+
  theme_classic()+scale_fill_manual(values = colors)+scale_x_continuous(labels=seq(1,12,1), breaks = seq(1,12,1))+
  labs(x="Month of Sample Collection", y="Density", fill="Cluster")+facet_wrap(~clust)+
  theme(text = element_text(size = 20), legend.position = "null")

ggplot(cdat, aes(x=year(date), fill=as.character(clust)))+geom_density(aes(group=as.character(clust)), alpha=0.5)+
  theme_classic()+scale_fill_manual(values = colors)+
  labs(x="Year of Sample Collection", y="Density", fill="Cluster")+
  theme(text = element_text(size = 20))

ggplot(cdat, aes(x=year(date), fill=as.character(clust)))+geom_density(aes(group=as.character(clust)), alpha=0.5)+
  theme_classic()+scale_fill_manual(values = colors)+
  labs(x="Year of Sample Collection", y="Density", fill="Cluster")+facet_wrap(~clust)+
  theme(text = element_text(size = 20), legend.position = "null")

dev.off()

#write.csv(cdat, "Data_wClusters.csv")

colors<-c("#E69F00","#56B4E9","#009E73","#0072B2","#D55E00","#CC79A7")

colors_months<-cet_pal(12, "c4s")
show_col(colors_years)

colors_years<-scales::viridis_pal()(8)

pdf("Umatrix_Clusters_Dec2024.pdf", width = 8, height = 5.5)

#plot umatrix with clusters on top
p1<-plotMatrix(umatrix, as.matrix(BestMatches_cluster[,c("obs","row","col")]), 
           Toroid=F, BmSize=4, DrawLegend=T, Clean=T, 
               Cls=BestMatches_cluster$clust, ClsColors = colors)+labs(tag="a")+
  theme(text = element_text(size=20))

dev.off()

#plot umatrix with clusters on top
p2<-plotMatrix(umatrix, as.matrix(BestMatches_cluster[,c("obs","row","col")]), Toroid=F, BmSize=4, DrawLegend=F, Clean=T, 
               Cls=BestMatches_cluster$clust, ClsColors = colors)+ggtitle("(b) Nodes Colored by Cluster")

pdf("Umatrix_Months.pdf", width = 8, height = 5.5)

#plot umatrix with clusters on top
p3<-plotMatrix(umatrix, as.matrix(cdat2[,c("obs","row","col")]), Toroid=F, BmSize=4, DrawLegend=F, Clean=T, 
           Cls=round(cdat2$mean_month), ClsColors = colors_months)+ggtitle("(c) Nodes Colored by Mean Month")

dev.off()

pdf("Umatrix_Years.pdf", width = 8, height = 5.5)

p4<-plotMatrix(umatrix, as.matrix(cdat2[,c("obs","row","col")]), Toroid=F, BmSize=4, DrawLegend=F, Clean=T, 
           Cls=round(cdat2$year_cluster), ClsColors = colors_years)+ggtitle("(d) Nodes Colored by Mean Decade")+
  labs(color="Mean Decade")+
  scale_color_manual(breaks = c(1,2,3,4,5,6,7,8), values = colors_years,
                                                       labels=c("1950", "1960", "1970","1980","1990", "2000", "2010", "2020"))

dev.off()

p3

pdf("ESOM_four_panel_nolegend.pdf", width = 12, height = 9)

ggarrange(p1, p2, p3, p4)

dev.off()

dis = dist(weights_clust)^2

#make silhouette plots for clusters - order largest to smallest
cluster1<-c("2")
cluster2<-c("4")
cluster3<-c("5")
cluster4 <- c("3")
cluster5 <- c("6")
cluster6<- c("1")

#rename clusters so they are ordered largest to smallest
kmeans_cluster$cluster <- ifelse(kmeans_cluster$cluster %in% cluster1,1,
                                    ifelse(kmeans_cluster$cluster %in% cluster2, 2,
                                           ifelse(kmeans_cluster$cluster %in% cluster3,3,
                                                  ifelse(kmeans_cluster$cluster %in% cluster4,4,
                                                         ifelse(kmeans_cluster$cluster %in% cluster5, 5,
                                                                ifelse(kmeans_cluster$cluster %in% cluster6,6, kmeans_cluster$cluster))))))

sil<-silhouette(kmeans_cluster$cluster,dis)

pdf("Sil_Width_Jan2024.pdf", width = 8, height = 4, family = "Times")

p2<-fviz_silhouette(sil)+
  labs(title="", fill="Clusters", col="Clusters", y="Silhouette width")+
  scale_color_manual(values=c("1" = colors[1], "2"= colors[2], "3" = colors[3],
                              "4" = colors[4], "5"=colors[5], "6"=colors[6]))+
  scale_fill_manual(values=c("1" = colors[1], "2"= colors[2], "3" = colors[3],
                             "4" = colors[4], "5"=colors[5], "6"=colors[6]))+
  labs(tag="b")+theme_classic()+theme(text = element_text(size=20), axis.text.x = element_blank(),
                                      axis.ticks.x = element_blank())

dev.off()

pdf("Umatrix_Clusters_Silhouette_Dec2024.pdf", width = 8, height = 10, family = "Times")

ggarrange(p1, p2, nrow = 2)

dev.off()

id.vars <- c("date", "siteno", "lat.x", "long.x", "River.x.x","joinkey", "area",
             "TDS", "Flow", "pH", "temp","MAP","Alkalinity",
             "Aluminum","Bromide","Nitrate","Iron", "Phosphorus",Rtypes[-4], LCtypes, "clust", "obs")
cdat_m <- reshape2::melt(cdat, id.vars=id.vars)
labs <- c(paste("Cluster", seq(1,6,1)))

cdat_m$variable <- as.character(cdat_m$variable)
cdat_m$value<-as.numeric(cdat_m$value)
cdat_m$clust <- factor(cdat_m$clust)

pdf("Solute_Distribution_Dec2024.pdf", width = 8, height = 7, family = "Times")

ggplot(data=cdat_m, aes(x=variable, y=value, color=clust))+
  geom_boxplot(position=position_dodge(width=0.8))+
  labs(x="Solute", y="Proportion of concentration")+  
  scale_color_manual(values=colors)+
  theme_bw() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle=45, hjust = 1),legend.position="none", text = element_text(size=20))+
  facet_wrap(~clust, nrow=3)

dev.off()

