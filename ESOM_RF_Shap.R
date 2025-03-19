require(lubridate)
require(dplyr)
require(reshape2)

setwd("/Users/keirajohnson/Box Sync/Keira_Johnson/GGM_Paper2Materials/Codes and Data")

lulc<-read.csv("LULC_filled_interpolated_joinkey_simplified.csv")

param <- c("Bicarbonate","Calcium", "Chloride","Potassium","Magnesium", "Sodium","Sulfate", "Silica")
Rtypes <- c("Sdep", "Conglomerate", "Carbonates", "Water", "Evaporites",
            "Mudstone", "Igmet", "Sandstone")
LCtypes <- colnames(lulc)[4:10]
attributes <- c("date", "siteno","River.x.x","joinkey","lat.x","long.x",
                "area","TDS","Flow", "pH", "temp",
                param, Rtypes)
id.vars <- c("date", "siteno", "lat.x", "long.x", "River.x.x","joinkey", "area",
             "TDS", "Flow", "pH", "temp")

mod.attributes <- c(Rtypes[-4], LCtypes, "MAP")

shap_df<-read.csv("Shap_values_ESOM_Texas.csv")

dat <- read.csv("CDAT_ESOM.csv")

dat_clust<-dat[,c(2,5,43)]

shap_df<-left_join(shap_df, dat_clust)

shap_df<-subset(shap_df, shap_df$class==shap_df$clust)

shap_mean<-shap_df %>%
  dplyr::group_by(feature) %>%
  dplyr::summarise(mean_abs_shap=mean(abs(phi), na.rm=T))

shap_mean_simple<-shap_df %>%
  dplyr::group_by(feature) %>%
  dplyr::summarise(mean_abs_shap=mean(abs(phi), na.rm=T))

features_df<-data.frame(c(unique(shap_df$feature)),
                        c("lithology", "lithology", "lithology", "lithology", "lithology", "lithology",
                          "lithology", "LULC", "LULC", "LULC", "LULC", "LULC", "LULC", "LULC", "climate"))

colnames(features_df)<-c("feature", "feature_class")

shap_mean<-left_join(shap_mean, features_df)

feature_pal=c("climate"="dodgerblue", "hydromod"="goldenrod", "lithology"="grey", "LULC"="forestgreen", "subsurface"="#ffb7a1",
              "topography"="#e85d04")

pdf("Shap_mean_ESOM.pdf", width = 16, height = 10)

SHAP_overall<-ggplot(shap_mean, aes(mean_abs_shap, feature, fill=feature_class))+geom_bar(stat="identity")+theme_classic()+
  scale_y_discrete(limits=unique(shap_mean_simple$feature[order(shap_mean_simple$mean_abs_shap)]))+
  theme(text = element_text(size=20, family = "Times"),legend.position = "bottom")+
  scale_fill_manual(values = feature_pal)+
  labs(x="Mean Absolute SHAP value", y="", fill="", tag = "b")

ggplot(shap_mean, aes(mean_abs_shap, as.character(class), fill=class))+geom_bar(stat="identity")+theme_classic()+
  #scale_y_discrete(limits=unique(shap_mean_simple$feature[order(shap_mean_simple$mean_abs_shap)]))+
  theme(text = element_text(size=20),legend.position = "bottom")+
  #scale_fill_manual(values = feature_pal)+
  labs(x="Mean Absolute SHAP value", y="", fill="")+
  facet_wrap(~feature)

ggplot(shap_mean, aes(mean_abs_shap, feature, fill=feature_class))+geom_bar(stat="identity")+theme_classic()+
  scale_y_discrete(limits=unique(shap_mean_simple$feature[order(shap_mean_simple$mean_abs_shap)]))+
  theme(text = element_text(size=20),legend.position = "bottom")+
  scale_fill_manual(values = feature_pal)+
  labs(x="Mean Absolute SHAP value", y="", fill="")+
  facet_wrap(~class)

dev.off()

dat <- read.csv("CDAT_ESOM.csv")
dat$Year<-year(dat$date)
dat<-dat[,-c(35:41)]

lulc<-read.csv("LULC_filled_interpolated_joinkey_simplified.csv")
dat<-left_join(dat, lulc[2:10])
dat$clust <- factor(dat$clust)

dat_sites <- dat[complete.cases(dat[,c(mod.attributes,"clust")]),]

dat<-dat[,c("joinkey","date", "clust", mod.attributes)]

dat_melt<-melt(dat, id.vars=c("joinkey","date", "clust"))

pdf("Feature_Vals_by_Cluster.pdf", width = 14, height = 10, family = "Times")

ggplot(dat_melt, aes(clust, value, fill=clust))+geom_boxplot(outliers = F)+theme_classic()+
  scale_fill_manual(values = color_pal)+
  facet_wrap(~variable, scales = "free")+
  labs(x="", y="Feature Value", fill="Clusters")+
  theme(text = element_text(size=20))

dev.off()

colnames(dat_melt)[4]<-"feature"

vals<-left_join(shap_df[,c(2,3,4,7,8)], dat_melt)

all_features<-unique(vals$feature)

color_pal<-c("#E69F00","#56B4E9","#009E73","#0072B2","#D55E00","#CC79A7")

pdf("Feature_Shap_plot_ESOM_V2.pdf", width = 10, height = 8)

for (i in 1:length(all_features)) {
  
  print(i)
  
  one_feature<-vals %>%
    dplyr::filter(feature==all_features[i])
  
  p1<-ggplot(one_feature, aes(value, phi, col=as.character(class)))+geom_point(alpha=0.6)+theme_classic()+
    geom_abline(y=0, slope = 0, col="red")+ggtitle(all_features[i])+
    labs(y=paste0("SHAP of ", all_features[i]), x=all_features[i], col="Cluster")+theme(text = element_text(size=20))+
    scale_color_manual(values = color_pal)
  
  print(p1)
  
  
}

dev.off()


normalize <- function(x, na.rm = TRUE) {
  return((x- min(x)) /(max(x)-min(x)))
}

vals<-vals %>%
  group_by(feature) %>%
  mutate(value_norm=normalize(value))

unique(vals$feature)

vals<-vals %>%
  mutate(feature=case_when(
    feature=="Shrub_grass"~"Shrub & Grass",
    feature=="Swamp_marsh"~"Swamp & Marsh",
    feature=="Igmet"~"Igneous & Metamor.",
    feature=="Barren_sparse"~"Barren or Sparse",
    feature=="Sdep"~"Sedimentary Dep.",
    .default = feature
  ))

vals$feature<-factor(vals$feature, levels = c("MAP", "Swamp & Marsh", "Evaporites", "Mudstone", "Carbonates", "Cropland",
                                              "Forest", "Sedimentary Dep.","Water","Igneous & Metamor.", "Sandstone",
                                              "Impervious","Barren or Sparse", "Conglomerate", "Shrub & Grass"))

pdf("SHAP_dot_plot_updatedAxis.pdf", width = 10, height = 10)

ggplot(vals, aes(phi, feature, col=value_norm))+geom_point()+facet_wrap(~clust)+
  geom_vline(xintercept = 0, col="grey")+theme_classic()+
  scale_color_gradientn(colors = c("#D7191C", "#FDAE61","goldenrod","#ABD9E9", "#2C7BB6"))+
  theme(text = element_text(size = 16), legend.position = "bottom", legend.text = element_text(angle = 45, hjust=1))+
  labs(x="SHAP value", y="", col="Normalized Feature Value")+scale_y_discrete(limits=rev(levels(vals$feature)))

dev.off()

