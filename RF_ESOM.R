####Random Forest on ESOM Clusters####
setwd("C:/Users/johnkeir/Box/Keira_Johnson/GGM_Paper2Materials/Codes and Data")
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

dat<- read.csv("CDAT_ESOM.csv")
dat$Year<-year(dat$date)
dat<-dat[,-c(35:41)]

lulc<-read.csv("LULC_filled_interpolated_joinkey_simplified.csv")
dat<-left_join(dat, lulc[2:10])

dat_sites <- dat[complete.cases(dat[,c(mod.attributes,"clust")]),]

dat <- dat[complete.cases(dat[,c(mod.attributes,"clust")]),c(mod.attributes, "clust")]

# write.csv(dat, "RF_input_data_ESOM.csv")
# 
# dat<-read.csv("RF_input_data_ESOM.csv")
# dat<-dat[,-1]
dat$clust <- factor(dat$clust)

# #split into test and train
set.seed(123)
ind <- sample(2, nrow(dat), replace = TRUE, prob = c(0.7, 0.3))
train <- dat[ind==1,]
test <- dat[ind==2,]

train_sites<-dat_sites[ind==1,]

#tune mtry
set.seed(123)
mtry <- tuneRF(y =train$clust,
               x=train[,mod.attributes],stepFactor=2,
               improve=0.01, ntree=500)


table(train$clust)

set.seed(123)
rf_test <- randomForest(clust~., mtry=1, data=train,
                        proximity=F, ntree=1000, sampsize=c(838,838,838,838,838,838), strata=train$clust,
                        replace=T)

pdf("RF_Ntree_Test.pdf", width = 8, height = 6)

plot(rf_test)

dev.off()

set.seed(123)
rf <- randomForest(clust~., mtry=1, data=train,
                   proximity=F, ntree=250, sampsize=c(838,838,838,838,838,838), strata=train$clust,
                   replace=T, importance=T)

rf

predictor <- Predictor$new(rf, data = train[, -16], y = train$clust)

#i=5

site_no_list<-train_sites$joinkey
wy_list<-train_sites$date

shap_list<-list()

for (i in 1:nrow(train)) {
  
  print(i)
  
  set.seed(123)
  shapley <- iml::Shapley$new(predictor, x.interest = train[i,-16])
  
  # Extract SHAP values
  shap_values <- shapley$results
  shap_values$joinkey<-site_no_list[i]
  shap_values$date<-wy_list[i]
  
  shap_list[[i]]<-shap_values
  
}

shap_df<-do.call(bind_rows, shap_list)

write.csv(shap_df, "Shap_values_ESOM_Texas.csv")


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
  theme_classic()+labs(x="Mean Decrease Gini Score", y="")+
  theme(text = element_text(size=20, family = "Times"))



pdf("RF_performance_gini_accuracy.pdf", width = 14, height = 7)

ggarrange(overall_imp_accuracy, overall_imp_gini)

dev.off()

pdf("RF_gini.pdf", width = 7, height = 7)

overall_imp_gini

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

