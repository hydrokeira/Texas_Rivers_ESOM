#load data
setwd("/Users/keirajohnson/Box Sync/Keira_Johnson/GGM_Paper2Materials/Codes and Data")

load("InputUMatrixData.RData")

#run RMSE code
errors=c()

for (i in 66:100) {
  
  print(i)
  
  set.seed(123)
  
  esom<-esomTrain(X1, InitMethod = "norm_mean_2std", Epochs=i, Key = 1:nrow(X1))
  
  bestmatchesX<-esom$BestMatches[,2]
  
  bestmatchesY<-esom$BestMatches[,3]
  
  closestWeights<-esom$Weights[(bestmatchesX-1)*esom$Columns+bestmatchesY,]
  
  rmse<-sqrt(mean(rowSums((closestWeights - X1)^2)))
  
  errors<-c(errors, rmse)
  
}

write.csv(errors, "ESOM_Epoch_RMSE.csv")

errors_df<-as.data.frame(errors)

pdf("ESOM_Epoch_RMSE.pdf", width = 7, height = 7)

ggplot(errors_df, aes(x=seq(1,100,1), y=errors))+geom_line()+labs(x="Epoch", y="RMSE")+
  theme_classic()+theme(text = element_text(size = 20))
dev.off()

