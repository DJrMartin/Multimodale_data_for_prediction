rm(list=ls())
library(yarrr)
library(caret)
library(randomForest)

load("data/data_preprocess/target_variables.rda")
meta_data = read.csv("data/raw_data/phenotype.csv", sep = ";", header=T)
dir("data/data_preprocess")
load("data/data_preprocess/B-splines-normalisation.rda")


## TEST d'un boosting entre frotti et serum =========================
# Combinaison of blood markers, serum, and metallomic for steatosis.
# Combinaison of serum, and frotti for inflammation.
X_1 = meta_data[,8:17]
X_2 = df$`mir-serum-mean-5months.csv`
X_3 = df$`metallomic-5months.csv`

# Synthetic variable
Y = S$STEATOSIS

# NA suppression
w = which(is.na(apply(cbind(X_1, X_2, X_3),1,sum)) == FALSE)
Y = Y[w]
X_1 = X_1[w,] ; X_2 = X_2[w,] ; X_3 = X_3[w,]

############### STRATEGY 1 ###############################################################################################
SEED = 123
cv = 30
### permu H0
set.seed(SEED)
RMSE0 = c()
for(i in 1:cv){
  train <- createDataPartition(Y, p=0.70, list=F)
  rf <- randomForest(sample(Y[train])~., X_1[train,], ntree = 1000, mtry=3, maxnodes = 10)
  RMSE0 = c(RMSE0, sqrt(mean((predict(rf, X_1[-train,])-Y[-train])^2)))
}
### Serum markers
set.seed(SEED)
RMSE1 = c()
for(i in 1:cv){
  train <- createDataPartition(Y, p=0.70, list=F)
  rf <- randomForest(Y[train]~., X_1[train,], ntree = 1500, mtry = 5, maxnodes = 5)
  RMSE1 = c(RMSE1, sqrt(mean((predict(rf, X_1[-train,])-Y[-train])^2)))
}
# MIR
set.seed(SEED)
RMSE2 = c()
for(i in 1:cv){
  train <- createDataPartition(Y, p=0.70, list=F)
  rf <- randomForest(Y[train]~., X_2[train,], ntree = 1500, mtry=150, maxnodes = 25)
  RMSE2 = c(RMSE2, sqrt(mean((predict(rf, X_2[-train,])-Y[-train])^2)))
}
# MIR_with_preselection
set.seed(SEED)
RMSE2.2 = all.Sel_RF = c()
for(i in 1:cv){
  train <- createDataPartition(Y, p=0.70, list=F)
  SelRF <- which(apply(X_2[train,], 2, function(x) cor.test(x, Y[train], method = "pearson")$p.value) < 0.3)
  rf <- randomForest(Y[train]~., X_2[train,SelRF], ntree = 1500, mtry = 150, maxnodes = 25)
  RMSE2.2 = c(RMSE2.2, sqrt(mean((predict(rf, X_2[-train,SelRF])-Y[-train])^2)))
}
# Metallome
set.seed(SEED)
RMSE3 = c()
for(i in 1:cv){
  train <- createDataPartition(Y, p=0.70, list=F)
  rf <- randomForest(Y[train]~., X_3[train,], ntree = 1500, mtry = 3, maxnodes = 5)
  RMSE3 = c(RMSE3, sqrt(mean((predict(rf, X_3[-train,])-Y[-train])^2)))
}
# Modèles mutlimodales
set.seed(SEED)
RMSE_all = c()
for(i in 1:cv){
  train <- createDataPartition(Y, p=0.70, list=F)
  # multi_models
  rf_1<- randomForest(Y[train]~., X_1[train,], ntree = 1500, mtry = 3, maxnodes = 5)
  rf_2 <- randomForest(Y[train]~., X_2[train,], ntree = 1500, mtry = 150, maxnodes = 25)
  rf_3 <- randomForest(Y[train]~., X_3[train,], ntree = 1500, mtry = 3, maxnodes = 5)
  # multi_predictions
  pred_1 = predict(rf_1, X_1) ;  pred_2 = predict(rf_2, X_2) ; pred_3 = predict(rf_3, X_3)
  df <- data.frame(pred_1, pred_2, pred_3)
  
  # Integration with a linear model
  model_integration <- lm(Y[train]~.,  df[train,])
  # plot(model_integration$fitted.values, Y[train])
  Y_pred <- predict(model_integration, df[-train,])
  
  RMSE_all = c(RMSE_all, sqrt(mean((Y_pred-Y[-train])^2)))
}
