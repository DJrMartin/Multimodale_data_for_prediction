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

SEED = 123
cv = 40
############### STRATEGY 2 ###############################################################################################
# Boosting with XtgradientBoosting
# $max_depth = 2 ; $eta = 0.001 ; $nthread = 2 ; $objective = "reg:linear" ; $validate_parameters = TRUE
library(xgboost)
# Gradient boosting with multimodal data

set.seed(SEED)
rmse.Early.XtremGB = c()

for(i in 1:cv){
  w <- createDataPartition(Y, p=0.70, list=F)
  
  Selvar <- which(apply(X_2[w,], 2, function(x) cor.test(x, Y[w], method = "pearson")$p.value) < 0.3)
  
  train_x <- cbind(X_1, X_2[,Selvar], X_3)[w,]
  train_y = Y[w]
  test_x <- cbind(X_1, X_2[,Selvar], X_3)[-w,]
  test_y = Y[-w]
  
  ### Xtrem gradient boosting ================================
  #put into the xgb matrix format
  dtrain = xgb.DMatrix(data =  as.matrix(train_x), label = train_y )
  dtest = xgb.DMatrix(data =  as.matrix(test_x), label = test_y)
  # train xgb, evaluating against the validation
  watchlist = list(train = dtrain, valid = dtest)
  
  model = xgb.train(data = dtrain, max.depth = 2, 
                    eta = 0.001, nthread = 2, watchlist = watchlist,
                    nround = 10000, objective = "reg:linear", 
                    early_stopping_rounds = 50, verbose=0)
  
  # Compute prdiction
  y_hat_valid = predict(model, dtest)
  # Compute RMSE
  test_mse = mean(((y_hat_valid - test_y)^2))
  test_rmse = sqrt(test_mse)
  
  rmse.Early.XtremGB = c(rmse.Early.XtremGB, test_rmse)
  print(cv)
}

# Random forest with multimodal data
set.seed(SEED)
rmse.Early.RF = c()

for(i in 1:cv){
  # Creation of the partition data
  w <- createDataPartition(Y, p=0.70, list=F)
  # Variables selection
  Selvar <- which(apply(X_2[w,], 2, function(x) cor.test(x, Y[w], method = "pearson")$p.value) < 0.3)
  # training and testing
  train_x <- data.frame(cbind(X_1, X_2[,Selvar], X_3)[w,])
  train_y = Y[w]
  test_x <- data.frame(cbind(X_1, X_2[,Selvar], X_3)[-w,])
  test_y = Y[-w]
  
  ### RF ================================
  model = randomForest(train_y~., train_x)
  
  y_hat_valid = predict(model, test_x)
  
  test_mse = mean(((y_hat_valid - test_y)^2))
  test_rmse = sqrt(test_mse)
  
  rmse.Early.RF = c(rmse.Early.RF, test_rmse)
  
  print(cv)
}

save(rmse.Early.XtremGB, rmse.Early.RF, file="data/results/EarlyIntegration.rda")

