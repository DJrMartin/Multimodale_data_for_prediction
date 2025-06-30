rm(list=ls())
library(yarrr)
library(caret)
library(xgboost)
load("data/data_preprocess/target_variables.rda")
meta_data = read.csv("data/raw_data/phenotype.csv", sep = ";", header=T)
load("data/data_preprocess/b-splines-normalisation.rda")
source("functions/boosting_RF_custom.R")

## TEST d'un boosting entre frotti et serum =========================
# Combinaison of blood markers, serum, and metallomic for steatosis.
# Combinaison of serum, and frotti for inflammation.
X_1 = meta_data[,8:17]
X_2 = df$`mir-serum-mean-5months.csv`[,]
X_3 = df$`metallomic-5months.csv`[,]
Y  = S$STEATOSIS

set.seed(123)
w <- createDataPartition(Y, p=0.70, list=F)
train_x <- cbind(X_1, X_2, X_3)[w,]
train_y = Y[w]
test_x <- cbind(X_1, X_2, X_3)[-w,]
test_y = Y[-w]

#put into the xgb matrix format
dtrain = xgb.DMatrix(data =  as.matrix(train_x), label = train_y )
dtest = xgb.DMatrix(data =  as.matrix(test_x), label = test_y)

# these are the datasets the rmse is evaluated for at each iteration
watchlist = list(train=dtrain, test=dtest)
###
# Grid search first principles 
###

max.depths = c(2, 5, 9)
etas = c(0.01, 0.001)

best_params = 0
best_score = 0

count = 1
for( depth in max.depths ){
  for( num in etas){
    
    bst_grid = xgb.train(data = dtrain, 
                         max.depth = depth, 
                         eta=num, 
                         nthread = 2, 
                         nround = 10000, 
                         watchlist = watchlist, 
                         objective = "reg:linear", 
                         early_stopping_rounds = 50, 
                         verbose=0)
    
    if(count == 1){
      best_params = bst_grid$params
      best_score = bst_grid$best_score
      count = count + 1
    }
    else if( bst_grid$best_score < best_score){
      best_params = bst_grid$params
      best_score = bst_grid$best_score
    }
  }
}

best_params
best_score

