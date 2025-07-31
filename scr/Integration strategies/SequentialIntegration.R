rm(list=ls())
library(yarrr)
library(caret)
library(randomForest)

load("data/data_preprocess/target_variables.rda")
meta_data = read.csv("data/raw_data/phenotype.csv", sep = ";", header=T)
dir("data/data_preprocess")
load("data/data_preprocess/B-splines-normalisation.rda")

SEED = 12
cv = 30
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

ntree = rep(1500, 3)
source("functions/boosting_RF_custom.R")

res.A = boosting_rf_cv(response = Y, p = 0.70, datasets = list(X_1, X_2), folds = cv, Selection = c(1, 0.3),
                       ntree = ntree, maxnodes = c(5, 25), mtry = c(3, 150), seed = SEED)
res.B = boosting_rf_cv(response = Y, p = 0.70, datasets = list(X_1, X_3), folds = cv, Selection = c(1, 1),
                       ntree = ntree, maxnodes = c(5, 5), mtry = c(3, 3), seed = SEED)
res.C = boosting_rf_cv(response = Y, p = 0.70, datasets = list(X_3, X_2), folds = cv, Selection = c(1, 0.3),
                       ntree = ntree, maxnodes = c(5, 25), mtry = c(3, 150), seed = SEED)
res.D = boosting_rf_cv(response = Y, p = 0.70, datasets = list(X_1, X_2, X_3), folds = cv, Selection = c(1, 0.3, 1),
                       ntree = ntree, maxnodes = c(5, 25, 5), mtry = c(3, 150, 3), seed = SEED)
res.D.order = boosting_rf_cv(response = Y, p = 0.70, datasets = list(X_1, X_3, X_2), folds = cv, Selection = c(1, 1, 0.3),
                             ntree = ntree, maxnodes = c(5, 5, 25), mtry = c(3, 3, 150), seed = SEED)

rmse.blood.MIR <- res.A$rmse_per_fold
rmse.blood.trace <- res.B$rmse_per_fold
rmse.trace.MIR <- res.C$rmse_per_fold
rmse.blood.MIR.trace <- res.D$rmse_per_fold
rmse.blood.trace.MIR <- res.D.order$rmse_per_fold

save(rmse.blood.MIR, rmse.blood.trace, rmse.trace.MIR, rmse.blood.MIR.trace, rmse.blood.trace.MIR, 
     file="data/results/SequentialIntegration.rda")
