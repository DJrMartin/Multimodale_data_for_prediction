rm(list=ls())
library(yarrr)
library(caret)
load("data/data_preprocess/target_variables.rda")
meta_data = read.csv("data/raw_data/phenotype.csv", sep = ";", header=T)
load("data/data_preprocess/b-splines-normalisation.rda")
source("functions/boosting_RF_custom.R")

## TEST d'un boosting entre frotti et serum =========================
# Combinaison of blood markers, serum, and metallomic for steatosis.
# Combinaison of serum, and frotti for inflammation.
X_1 = meta_data[,8:17]
X_2 = df$`mir-serum-mean-5months.csv`
X_3 = df$`metallomic-5months.csv`
Y = as.factor(meta_data$Triglycerides.Hp>40)

w = which(is.na(apply(cbind(X_1, X_2, X_3),1,sum)) == FALSE)
Y = Y[w]
X_1 = X_1[w,] ; X_2 = X_2[w,] ; X_3 = X_3[w,]

SEED = 1
cv = 30
# Condition de base
### Blood markers (permu H0)
set.seed(SEED)
pred = y_test = roc0 = c()
for(i in 1:cv){
  train <- createDataPartition(Y, p=0.70, list=F)
  rf <- randomForest(sample(Y[train])~., X_1[train,], ntree = 1000, mtry=3, maxnodes = 10)
  pred = c(pred, predict(rf, X_1[-train,], type='prob')[,1])
  y_test = c(y_test, Y[-train])
  roc0 = c(roc0, pROC::auc(Y[-train],predict(rf, X_1[-train,], type='prob')[,1]))
}
### Blood markers
set.seed(SEED)
pred = y_test = roc01 = c()
for(i in 1:cv){
  train <- createDataPartition(Y, p=0.70, list=F)
  rf <- randomForest(Y[train]~., X_1[train,], ntree = 1500, mtry = 5, maxnodes = 5)
  pred = c(pred, predict(rf, X_1[-train,], type='prob')[,1])
  y_test = c(y_test, Y[-train])
  roc01 = c(roc01, pROC::auc(Y[-train],predict(rf, X_1[-train,], type='prob')[,1]))
}
# X_2
set.seed(SEED)
pred = y_test = roc02 = c()
for(i in 1:cv){
  train <- createDataPartition(Y, p=0.70, list=F)
  rf <- randomForest(Y[train]~., X_2[train,], ntree = 1500, mtry=150, maxnodes = 25)
  pred = c(pred, predict(rf, X_2[-train,], type='prob')[,1])
  y_test = c(y_test, Y[-train])
  roc02 = c(roc02, pROC::auc(Y[-train],predict(rf, X_2[-train,], type='prob')[,1]))
}
# X_2_with_preselection
set.seed(SEED)
pred = y_test = roc02.2 = all.Sel_RF = c()
for(i in 1:cv){
  train <- createDataPartition(Y, p=0.70, list=F)
  SelRF <- which(apply(X_2[train,], 2, function(x) wilcox.test(x~Y[train], method = "pearson")$p.value)<0.2)
  rf <- randomForest(Y[train]~., X_2[train,SelRF], ntree = 1500, mtry = 150, maxnodes = 25)
  
  pred = c(pred, predict(rf, X_2[-train,SelRF], type='prob')[,1])
  y_test = c(y_test, Y[-train])
  roc02.2 = c(roc02.2, pROC::auc(Y[-train],predict(rf, X_2[-train,SelRF], type='prob')[,1]))
}
# X_3
set.seed(SEED)
pred = y_test = roc03 = c()
for(i in 1:cv){
  train <- createDataPartition(Y, p=0.70, list=F)
  rf <- randomForest(Y[train]~., X_3[train,], ntree = 1500, mtry = 3, maxnodes = 5)
  
  pred = c(pred, predict(rf, X_3[-train,], type='prob')[,1])
  y_test = c(y_test, Y[-train])
  roc03 = c(roc03, pROC::auc(Y[-train],predict(rf, X_3[-train,], type='prob')[,1]))
}
boxplot(roc0, roc01, roc02, roc02.2, roc03)
t.test(roc01, roc03)

# Boosting with SelRF.

ntree = rep(1500, 4)
SelRF <- which(apply(X_2[train,], 2, function(x) cor.test(x, Y[train], method = "pearson")$p.value)<0.2)
res.A = boosting_rf_cv(response = Y, p = 0.70, datasets = list(X_1, X_2[,SelRF]), folds = cv, 
                       ntree = rep(1500, 3), maxnodes = c(5, 25, 5), mtry = c(3, 150, 3), seed = SEED)
res.B = boosting_rf_cv(response = Y, p = 0.70, datasets = list(X_1, X_3), folds = cv, 
                       ntree = rep(1500, 3), maxnodes = c(5, 5), mtry = c(3, 3), seed = SEED)
res.C = boosting_rf_cv(response = Y, p = 0.70, datasets = list(X_2[,SelRF], X_3), folds = cv, 
                       ntree = rep(1500, 3), maxnodes = c(25, 5), mtry = c(150, 3), seed = SEED)
res.D = boosting_rf_cv(response = Y, p = 0.70, datasets = list(X_1, X_2[,SelRF], X_3), folds = cv, 
                       ntree = rep(1500, 3), maxnodes = c(5, 25, 5), mtry = c(3, 150, 3), seed = SEED)

layout(matrix(c(1,2), nrow=1))
par(mar=c(5,4,2,1))
ylim = range(data.frame(RMSE0, res.D$rmse_per_fold))
boxplot(data.frame(RMSE0, RMSE1, RMSE2, RMSE2.2, RMSE3), axes=F, ylab = "RMSE", col = yarrr::piratepal("pony"), 
        ylim = ylim, main = "Predictive performances from\nrandomForests", cex.main = 0.8)
box()
axis(1, at=1:5, labels=rep("",5))
text(1.3:5.3, 0.2, c("Negative control","Blood only","Mir-serum","Mir-serum (SelVar)","Metallomic"), 
     xpd = NA, cex = 0.8, srt = 30, pos = 2)
axis(2)
par(mar=c(5,1,2,2))
boxplot(data.frame(RMSE1,  res.A$rmse_per_fold, res.B$rmse_per_fold, 
                   res.C$rmse_per_fold, res.D$rmse_per_fold ), 
        axes=F, ylab = "RMSE", col = c(yarrr::piratepal("pony")[2], yarrr::piratepal("info")),
        ylim = ylim, main = "Predictive performances from\nboosted randomForests", cex.main = 0.8)
box()
abline(v = 1.5, lty='dashed')
axis(1, at=1:5, labels=rep("",5))
text(1.3:5.3, 0.2, 
     c("Blood only","Blood + mir-serum", "Blood + metallomic", "Mir-serum + metallomic", "Blood + mir-serum\n+ metallomic"),
     xpd = NA, cex = 0.8, srt = 30, pos = 2)
text(c(5),y = 1.5, "**", cex=1)
t.test(RMSE1, res.D$rmse_per_fold)