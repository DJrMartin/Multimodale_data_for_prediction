rm(list=ls())
library(yarrr)
library(caret)
load("data/data_preprocess/target_variables.rda")
meta_data = read.csv("data/raw_data/phenotype.csv", sep = ";", header=T)
dir("data/data_preprocess")
load("data/data_preprocess/B-splines-normalisation.rda")
source("functions/boosting_RF_custom.R")

## TEST d'un boosting entre frotti et serum =========================
# Combinaison of blood markers, serum, and metallomic for steatosis.
# Combinaison of serum, and frotti for inflammation.
X_1 = meta_data[,8:17]
X_2 = df$`mir-serum-mean-5months.csv`
X_3 = df$`metallomic-5months.csv`
Y = S$STEATOSIS

w = which(is.na(apply(cbind(X_1, X_2, X_3),1,sum)) == FALSE)
Y = Y[w]
X_1 = X_1[w,] ; X_2 = X_2[w,] ; X_3 = X_3[w,]

SEED = 1
cv = 20
# Condition de base
### Blood markers (permu H0)
set.seed(SEED)
pred = y_test = RMSE0 = c()
for(i in 1:cv){
  train <- createDataPartition(Y, p=0.70, list=F)
  rf <- randomForest(sample(Y[train])~., X_1[train,], ntree = 1000, mtry=3, maxnodes = 10)
  pred = c(pred, predict(rf, X_1[-train,]))
  y_test = c(y_test, Y[-train])
  RMSE0 = c(RMSE0, sqrt(mean((predict(rf, X_1[-train,])-Y[-train])^2)))
}
### Blood markers
set.seed(SEED)
pred = y_test = RMSE1 = c()
for(i in 1:cv){
  train <- createDataPartition(Y, p=0.70, list=F)
  rf <- randomForest(Y[train]~., X_1[train,], ntree = 1500, mtry = 5, maxnodes = 5)
  pred = c(pred, predict(rf, X_1[-train,]))
  y_test = c(y_test, Y[-train])
  RMSE1 = c(RMSE1, sqrt(mean((predict(rf, X_1[-train,])-Y[-train])^2)))
}

# X_2
set.seed(SEED)
pred = y_test = RMSE2 = c()
for(i in 1:cv){
  train <- createDataPartition(Y, p=0.70, list=F)
  rf <- randomForest(Y[train]~., X_2[train,], ntree = 1500, mtry=150, maxnodes = 25)
  pred = c(pred, predict(rf, X_2[-train,]))
  y_test = c(y_test, Y[-train])
  RMSE2 = c(RMSE2, sqrt(mean((predict(rf, X_2[-train,])-Y[-train])^2)))
}
# X_2_with_preselection
set.seed(SEED)
pred = y_test = RMSE2.2 = all.Sel_RF = c()
for(i in 1:cv){
  train <- createDataPartition(Y, p=0.70, list=F)
  SelRF <- which(apply(X_2[train,], 2, function(x) cor.test(x, Y[train], method = "pearson")$p.value)<0.2)

  rf <- randomForest(Y[train]~., X_2[train,SelRF], ntree = 1500, mtry = 150, maxnodes = 25)
  pred = c(pred, predict(rf, X_2[-train,SelRF]))
  y_test = c(y_test, Y[-train])
  RMSE2.2 = c(RMSE2.2, sqrt(mean((predict(rf, X_2[-train,SelRF])-Y[-train])^2)))
}
wilcox.test(RMSE2, RMSE2.2)

# X_3
set.seed(SEED)
pred = y_test = RMSE3 = c()
for(i in 1:cv){
  train <- createDataPartition(Y, p=0.70, list=F)
  rf <- randomForest(Y[train]~., X_3[train,], ntree = 1500, mtry = 3, maxnodes = 5)
  pred = c(pred, predict(rf, X_3[-train,]))
  y_test = c(y_test, Y[-train])
  RMSE3 = c(RMSE3, sqrt(mean((predict(rf, X_3[-train,])-Y[-train])^2)))
}

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

# Boosting with XtgradientBoosting
# $max_depth = 2 ; $eta = 0.001 ; $nthread = 2 ; $objective = "reg:linear" ; $validate_parameters = TRUE
library(xgboost)
set.seed(SEED)
pred = y_test = RMSE4 = c()
for(i in 1:cv){
  w <- createDataPartition(Y, p=0.70, list=F)
  train_x <- cbind(X_1, X_2, X_3)[w,]
  train_y = Y[w]
  test_x <- cbind(X_1, X_2, X_3)[-w,]
  test_y = Y[-w]
  
  #put into the xgb matrix format
  dtrain = xgb.DMatrix(data =  as.matrix(train_x), label = train_y )
  dtest = xgb.DMatrix(data =  as.matrix(test_x), label = test_y)
  # train xgb, evaluating against the validation
  watchlist = list(train = dtrain, valid = dtest)
  
  model = xgb.train(data = dtrain, max.depth = 2, 
                    eta = 0.001, nthread = 2, watchlist = watchlist,
                    nround = 10000, objective = "reg:linear", 
                    early_stopping_rounds = 50, verbose=0)
  
  y_hat_valid = predict(model, dtest)
  
  test_mse = mean(((y_hat_valid - test_y)^2))
  test_rmse = sqrt(test_mse)
  RMSE4 = c(RMSE4, test_rmse)
}

layout(matrix(c(1,2), nrow=1))
par(mar=c(8,4,3,1))
ylim = c(0.6, 2.1)

boxplot(data.frame(RMSE0, RMSE1, RMSE2, RMSE2.2, RMSE3), axes=F, ylab = "RMSE", col = c("grey", "white", "firebrick", "firebrick", "white"), 
        ylim = ylim, main = "Predictive performances from\nSingle model", cex.main = 1)
box()
axis(1, at=1:5, labels=rep("",5))
text(1.3:5.3, 0.4, c("Null distribution","Blood markers","MIR","MIR (SelVar)","Metallome"), 
     xpd = NA, cex = 0.8, srt = 45, pos = 2)
axis(2)
segments(x0 = c(2,2,3), x1 = c(3,4,5), y0= c(1.8, 1.7, 1.9))
text(x = c(2.5, 3, 4), y=c(1.8, 1.7, 1.9)+0.05, "*", cex=1)

par(mar=c(8,1,3,2))
boxplot(data.frame(RMSE1,  res.A$rmse_per_fold, res.B$rmse_per_fold, 
                   res.C$rmse_per_fold, res.D$rmse_per_fold, RMSE4), 
        axes=F, ylab = "RMSE", col = c("white","white", "cornflowerblue", "white", "cornflowerblue","white"),
        ylim = ylim, main = "Predictive performances from\nAggregated models", cex.main = 1)
box()
abline(v = 1.5, lty='dashed')
axis(1, at=1:6, labels=rep("",6))
text(1.3:6.3, 0.4, 
     c("Blood markers","Blood markers\n+ MIR (SelVar)", "Blood markers + metallome", "MIR (SelVar) + metallome", "Blood markers + MIR (SelVar)\n+metallome (multiple models)", "Blood markers + MIR (SelVar)\n+metallome (single model)"),
     xpd = NA, cex = 0.8, srt = 45, pos = 2)
t.test(res.C$rmse_per_fold, RMSE1)

segments(x0 = c(5,1,1,3), x1 = c(6,5,3,6), y0= c(2,1.8, 1.7, 1.9))

text(x = c(5.5), y=c(2.05), "***", cex=1)
text(x = c(2), y=c(1.75), "*", cex=1)
text(c(3, 4.5) , y = c(1.85, 1.95), "**", cex=1)

t.test(RMSE2.2, RMSE3)
t.test(res.D$rmse_per_fold, RMSE1)
t.test(RMSE4, RMSE1)

plot.new()
legend("center", legend = c("Random" , "Neutral" ,"Increase", "Decrease"), 
       fill = c( "grey","white","cornflowerblue", "firebrick"), bty="n", ncol=4)


grep("X1522.", colnames(derivates$`mir-serum-mean-5months.csv`))
plot(S$STEATOSIS, derivates$`mir-serum-mean-5months.csv`[,1198])
