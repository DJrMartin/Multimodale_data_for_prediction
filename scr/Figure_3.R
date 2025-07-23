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

############### STRATEGY 2 ###############################################################################################
# Boosting with XtgradientBoosting
# $max_depth = 2 ; $eta = 0.001 ; $nthread = 2 ; $objective = "reg:linear" ; $validate_parameters = TRUE
library(xgboost)
# Gradient boosting with multimodal data
set.seed(SEED)
pred = y_test = RMSE4 = c()
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
  
  y_hat_valid = predict(model, dtest)
  
  test_mse = mean(((y_hat_valid - test_y)^2))
  test_rmse = sqrt(test_mse)
  RMSE4 = c(RMSE4, test_rmse)
}

# Random forest with multimodal data
set.seed(SEED)
pred = y_test = RMSE5 = c()
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
  RMSE5 = c(RMSE5, test_rmse)
}

############### STRATEGY 3 ############################################################################
# Boosting with SelRF.
ntree = rep(1500, 3)
source("functions/boosting_RF_custom.R")

res.A = boosting_rf_cv(response = Y, p = 0.70, datasets = list(X_1, X_2), folds = cv, Selection = c(1, 0.3),
                       ntree = ntree, maxnodes = c(5, 25), mtry = c(3, 150), seed = SEED)
res.B = boosting_rf_cv(response = Y, p = 0.70, datasets = list(X_1, X_3), folds = cv, Selection = c(1, 1),
                       ntree = ntree, maxnodes = c(5, 5), mtry = c(3, 3), seed = SEED)
res.C = boosting_rf_cv(response = Y, p = 0.70, datasets = list(X_2, X_3), folds = cv, Selection = c(0.3, 1),
                       ntree = ntree, maxnodes = c(25, 5), mtry = c(150, 3), seed = SEED)
res.D = boosting_rf_cv(response = Y, p = 0.70, datasets = list(X_1, X_2, X_3), folds = cv, Selection = c(1, 0.3, 1),
                       ntree = ntree, maxnodes = c(5, 25, 5), mtry = c(3, 150, 3), seed = SEED)
res.D.order = boosting_rf_cv(response = Y, p = 0.70, datasets = list(X_1, X_3, X_2), folds = cv, Selection = c(1, 1, 0.3),
                       ntree = ntree, maxnodes = c(5, 5, 25), mtry = c(3, 3, 150), seed = SEED)

pdf(file = "figures/Figure_4.pdf", width = 8.5, height = 4)
############### FIGURE 4A ###################
layout(matrix(c(1,1,2,3,3), nrow=1))
par(mar=c(10,4,3,1))
ylim = c(range(res.A$rmse_per_fold, res.D.order$rmse_per_fold)[1], 2.2)

boxplot(data.frame(RMSE0, RMSE1, RMSE2, RMSE2.2, RMSE3, RMSE_all), axes=F, ylab = "RMSE", 
        col = c("grey", "white", "firebrick", "firebrick", "firebrick", "#E1ECD8"), 
        ylim = ylim, main = "Strategy n°1 - Pooled predictions", cex.main = 1)
box()
axis(1, at=1:6, labels=rep("",6))
text(1.3:6.3, 0.3, c("Null distribution","Serum markers","MIR","MIR (SelVar)","Metallome", "Pooled predictions"), 
     xpd = NA, cex = 0.8, srt = 45, pos = 2)
axis(2)
# segments(x0 = c(2,2,3), x1 = c(3,4,5), y0 = c(1.8, 1.7, 1.9))
# text(x = c(2.5, 3, 4), y = c(1.8, 1.7, 1.9) + 0.05, "*", cex = 1)

############### FIGURE 4B ###################
par(mar=c(10,1,3,2))
boxplot(data.frame(RMSE_all,  RMSE4, RMSE5), 
        axes = F, ylab = "RMSE", col = c("#E1ECD8", "white", "white"),
        ylim = ylim, main = "Strategy n°2 - Pooled data", cex.main = 1)
box()
abline(v = 1.5, lty = 2)
axis(1, at = 1:3, labels = rep("", 3))
text(1.3:3.3, 0.3, c("Serum markers","Gradient boosting","RF"), 
     xpd = NA, cex = 0.8, srt = 45, pos = 2)

############### FIGURE 4C ###################
par(mar=c(10,1,3,2))
boxplot(data.frame(RMSE_all,  res.A$rmse_per_fold, res.B$rmse_per_fold, 
                   res.C$rmse_per_fold, res.D$rmse_per_fold, res.D.order$rmse_per_fold), 
        axes=F, ylab = "RMSE", col = c("#E1ECD8", "cornflowerblue", "cornflowerblue", "white", "cornflowerblue", "cornflowerblue"),
        ylim = ylim, main = "Strategy n°3 - Aggregated models", cex.main = 1)
box()
abline(v = 1.5, lty = 'dashed')
axis(1, at = 1:6, labels = rep("", 6))
text(1.3:6.3, 0.3, 
     c("Serum markers","Serum markers\n+ MIR (SelVar)", "Serum markers + metallome", "MIR (SelVar) + metallome", 
       "Serum markers + MIR (SelVar)\n+ metallome", "Serum markers + metallome\n+ MIR (SelVar)"),
     xpd = NA, cex = 0.8, srt = 45, pos = 2)
# segments(x0 = c(5,1,1,3), x1 = c(6,5,3,6), y0= c(2,1.8, 1.7, 1.9))
# text(x = c(5.5), y=c(2.05), "***", cex=1)
# text(x = c(2), y=c(1.75), "*", cex=1)
# text(c(3, 4.5) , y = c(1.85, 1.95), "**", cex=1)
dev.off()

# t.test(RMSE1, RMSE3)
pdf(file = "figures/legend_Fig_4.pdf", width = 8.5, height = 4.5)
plot.new()
legend("center", legend = c("Random" , "Neutral" ,"Decrease", "Increase"),title = "RMSE compared to RMSE from models fitted on Serum markers", 
       fill = c( "grey","white","cornflowerblue", "firebrick"), bty="n", ncol=4)
dev.off()