rm(list=ls())
library(caret)
library(FactoMineR)
library(randomForest)
library(rospca)

# data importation.
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
y = Y[w]
X_1 = data.frame(X_1[w,]) ; X_2 = data.frame(X_2[w,]) ; X_3 = data.frame(X_3[w,])

SEED = 123

# Nombre d’échantillons
n_samples <- nrow(X_1)

# Créer des noms d’échantillons communs
common_names <- paste0("Sample_", 1:n_samples)
rownames(X_1) <- rownames(X_2) <- rownames(X_3) <- common_names

### TEST with MIR and Blood markers.
rmse.MixedR = c()
cv = 40
set.seed(SEED)
for (cv in 1:cv){
  intraining = createDataPartition(y, p = 0.7, list = FALSE)
  
  # Variables selections of the MIR.
  SelRF <- which(apply(X_2[intraining,], 2, function(x) cor.test(x, y[intraining], method = "kendall")$p.value) < 0.3)
  # Sets
  X_2.train = X_2[intraining,SelRF]
  X_2.test = X_2[-intraining,SelRF]
  
  # Robust PCA
  rpca1 <- rospca(X_1[intraining,], k = 10) 
  rpca2 <- rospca(X_2.train, k = 10)
  rpca3 <- rospca(X_3[intraining,], k = 10) 
  
  # Score of the train set
  res.rpca1 = rpca1$scores
  res.rpca2 = rpca2$scores
  res.rpca3 = rpca3$scores
  rpca.train = data.frame(res.rpca1, res.rpca2, res.rpca3)
  
  # Score projection of the test set
  res.rpca1_test = as.matrix(X_1[-intraining,]) %*% rpca1$loadings
  res.rpca2_test = as.matrix(X_2.test) %*% rpca2$loadings
  res.rpca3_test = as.matrix(X_3[-intraining,]) %*% rpca3$loadings
  rpca.test = data.frame(res.rpca1_test, res.rpca2_test, res.rpca3_test)
  
  ## Algorithm RF
  rf <- randomForest(y[intraining]~., rpca.train)
  ## Prediction
  y_pred <- predict(rf, rpca.test)
  
  # Compute the RMSE.
  rmse.MixedR <- c(rmse.MixedR, sqrt(mean((y_pred-y[-intraining])^2)))
  print(cv)

}

boxplot(rmse.MixedR, rmse.MixedK)
# save(rmse.MixedR, rmse.MixedK, file = 'data/res_ML.rda')
