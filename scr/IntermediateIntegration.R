rm(list=ls())
library(caret)
library(FactoMineR)
library(randomForest)
library("rospca")
library(FactoMineR)

# data importation.
load("data/data_preprocess/target_variables.rda")
load("data/data_preprocess/B-splines-normalisation.rda")
# meta data.
meta_data = read.csv("data/raw_data/phenotype.csv", sep = ";", header=T)

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

# Nombre d’échantillons
n_samples <- nrow(X_1)

# Créer des noms d’échantillons communs
common_names <- paste0("Sample_", 1:n_samples)
rownames(X_1) <- rownames(X_2) <- rownames(X_3) <- common_names

### TEST with MIR and Blood markers. ============================================
rmse_mir.blood = c()
cv = 30
k = 10
for (cv in 1:cv){
  
  intraining = createDataPartition(y, p = 0.7, list = FALSE)
  
  # Variables selections of the MIR.
  SelRF <- which(apply(X_2[intraining,], 2, function(x) cor.test(x, y[intraining], method = "kendall")$p.value) < 0.3)
  # Sets
  X_2.train = X_2[intraining,SelRF]
  X_2.test = X_2[-intraining,SelRF]
  
  # Robust PCA
  res.rpca <- rospca(X_2.train, k = k)
  
  # CCA
  CCA.train = matrix(NA, nrow = length(intraining), ncol = k)
  for(p in 1:k){
    CCA.train[,p] <- lm(res.rpca$scores[,p]~., X_1[intraining,])$fitted
  }
  Exp.variance <- sd(CCA.train)^2/sd(res.rpca$scores)^2
  
  # Projection of the individuals in the test set
  scores.test <- as.matrix(X_2.test) %*% res.rpca$loadings
  
  # CCA on the test set.
  CCA.test = matrix(NA, nrow = n_samples-length(intraining), ncol = k)
  for(p in 1:k){
    CCA.test[,p] <- lm(scores.test[,p]~., X_1[-intraining,])$fitted
  }
  
  ## Algorithm RF
  rf <- randomForest(y[intraining]~., CCA.train)
  
  y_pred <- predict(rf, CCA.test)
  rmse_mir.blood <- c(rmse_mir.blood, sqrt(mean((y_pred-y[-intraining])^2)))

}

### TEST with MIR + Blood markers /  ============================================
rmse_all.integrate = c()
cv = 30
k = 10
for (cv in 1:cv){
  
  intraining = createDataPartition(y, p = 0.7, list = FALSE)
  
  # Variables selections of the MIR.
  SelRF <- which(apply(X_2[intraining,], 2, function(x) cor.test(x, y[intraining], method = "kendall")$p.value) < 0.3)
  # Sets
  X_2.train = X_2[intraining,SelRF]
  X_2.test = X_2[-intraining,SelRF]
  
  # Robust PCA
  res.rpca <- rospca(X_2.train, k = k)
  res.pca.X1 <- PCA(X_1[intraining,], graph = F, ncp = k)
  
  # CCA
  CCA.train.1 = CCA.train.2 = CCA.train.3 = matrix(NA, nrow = length(intraining), ncol = k)
  for(p in 1:k){
    CCA.train.1[,p] <- lm(res.rpca$scores[,p]~., X_1[intraining,])$fitted
    CCA.train.2[,p] <- lm(res.pca.X1$ind$coord[,p]~., X_3[intraining,])$fitted
    CCA.train.3[,p] <- lm(res.rpca$scores[,p]~., X_3[intraining,])$fitted
  }
  
  # Projection of the individuals in the test set
  scores.test.1 <- as.matrix(X_2.test) %*% res.rpca$loadings
  scores.test.2 <- as.matrix(X_1[-intraining,]) %*% res.pca.X1$var$coord
  
  # CCA on the test set.
  CCA.test.1 = CCA.test.2 = CCA.test.3 = matrix(NA, nrow = n_samples-length(intraining), ncol = k)
  
  for(p in 1:k){
    CCA.test.1[,p] <- lm(scores.test.1[,p]~., X_1[-intraining,])$fitted
    CCA.test.2[,p] <- lm(scores.test.2[,p]~., X_3[-intraining,])$fitted
    CCA.test.3[,p] <- lm(scores.test.1[,p]~., X_3[-intraining,])$fitted
  }
  
  df = data.frame(CCA.train.1, CCA.train.2, CCA.train.3)
  df.new = data.frame(CCA.test.1, CCA.test.2, CCA.test.3)
  
  ## Algorithm RF
  rf <- randomForest(y[intraining]~., df)
  y_pred <- predict(rf, df.new)
  rmse_all.integrate <- c(rmse_all.integrate, sqrt(mean((y_pred-y[-intraining])^2)))
}

