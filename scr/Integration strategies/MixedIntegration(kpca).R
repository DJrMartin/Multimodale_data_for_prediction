rm(list=ls())
library(caret)
library(FactoMineR)
library(randomForest)
library(kernlab)

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
rmse.MixedK = c()
cv = 40
k = 50
set.seed(SEED)
for (cv in 1:cv){
  
  intraining = createDataPartition(y, p = 0.7, list = FALSE)
  
  # Variables selections of the MIR.
  SelRF <- which(apply(X_2[intraining,], 2, function(x) cor.test(x, y[intraining], method = "kendall")$p.value) < 0.3)
  # Sets
  X_2.train = X_2[intraining,SelRF]
  X_2.test = X_2[-intraining,SelRF]
  
  # Kernel PCA
  kpca1 <- kpca(~.,X_1[intraining,], 
                    kernel = "rbfdot",
                    kpar = list(sigma = 1e-4),
                    features = 10, th = 1e-4) 
  kpca2 <- kpca(~.,X_2.train, 
                    kernel = "rbfdot",
                    kpar = list(sigma = 100),
                    features = 10, th = 1e-4) 
  kpca3 <- kpca(~.,X_3[intraining,], 
                    kernel = "rbfdot",
                    kpar = list(sigma = 0.1),
                    features = 10, th = 1e-4) 
  
  res.kpca1 = pcv(kpca1)
  res.kpca2 = pcv(kpca2)
  res.kpca3 = pcv(kpca3)
  kpca.train = cbind(res.kpca1,res.kpca2,res.kpca3)
  
  res.kpca1_test = predict(kpca1, X_1[-intraining,])
  res.kpca2_test = predict(kpca2, X_2.test)
  res.kpca3_test = predict(kpca3, X_3[-intraining,])
  kpca.test = cbind(res.kpca1_test,res.kpca2_test, res.kpca3_test)
  
  ## Algorithm RF
  rf <- randomForest(y[intraining]~., kpca.train)
  
  y_pred <- predict(rf, kpca.test)
  rmse.MixedK <- c(rmse.MixedK, sqrt(mean((y_pred-y[-intraining])^2)))
  
  print(cv)
}

boxplot(rmse.MixedK)
