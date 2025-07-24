rm(list=ls())
library(caret)
library(kernlab)
library(quadprog)

load("data/data_preprocess/target_variables.rda")
meta_data = read.csv("data/raw_data/phenotype.csv", sep = ";", header=T)
dir("data/data_preprocess")
load("data/data_preprocess/B-splines-normalisation.rda")
source("functions/Preselection_of_bsplines.R")

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
X_1 = as.matrix(X_1[w,]) ; X_2 = X_2[w,] ; X_3 = X_3[w,]


# --------- 1. Sets of training ---------
set.seed(123)

trainIndex <- createDataPartition(y, p = 0.7, list = FALSE)
SelRF <- which(apply(X_2[trainIndex,], 2, function(x) cor.test(x, y[trainIndex], method = "pearson")$p.value) < 0.3)

X1_train <- X_1[trainIndex, ]; X2_train <- X_2[trainIndex, SelRF]; X3_train <- X_3[trainIndex, ]
X1_test  <- X_1[-trainIndex, ]; X2_test  <- X_2[-trainIndex, SelRF]; X3_test  <- X_3[-trainIndex, ]

y_train <- y[trainIndex]
y_test <- y[-trainIndex]
n_train <- length(y_train)
n_test <- length(y_test)

# --------- 2. Computed the 3 kernel for the training set ---------
# Parameters
sigma <- 1

K1 <- kernelMatrix(rbfdot(sigma = sigma), X1_train)
K2 <- kernelMatrix(rbfdot(sigma = sigma), X2_train)
K3 <- kernelMatrix(rbfdot(sigma = sigma), X3_train)

# Fonction de centrage
center_kernel <- function(K) {
  N <- nrow(K)
  H <- diag(N) - matrix(1, N, N) / N
  return(H %*% K %*% H)
}

# Centrer et normaliser
K1 <- center_kernel(K1); K1 <- K1 / sqrt(sum(K1^2))
K2 <- center_kernel(K2); K2 <- K2 / sqrt(sum(K2^2))
K3 <- center_kernel(K3); K3 <- K3 / sqrt(sum(K3^2))

# Initialisation des poids de noyau
d <- rep(1/3, 3)

# --------- 3. Boucle d'entraînement MKL ---------

for (iter in 1:10) {
  # Combinaison linéaire des noyaux
  K_combined <- d[1] * K1 + d[2] * K2 + d[3] * K3
  
  # SVM entraînement
  model <- ksvm(as.kernelMatrix(K_combined), as.factor(y_train), type = "C-svc", kernel = "matrix", C = 1)
  
  # Extraction des alpha
  alpha <- alpha(model)
  SV_index <- alphaindex(model)[[1]]
  alpha_full <- rep(0, n_train)
  alpha_full[SV_index] <- alpha[[1]] * y_train[SV_index]
  
  # Objectif quadratique pour chaque noyau
  obj <- sapply(list(K1, K2, K3), function(K) sum((alpha_full %*% t(alpha_full)) * K))
  Dmat <- matrix(obj, nrow = 3)
  dvec <- rep(0, 3)
  Amat <- cbind(rep(1, 3), diag(3))
  bvec <- c(1, 0, 0, 0)
  
  # Résolution du QP pour obtenir nouveaux poids
  sol <- solve.QP(Dmat, dvec, Amat, bvec, meq = 1)
  d <- sol$solution
  d <- d / sum(d)
  cat(sprintf("Itération %d : poids transcriptome = %.3f, métabolome = %.3f, microbiote = %.3f\n", iter, d[1], d[2], d[3]))
}

# --------- 4. Noyaux test (validation) ---------

K1_test <- kernelMatrix(rbfdot(sigma = sigma), X1_test, X1_train)
K2_test <- kernelMatrix(rbfdot(sigma = sigma), X2_test, X2_train)
K3_test <- kernelMatrix(rbfdot(sigma = sigma), X3_test, X3_train)

# Centrage test (par rapport au train)
center_test_kernel <- function(K_test, K_train) {
  N <- nrow(K_train)
  mean_train <- colMeans(K_train)
  K_test_centered <- K_test - matrix(mean_train, nrow = nrow(K_test), ncol = N, byrow = TRUE)
  return(K_test_centered)
}

K1_test <- center_test_kernel(K1_test, K1); K1_test <- K1_test / sqrt(sum(K1^2))
K2_test <- center_test_kernel(K2_test, K2); K2_test <- K2_test / sqrt(sum(K2^2))
K3_test <- center_test_kernel(K3_test, K3); K3_test <- K3_test / sqrt(sum(K3^2))

# Combinaison des noyaux test
K_test_combined <- d[1] * K1_test + d[2] * K2_test + d[3] * K3_test

# --------- 5. Prédiction et évaluation ---------

pred <- predict(model, as.kernelMatrix(K_test_combined))

accuracy <- sum(diag(conf_mat)) / sum(conf_mat)
cat(sprintf("Taux de bonne classification (validation) : %.2f%%\n", 100 * accuracy))