###############################################################
################# TESTING OF THE FUNCTIONS ####################
###############################################################
source("functions/preprocess.R")
## TEST n°1: representative.mean() ####################
mu = c(1,1, rep(4, 2) ,rep(7, 3)) # Situation where the representative spectrum ahve a mean close to 4.
# Generate random variables whit a high level of noise.
data = NULL
for(i in 1:length(mu)){
  data = rbind(data, rnorm(500, mean = mu[i], sd = 3))
}
# Visualisation of the signal
matplot(t(data), type='l')
# Visualisation of the clustering.
plot(hclust(dist(data)))

# The expected result is zero.
# sum(representative.mean(data, 0.4)$mean.spectra - apply(data[mu==4,], 2, mean))

## TEST n°2: boosting_rf_custom() ####################
source("functions/boosting_RF_custom.R")
# Simuler des données
set.seed(42)
n <- 200
X1 <- data.frame(x = rnorm(n))
X2 <- data.frame(z = rnorm(n))
y <- 3 * X1$x + 2 * X2$z + rnorm(n)

# Appliquer la cross-validation
result_cv <- boosting_rf_cv(response = y, datasets = list(X1, X2), folds = 5, ntree = c(100, 100), mtry = c(1, 1), maxnodes = c(5,5))
plot(result_cv$res.test)
abline(0,1)
# Afficher les résultats
print(result_cv$rmse_per_fold)
print(paste("RMSE moyenne :", round(result_cv$mean_rmse, 3)))
