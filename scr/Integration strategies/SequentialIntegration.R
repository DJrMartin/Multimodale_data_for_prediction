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

SEED = 123
cv = 30
ntree = rep(1500, 3)
source("functions/boosting_RF_custom.R")

res.A = boosting_rf_cv(response = Y, p = 0.70, datasets = list(X_1, X_2), folds = cv, Selection = c(1, 0.3),
                       ntree = ntree, maxnodes = c(5, 25), mtry = c(3, 150), seed = SEED)
res.B = boosting_rf_cv(response = Y, p = 0.70, datasets = list(X_1, X_3), folds = cv, Selection = c(1, 1),
                       ntree = ntree, maxnodes = c(5, 5), mtry = c(3, 3), seed = SEED)
res.C = boosting_rf_cv(response = Y, p = 0.70, datasets = list(X_3, X_2), folds = cv, Selection = c(1, 0.3),
                       ntree = ntree, maxnodes = c(5, 25), mtry = c(3, 150), seed = SEED)
res.D = boosting_rf_cv(response = Y, p = 0.70, datasets = list(X_1, X_2, X_3), folds = cv, Selection = c(1, 1, 1),
                       ntree = ntree, maxnodes = c(5, 25, 5), mtry = c(3, 150, 3), seed = SEED)
res.D.order = boosting_rf_cv(response = Y, p = 0.70, datasets = list(X_1, X_3, X_2), folds = cv, Selection = c(1, 1, 1),
                             ntree = ntree, maxnodes = c(5, 5, 25), mtry = c(3, 3, 150), seed = SEED)
# res.D.order.2 = boosting_rf_cv(response = Y, p = 0.70, datasets = list(X_2, X_3, X_1), folds = cv, Selection = c(1, 1, 1),
#                              ntree = ntree, maxnodes = c(25, 5, 5), mtry = c(150, 3, 3), seed = SEED)

# rmse.blood.MIR <- res.A$rmse_per_fold
# rmse.blood.trace <- res.B$rmse_per_fold
# rmse.trace.MIR <- res.C$rmse_per_fold
# rmse.blood.MIR.trace <- res.D$rmse_per_fold
# rmse.blood.trace.MIR <- res.D.order$rmse_per_fold
# rmse.MIR.trace.blood <- res.D.order.2$rmse_per_fold
# 
# save(rmse.blood.MIR, rmse.blood.trace, rmse.trace.MIR, rmse.blood.MIR.trace, rmse.blood.trace.MIR, rmse.MIR.trace.blood,
#      file="data/results/SequentialIntegration.rda")

#### INTERPRETATION Biomarkers and Metalome ==============================================================================
imp_trace = imp_MIR = imp_markers = NULL
for(i in 1:cv){
  imp_markers <- cbind(imp_markers, res.D.order$imp[[i]][[1]])
  imp_trace <- cbind(imp_trace, res.D.order$imp[[i]][[2]])
  imp_MIR <- cbind(imp_MIR, res.D.order$imp[[i]][[3]])
}

which.max(rowSums(imp_markers))
which.max(rowSums(imp_trace))
which.max(rowSums(imp_MIR))
cor(cbind(X_1[,4], X_2[,1204], X_3[, 14]))
summary(model <- lm(Y~X_1[,4]+X_2[,1204]+X_3[,14]))
plot(model$fitted.values, Y)

#### INTERPRETATION MIR ==============================================================================
Spec <- read.csv('data/data_clean/mir-serum-mean-5months.csv', sep = ",")
p = ncol(Spec)

imp = NULL
for(i in 1:cv){
  imp <- cbind(imp, res.D.order$imp[[i]][[3]])
}

IMP = rowSums(imp)

VI_rf_default = matrix(0,length(c(7,5,3)),p)
i_deb = 0
cnt = 1
for (i in c(7,5,3)){
  nBasis = round(p/i)
  end = min(nBasis*i,p)
  VI_rf_default[cnt,1:end] = c(t(matrix(IMP[i_deb+(1:nBasis)],nBasis,i)))[1:end]
  i_deb = i_deb+nBasis
  cnt = cnt+1
}
matrix_plot = t(VI_rf_default)
rownames(matrix_plot) = substr(colnames(Spec),2,5)
colnames(matrix_plot) = c(7,5,3)

colnames(X_2)[order(IMP)][1:5]

freq = as.numeric(rownames(matrix_plot))
w.freq = which(freq < 1800 & freq > 800)

which(apply(matrix_plot[w.freq,], 1, max)>quantile(apply(matrix_plot[w.freq,], 1, max),0.95))

pdf(file = "figures/wavenumbers_A.pdf", width = 8.5, height = 4)
par(mar = c(4,4,2,2))
image(matrix_plot[w.freq,], col = hcl.colors(200, "Spectral", rev = TRUE),  axes=FALSE, xlab = "Wavenumbers (cm-1)")
axis(2, at=seq(0,1, length=3), labels = colnames(matrix_plot), lwd=0, pos=0)
axis(3, at=seq(0,1, length = length(w.freq)), labels = freq[w.freq], lwd=0, pos=-0.65)
dev.off()

pdf(file = "figures/wavenumbers_B.pdf", width = 8.5, height = 4)
plot(as.numeric(Spec[1,-1][,w.freq]), type = "l", axes = F, xlab = "",
     ylab = "Absorbances")
axis(2)
I.wave <- which(freq[w.freq] == "1231" | freq[w.freq] == "1180" | 
                  freq[w.freq] == "951" | freq[w.freq] == "836")
arrows(x0 = I.wave, y0 = rep(max(Spec[1,-1][,w.freq]), 4), y1 = rep(max(Spec[1,-1][,w.freq])-0.1, 4), length = 0.1, angle = 45)
dev.off()
