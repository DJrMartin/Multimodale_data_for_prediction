load('data/data_preprocess/b-splines-normalisation.rda')
load("data/data_preprocess/target_variables.rda")
Spec <- read.csv('data/data_clean/mir-serum-mean-5months.csv', sep = ",")

library(randomForest)
library(caret)

df_splines <- data.frame(df$`mir-serum-mean-5months.csv`)
w <- which(is.na(df_splines[,1]))

set.seed(123)
imp = NULL
for(i in 1:30){
  int <- createDataPartition(S$STEATOSIS[-w], p = 0.7, list = F)
  rf <- randomForest(S$STEATOSIS[-w][int]~., df_splines[-w,][int,])
  imp <- cbind(imp, rf$importance)
}

IMP = rowSums(imp)

p = ncol(Spec)

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

image(matrix_plot,col = hcl.colors(200, "Spectral", rev = TRUE),  axes=FALSE,)
axis(2, at=seq(0,1, length=3), labels=colnames(matrix_plot), lwd=0, pos=0)
axis(3, at=seq(0,1, length = length(rownames(matrix_plot))), labels = rownames(matrix_plot), lwd=0, pos=-0.65)
