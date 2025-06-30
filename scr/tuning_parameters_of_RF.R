rm(list=ls())
library("randomForest")
library("caret")

load("data/data_preprocess/D2-normalisation.rda")
load("data/data_preprocess/target_variables.rda")
source("functions/tuning_parameters_RF.R")

X = derivates$`mir-serum-mean-2months.csv`
Y = S$STEATOSIS

# Determination of the tuning parameters of RFs for blood smear and serum (MIR) at 2 months.  =========================

set.seed(12)
# res.mir.frotti_2M <- tuning_RF(X, Y, nb.cv = 5)
# which.min(apply(res.mir.frotti_2M, 2, mean))
# Parameters for mir-frotti_D2-normalisation : maxnodes = 25, mtry = 75, ntree = 1500.

set.seed(12)
# res.mir.serum_2M <- tuning_RF(X, Y, nb.cv = 5)
# which.min(apply(res.mir.serum_2M, 2, mean))
# Parameters for mir-serum_D2-normalisation : maxnodes = 25, mtry = 150, ntree = 1500.



