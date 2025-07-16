library(randomForest)
library(Metrics)

boosting_rf_cv <- function(response, datasets, folds = 5, p = 75, ntree = 100, seed = 123, maxnodes, mtry, 
                           Selection = 1) {
  set.seed(seed)
  
  n <- length(response)
  if (any(sapply(datasets, nrow) != n)) stop("All the datasets need to present the same number of lines.")
  if (length(datasets) < 2) stop("At least 2 datasets are required !")
  
  rmse_list <- c()
  res.test <- c()
  all_feature_importances <- list()
  
  for (k in 1:folds) {
    train_idx <- createDataPartition(response, p = p, list=F)
    test_idx <- setdiff(1:n, train_idx)  

    y_train <- response[train_idx]
    y_test <- response[test_idx]
    
    train_sets <- lapply(datasets, function(d) d[train_idx, , drop = FALSE])
    test_sets  <- lapply(datasets, function(d) d[test_idx, , drop = FALSE])
    
    ### === Étape 1 : entraînement modèle de boosting === ###
    rf_models <- list()
    predictions_train <- list()
    feature_importances_fold <- list()
    Selvar = list()
    
    Selvar[[1]] = which(apply(train_sets[[1]], 2, function(x) cor.test(x, y_train, method = "pearson")$p.value) < Selection[1])
    
    # Modèle 1 sur le premier jeu de données
    rf_models[[1]] <- randomForest(x = train_sets[[1]][, Selvar[[1]]], y = y_train, ntree = ntree[1], maxnodes = maxnodes[1], mtry = mtry[1])
    
    pred_train <- predict(rf_models[[1]], train_sets[[1]][, Selvar[[1]]])
    residuals <- y_train - pred_train
    feature_importances_fold[[1]] <- importance(rf_models[[1]])
    
    # Résidus par les autres jeux de données
    for (i in 2:length(datasets)) {
      
      Selvar[[i]] = which(apply(train_sets[[i]], 2, function(x) cor.test(x, y_train, method = "pearson")$p.value) < Selection[i])
      
      rf_models[[i]] <- randomForest(x = train_sets[[i]][, Selvar[[i]]], y = residuals, 
                                     ntree = ntree[i], maxnodes = maxnodes[i], mtry = mtry[i])
      pred_res <- predict(rf_models[[i]], train_sets[[i]][, Selvar[[i]]])
      residuals <- residuals - pred_res
      
      feature_importances_fold[[i]] <- importance(rf_models[[i]])
    }
    
    ### === Étape 2 : prédiction sur le jeu de test === ###
    final_pred <- rep(0, length(test_idx))
    for (i in 1:length(datasets)) {
      final_pred <- final_pred + predict(rf_models[[i]], test_sets[[i]][,Selvar[[i]]])
    }
    
    all_feature_importances[[k]] <- feature_importances_fold
    res.test <- rbind(res.test, cbind(y_test, final_pred))
    rmse_fold <- rmse(y_test, final_pred)
    rmse_list <- c(rmse_list, rmse_fold)
  }
  
  return(list(
    res.test = res.test,
    rmse_per_fold = rmse_list,
    mean_rmse = mean(rmse_list),
    imp = all_feature_importances
  ))
}