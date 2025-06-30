tuning_RF <- function(X, Y, p=0.8, maxnodes = c(5, 25, 50), 
                      mtry = c(20, 75, 150), ntree = c(500, 1500, 3000), nb.cv=10){
  
  w = which(is.na(X[,1]) == FALSE)
  X = X[w,]
  Y = Y[w]
  ## Tuning parameters.
  
  training <- createDataPartition(Y, p = 0.8, times = nb.cv)
  
  res = matrix(NA, nrow = nb.cv, ncol = length(maxnodes)*length(mtry)*length(ntree))
  
  for(i in seq(length(training))){
    
    y_train = as.numeric(Y[training[[i]]])
    train = data.frame(X[training[[i]],])
    test = data.frame(X[-training[[i]],])
    
    RMSE = NULL
    condition = NULL
    
    for(a in maxnodes){
      for(b in mtry){
        for(c in ntree){
          rf <- randomForest(y_train~., data=train, maxnodes = a, mtry = b, ntree = c)
          RMSE = c(RMSE, sqrt(mean((predict(rf, test)-Y[-training[[i]]])^2)))
          condition = c(condition , paste0(a,"_",b,"_",c))
        }
      }
    }
    res[i,] = RMSE
  }
  colnames(res) = condition
  return(model.results.RMSE = res)
}

