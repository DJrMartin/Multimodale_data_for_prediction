
representative.mean <- function(x, depth = 0.35){
  ## Initialisation kmeans
  res.K <- kmeans(dist(x), centers = 2)
  condition <- res.K$betweenss/res.K$totss
  ## Iterations of multiple kmeans
  nb = 0
  while(condition > 0.2){
    if(max(table(res.K$cluster))>2){
      nb = nb+1
      if(sum(res.K$cluster==1)/length(res.K$cluster)==0.5){
        condition <- 0
      }else{
        x <- x[res.K$cluster==which.max(table(res.K$cluster)),]
        res.K <- kmeans(dist(x), centers = 2)
        condition <- res.K$betweenss/res.K$totss
      }
    }else{
      x <- x[res.K$cluster==which.max(table(res.K$cluster)),]
      condition <- 0.0000000001
    }
  }
  if(condition == 0){
    spectr.mu <- rbind(apply(x[res.K$cluster==1,], 2, mean),
                       apply(x[res.K$cluster==2,], 2, mean))
    message("Two spectra are representative of the samples.")
  }else{
    spectr.mu <- apply(x, 2, mean)
  }
  return(list(mean.spectra = spectr.mu, iter.nb = nb))
}

