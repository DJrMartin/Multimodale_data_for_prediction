# Function to transform the spectral signal
library(fda)
# [1]
projRecomp = function(xdata, nBasis, t = 1:dim(xdata)[2], basis = "splines"){
  t = sort((t-min(t))/(max(t)-min(t)))
  if (basis == "fourier") {basisobj = create.fourier.basis(nbasis = nBasis)}
  if (basis == "splines") {basisobj = create.bspline.basis(norder = 3, breaks = seq(head(t,1), tail(t,1), length = nBasis-1))}
  BFunction = getbasismatrix(t, basisobj) 
  Fdata = t(sapply(1:nrow(xdata), 
                   FUN = function(i) t(solve(t(BFunction)
                                             %*% BFunction)%*% t(BFunction) %*% xdata[i,])))
  FdataRec = t(sapply(1:nrow(xdata), 
                      FUN = function(i) t(BFunction %*% solve(t(BFunction)%*% BFunction)%*% 
                                            t(BFunction) %*% xdata[i,])))
  return(list(coeffProj = Fdata, foncRecon = FdataRec, BFunction = BFunction, basisobj = basisobj))
}
# [2]
Spectral_transformation <- function(X, basis = "splines", b = c(7,5,2)){

  p = ncol(X)
  VI = matrix(0,length(b),p)
  x_all = NULL
  cnt = 1
  L = 0
  
  x_recom=NULL
  for (i in b){
    nBasis = round(p/i)
    L = L + nBasis
    x = projRecomp(as.matrix(X), nBasis = nBasis, basis = basis)
    if (cnt==1){
      x_all = x$coeffProj
      x_recom= x$foncRecon
    } else {
      x_all = cbind(x_all,x$coeffProj)
      x_recom= cbind(x_recom, x$foncRecon)
    }
    
    cnt = cnt+1
  
  }
  return(list(data.projection = x_all, data.reconstruction = x_recom))
}

# [3]
read.csv.transform <- function(path, m = 2, transformation = "none", row.names = 1, sep = ","){
  if(transformation != "fourier" & transformation !=  "none" & transformation !=  "splines" & transformation !=  "Z-norm"){
    stop("Please provide a correct transformation : 'none', fourier', 'splines' or 'Z-norm'.")
  }
  # path = "data/data_clean/mir-frotti-mean-2months.csv"
  x <- read.csv(path, row.names = row.names, sep = sep)
  x <- prospectr::savitzkyGolay(x, m = m, p = 3, w = 11)
  if(transformation == "splines"){
    x <- Spectral_transformation(x, basis = "splines", b = c(7,5,2))$data.projection
  }
  if(transformation == "fourier"){
    x <- Spectral_transformation(x, basis="fourier")$data.projection
  }
  if(transformation == "Z-norm"){
    x <- scale(x)
  }
  return(x)
}

