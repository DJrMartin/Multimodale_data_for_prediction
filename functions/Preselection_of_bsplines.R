#' Selection of most dicriminant pathway
#'
#' @param Splines is an object resulting from MIR2splines function
#' @param labels  target variable
#' @param method lasso (default), ttest
#' @param nb_repeat number of repetitions
#' @param pvalue_threshold a primary selection is done by a ttest with pvalue threshold equal to .2 (default)
#' @param qt a remplir
#' @param balance a remplir
#' @param graph a remplir
#' @param save_res TRUE is results have to be save in a file
#' @param filename name file for saving
#' @param verbose if true some display (default FALSE)
#'
#' @return
#'
#' @importFrom stats cor.test
#' @importFrom Matrix image
#' @importFrom utils str
#' @importFrom openxlsx write.xlsx
#' @importFrom graphics axis abline
#' @importFrom glmnet cv.glmnet
#' @importFrom glmnet predict
#' @import ggplot2
#' @import grDevices
#' @export
#'
#' @examples
SplinesSelection_regression = function(Splines,
                                       labels,
                                       method="lasso",
                                       nb_repeat = 1,
                                       pvalue_threshold = 0.2,
                                       qt=0.9,
                                       balance = FALSE,
                                       graph=FALSE,
                                       save_res = FALSE,
                                       verbose = FALSE,
                                       filename = "donnees.xlsx"
){
  n = length(labels)
  p = dim(Splines$Coefs)[2]
  Splines_sc = scale(Splines$Coefs)
  pds_matrix = (as.data.frame(cbind(Splines_sc , Label = labels))) # scale
  colnames(pds_matrix)[dim(pds_matrix)[2]] = "Label"
  
  # build training set
  SignDif = list()
  selection_ttest = selection = ypred = yobs = min_mse = NULL
  matStud <- matrix(1, p,nb_repeat) # to store p values
  
  for (b in 1:nb_repeat){
    subset = sample(1:n,round(.8*n),replace=TRUE)
    
    
    matMean <- matrix(0, p,1) # to store means
    for (j in 1:(dim(pds_matrix)[2]-1))  # t-test performed on each spectral variable
    {
      if (length(subset)>20){
        test<- cor.test(pds_matrix[subset,j],pds_matrix$Label[subset],methoe="pearson")
      } else {
        test<- cor.test(pds_matrix[subset,j],pds_matrix$Label[subset],methoe="kendall")
      }
      
      matStud[j,b] <- test$p.value
      matMean[j] <- test$estimate
    }
    
    
    ### p value threshold
    select <- which(matStud[,b] < pvalue_threshold)
    selcol <- c(select)
    spbin2 <- Splines$Coefs[,selcol]
    SignDif[[b]] <- cbind(matMean[selcol], matStud[selcol,b])
    colnames(SignDif[[b]]) = c("StatTest","pvalue")
    rownames(SignDif[[b]]) = colnames(Splines$Coefs)[selcol]
    selection_ttest = c(selection_ttest,select)
    
    if (method=="lasso"){
      cv = cv.glmnet(Splines_sc[subset, selcol],
                     labels[subset],
                     alpha=1,
                     nfolds=5)
      label_pred = predict(cv, Splines_sc[-subset, selcol], s = "lambda.1se")
      ypred = c(ypred,label_pred)
      yobs = c(yobs,labels[-subset])
      coef_selected = glmnet::coef.glmnet(cv,s=cv$lambda.1se)
      if (verbose){print(str(coef_selected))}
      select = coef_selected@i[-1]
      selection = c(selection,selcol[select])
      min_mse = c(min_mse,min(cv$cvm))
    } else {selection = selection_ttest}
  }
  
  if (nb_repeat>1){
    tab_select = table(selection)
    
    if (qt<1){
      keep = as.numeric(attributes(which(tab_select/nb_repeat>qt))$names)
      selected = Splines$Coefs[,keep]
    } else {
      keep = as.numeric(attributes(sort(tab_select,decreasing=TRUE)[1:qt])$dimnames[[1]])
      selected = Splines$Coefs[,keep]
    }
    
  } else {
    selected = Splines$Coefs[,selcol]
  }
  
  
  if (graph){
    pMIR = dim(Splines$MIR)[2]
    resolution = Splines$resolution #default
    VI_default = matrix(0,length(resolution),pMIR)
    
    i_deb = 0
    cnt = 1
    for (i in resolution){
      nBasis = round(pMIR/i)
      end = min(nBasis*i,pMIR)
      #     VI_default[cnt,1:end] = c(t(matrix(apply(matStud[i_deb+(1:nBasis),]<pvalue_threshold,1,sum),nBasis,i)))[1:end]
      VI_default[cnt,1:end] = c(t(matrix(apply(1-matStud[i_deb+(1:nBasis),],1,median),nBasis,i)))[1:end]
      i_deb = i_deb+nBasis
      cnt = cnt+1
    }
    importance=t(VI_default)
    rownames(importance)=substr(colnames(Splines$MIR)[1:pMIR],1,4)
    colnames(importance)=resolution
    
    
    
    image(importance,col = hcl.colors(64, "inferno", rev = TRUE),  axes=FALSE)
    axis(2, at=seq(0,1, length=3), labels=colnames(importance), lwd=0, pos=0)
    axis(3, at=seq(0,1, length=11),
         labels=rownames(importance)[round(seq(1,pMIR-1,length=11))], lwd=0,cex.axis=.6)
    abline(v = seq(0,1, length=11),col="gray",lty=3)
    xf = seq(0,1, length=pMIR-1)
    abline(v = xf[which(abs(diff(as.numeric(rownames(importance))))>100)],lty=1,lwd=2,col="green")
  }
  
  if (save_res){ #### =============== TO DO !!! =================== #####
    filesignature <- paste(filename)
    write.xlsx(as.data.frame(selected),
               filesignature,
               colNames = TRUE,
               rowNames = TRUE,
               overwrite = TRUE,
               sheetName="Splines Selection")
  }
  
  return(list(
    selected_splines_coef = selected,
    selected_splines_coef_names = colnames(selected),
    pvalues = matStud,
    importance=importance,
    selection = selcol[selection],
    ypred = ypred,
    yobs = yobs,
    min_mse = min_mse))
}


