#' @title K clust compute selection V2
#' @author Emilie Poisson Caillault and Erwan Vincent
#' @description Function which select the number of cluster to compute thanks to a selected method
#' @param eigenValues The eigenvalues of the laplacian matrix.
#' @param method The method that will be used. "default" to let the function choose the most suitable method. "PEV" for the Principal EigenValue method. "GAP" for the GAP method.
#' @param Kmax The maximum number of cluster which is allowed.
#' @param tolerence The tolerance allowed for the Principal EigenValue method.
#' @param threshold The threshold to select the dominant eigenvalue for the GAP method.
#' @param verbose To output the verbose in the terminal.
#' @return a vector which contain the number of cluster to compute.
compute.kclust2 <- function(eigenValues, method = "default", Kmax = 20, tolerence =  1,threshold = 0.9,verbose = FALSE){
  
  if(verbose){message("CALCULATION OF THE NUMBER OF CLUSTER")}
  
  eigValOne = sum(floor(eigenValues+ 1e-3) == 1) 
  if(eigValOne >= 2){
    K = eigValOne
  }else{
    K <- sum(eigenValues >= tolerence)
  }
  
  if(method == "default" || method =="default2"){
    if(K > 1){
      method = "PEV"
    }else if(method =="default2"){
      K = 2
    }else{
      method = "GAP"
    }
  }
  
  if(verbose){message(method)}
  message(method)
  if(method == "GAP"){
    selectedEigVal <- eigenValues[eigenValues>=threshold]
    if(length(selectedEigVal) == 1){
      K = 1
    }else if(length(selectedEigVal)<3){
      K = 2
    }else{
      ecart <- diff(selectedEigVal)
      if(verbose){plot(ecart)}
      K <- which.max(abs(ecart))
    }
  }
  
  if(K > Kmax){K <- Kmax}
  if(K == 0){K <- 1}
  
  out <- K
}