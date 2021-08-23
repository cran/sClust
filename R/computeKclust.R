#' @title Gram similarity matrix checker
#' @author Emilie Poisson Caillault and Erwan Vincent
#' @description Function which select the number of cluster to compute thanks to a selected method
#' @param eigenValues The eigenvalues of the laplacian matrix.
#' @param method The method that will be used. "default" to let the function choose the most suitable method. "PEV" for the Principal EigenValue method. "GAP" for the GAP method.
#' @param Kmax The maximum number of cluster which is allowed.
#' @param tolerence The tolerance allowed for the Principal EigenValue method.
#' @param threshold The threshold to select the dominant eigenvalue for the GAP method.
#' @param verbose To output the verbose in the terminal.
#' @return a vector which contain the number of cluster to compute.
#' @examples
#' ### Example 1: 2 disks of the same size
#' n<-100 ; r1<-1
#' x<-(runif(n)-0.5)*2;
#' y<-(runif(n)-0.5)*2
#' keep1<-which((x*2+y*2)<(r1*2))
#' disk1<-data.frame(x+3*r1,y)[keep1,]
#' disk2 <-data.frame(x-3*r1,y)[keep1,]
#' sameTwoDisks <- rbind(disk1,disk2)
#' W <- compute.similarity.ZP(scale(sameTwoDisks))
#' W <- checking.gram.similarityMatrix(W)
#' eigVal <- compute.laplacian.NJW(W,verbose = TRUE)$eigen$values
#' K <- compute.kclust(eigVal, method="default", Kmax=20, tolerence=0.99, threshold=0.9, verbose=TRUE)
#' 
#' ### Example 2: Speed and Stopping Distances of Cars
#' W <- compute.similarity.ZP(scale(cars))
#' W <- checking.gram.similarityMatrix(W)
#' eigVal <- compute.laplacian.NJW(W,verbose = TRUE)$eigen$values
#' K <- compute.kclust(eigVal, method="default", Kmax=20, tolerence=0.99, threshold=0.9, verbose=TRUE)
compute.kclust <- function(eigenValues, method = "default", Kmax = 20, tolerence =  1,threshold = 0.9,verbose = FALSE){
  
  if(verbose){message("CALCULATION OF THE NUMBER OF CLUSTER")}
  
  K <- sum(eigenValues >= tolerence)
  
  if(method == "default"){
    if(K > 1){
      method = "PEV"
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
