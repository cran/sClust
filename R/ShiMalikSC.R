#' @title Bi-parted Spectral Clustering. Shi and Malik.
#' @author Emilie Poisson Caillault and Erwan Vincent
#' @description Bi-parted spectral clustering based on Shi and Malik algorithm, which separates the data into two distinct clusters
#' @param W Gram Similarity Matrix.
#' @param flagDiagZero if True, Put zero on the similarity matrix W.
#' @param verbose To output the verbose in the terminal.
#' @return returns a list containing the following elements:
#' \itemize{
#'  \item{cluster: }{a vector containing the cluster}
#'  \item{eigenVect: }{a vector containing the eigenvectors}
#'  \item{eigenVal: }{a vector containing the eigenvalues}
#' }
#' @references Shi, J and Malik, J. (2000). Normalized cuts and image segmentation. In PAMI, Transactions on Pattern Analysis and Machine Intelligence, pages 888-905
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
#' res <- ShiMalikSC(W,flagDiagZero=TRUE,verbose=FALSE)
#' plot(sameTwoDisks, col = res$cluster)
#' plot(res$eigenVect[,1:2], col = res$cluster, main="spectral space",
#'      xlim=c(-1,1),ylim=c(-1,1)); points(0,0,pch='+');
#' plot(res$eigenVal, main="Laplacian eigenvalues",pch='+'); 
#' 
#' ### Example 2: Speed and Stopping Distances of Cars
#' W <- compute.similarity.ZP(scale(iris[,-5]))
#' res <- ShiMalikSC(W,flagDiagZero=TRUE,verbose=TRUE)
#' plot(iris, col = res$cluster)
#' plot(res$eigenVect[,1:2], col = res$cluster, main="spectral space",
#'      xlim=c(-1,1),ylim=c(-1,1)); points(0,0,pch='+');
#' plot(res$eigenVal, main="Laplacian eigenvalues",pch='+'); 
ShiMalikSC <- function(W, flagDiagZero=FALSE, verbose = FALSE){
  
  #Checking the similarity matrix
  W <- checking.gram.similarityMatrix(W, flagDiagZero=flagDiagZero, verbose = verbose)
  
  #Calculation of the degree matrix
  if(verbose){message("CALCULATION OF THE DEGREE MATRIX")}
  degrees <- rowSums(W)
  Ds <- diag(1/sqrt(degrees))
  
  #Calculation of the Laplacian matrix
  if(verbose){message("CALCULATION OF THE LAPLACIAN MATRIX")}
  L2 <- diag(1,nrow(Ds)) - Ds %*% W %*% Ds
  
  #Calculation of the eigen vectors and values
  if(verbose){message("CALCULATION OF THE EIGEN VECTORS AND VALUES")}
  U <- eigen(L2, symmetric=TRUE)
  if(verbose){
    message("EIGEN VALUES = ")
    print(U$values)
    message("EIGEN VECTOR = ")
    print(U$vector)
  }
  
  #Recuperation of the second smallest eigenvector
  if(verbose){message("GETTING THE SECOND SMALLEST EIGENVECTOR")}
  
  V <- U$vectors[, ncol(L2)-1]
  if(verbose){print(V)}
  print(V)
  
  #Calculation of the weight vector
  if(verbose){message("CALCULATION OF THE WEIGHT VECTOR")}
  g <- as.matrix(sort(Ds %*% V, decreasing = TRUE)) 
  if(verbose){cat("g = ",g,"\n")}
  
  #Calculation of the cluster
  if(verbose){message("CALCULATION OF THE CLUSTER")}
  cluster <- apply(as.matrix(g), MARGIN = 1, FUN = function(x) if(x>=0){x=1}else{x=2})
  if(verbose){print(cluster)}
  
  out <- list(cluster = cluster, eigenVect = U$vector, eigenVal = U$values)
  
}