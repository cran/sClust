#' @title Spectral-PAM clustering
#' @author Emilie Poisson Caillault and Erwan Vincent
#' @description The function, for a given similarity matrix, will separate the data using a spectral space.It is based on the Jordan and Weiss algorithm. This version uses K-medoid to split the clusters.
#' @param W Gram Similarity Matrix.
#' @param K number of cluster to obtain.
#' @param flagDiagZero if True, Put zero on the similarity matrix W.
#' @param verbose To output the verbose in the terminal.
#' @return returns a list containing the following elements:
#' \itemize{
#'  \item{cluster: }{a vector containing the cluster}
#'  \item{eigenVect: }{a vector containing the eigenvectors}
#'  \item{eigenVal: }{a vector containing the eigenvalues}
#' }
#' @import cluster
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
#' res <- spectralPAM(W,K=2,flagDiagZero=TRUE,verbose=TRUE)
#' plot(sameTwoDisks, col = res$cluster)
#' plot(res$eigenVect[,1:2], col = res$cluster, main="spectral space",
#'      xlim=c(-1,1),ylim=c(-1,1)); points(0,0,pch='+');
#' plot(res$eigenVal, main="Laplacian eigenvalues",pch='+'); 
#' abline(h=1,lty="dashed",col="red")
#' 
#' ### Example 2: Speed and Stopping Distances of Cars
#' W <- compute.similarity.ZP(scale(iris[-5]))
#' res <- spectralPAM(W,K=2,flagDiagZero=TRUE,verbose=TRUE)
#' plot(iris, col = res$cluster)
#' plot(res$eigenVect[,1:2], col = res$cluster, main="spectral space",
#'      xlim=c(-1,1),ylim=c(-1,1)); points(0,0,pch='+');
#' plot(res$eigenVal, main="Laplacian eigenvalues",pch='+'); 
#' abline(h=1,lty="dashed",col="red")
spectralPAM <- function(W, K, flagDiagZero=FALSE, verbose = FALSE){
  #case K=1
  out<-NULL
  if (K==1) {
    if(verbose){message("K=1 NO clustering")}
    return (list(cluster = rep(1,nrow(W)), eigenVect = NULL, eigenVal = NULL))
  }
  
  #Checking the similarity matrix
  W <- checking.gram.similarityMatrix(W, flagDiagZero=flagDiagZero, verbose = verbose)
  
  #Compute the laplacian
  laplacianResult <- compute.laplacian.NJW(W, verbose = verbose)
  
  #K Eigenvectors selection
  if(verbose){message(paste(K, " EIGEN VECTORS SELECTION"))}
  V <- laplacianResult$eigen$vector[,1:K]
    
  #normalization with sphere projection if more than one cluster
  if(verbose){message("NORMALIZATION OF EIGENVECTOR")}
  Y <- V/apply(V,MARGIN=1,FUN=function(x) norm(matrix(x),"f"))
  if(verbose){
    message("NORMALIZED K EIGEN VECTORS = ")
    print(Y)
  }
  #PAM clustering
  if(verbose){message("PAM CLUSTERING ON THE EIGEN VECTORS")}
  cluster <- cluster::pam(Y,K)$clustering
  if(verbose){cat("cluster = ", cluster, "\n")}
  
  out <- list(cluster = cluster, eigenVect = V, eigenVal = laplacianResult$eigen$values, eigenVectNorm = Y)
  
}



