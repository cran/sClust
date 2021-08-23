#' @title Spectral Clustering based on the Von Luxburg algorithm
#' @author Emilie Poisson Caillault and Erwan Vincent
#' @description The function, for a given similarity matrix, will separate the data using a spectral space. It uses the Von Luxburg algorithm to do this
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
#' @references Von Luxburg, U. (2007). A Tutorial on Spectral Clustering. Statistics and  Computing,  Volume  17(4), pages 395-416
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
#' res <- VonLuxburgSC(W,K=2,flagDiagZero=TRUE,verbose=TRUE)
#' plot(sameTwoDisks, col = res$cluster)
#' plot(res$eigenVect[,1:2], col = res$cluster, main="spectral space",
#'      xlim=c(-1,1),ylim=c(-1,1)); points(0,0,pch='+');
#' plot(res$eigenVal, main="Laplacian eigenvalues",pch='+'); 
#' 
#' ### Example 2: Speed and Stopping Distances of Cars
#' W <- compute.similarity.ZP(scale(iris[,-5]))
#' res <- VonLuxburgSC(W,K=2,flagDiagZero=TRUE,verbose=TRUE)
#' plot(iris, col = res$cluster)
#' plot(res$eigenVect[,1:2], col = res$cluster, main="spectral space",
#'      xlim=c(-1,1),ylim=c(-1,1)); points(0,0,pch='+');
#' plot(res$eigenVal, main="Laplacian eigenvalues",pch='+'); 
VonLuxburgSC <- function(W, K=5, flagDiagZero=FALSE, verbose = FALSE){
  
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
  
  #K Eigenvectors selection
  if(verbose){message(paste(K, " EIGEN VECTORS SELECTION"))}
  V <- U$vector[,(ncol(U$vector)-K):ncol(U$vector)]
  
  #Normalisation of the eigenvector matrix
  if(verbose){message("CALCULATION OF THE NORMALIZED EIGENVECTOR MATRIX")}
  Fmat <- Ds %*% V
  if(verbose){print(Fmat)}
  
  if(verbose){message("CUTTING THE CLUSTER")}
  cluster <- kmeans(Fmat,K)$cluster
  if(verbose){cat("cluster = ", cluster, "\n")}
  
  out <- list(cluster = cluster, eigenVect = U$vector, eigenVal = U$values)
}