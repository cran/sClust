#' @title Bi-parted Spectral Clustering. Peronna and Freeman.
#' @author Emilie Poisson Caillault and Erwan Vincent
#' @description Bi-parted spectral clustering based on Peronna and Freeman algorithm, which separates the data into two distinct clusters
#' @param W Gram Similarity Matrix.
#' @param flagDiagZero if True, Put zero on the similarity matrix W.
#' @param verbose To output the verbose in the terminal.
#' @return returns a list containing the following elements:
#' \itemize{
#'  \item{cluster: }{a vector containing the cluster}
#'  \item{eigenVect: }{a vector containing the eigenvectors}
#'  \item{eigenVal: }{a vector containing the eigenvalues}
#' }
#' @references Perona, P. and Freeman, W. (1998). A factorization approach to grouping. In European Conference on Computer Vision, pages 655-670
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
#' res <- PeronaFreemanSC(W,flagDiagZero=TRUE,verbose=TRUE)
#' plot(sameTwoDisks, col = res$cluster)
#' plot(res$eigenVect[,1:2], col = res$cluster, main="spectral space",
#'      xlim=c(-1,1),ylim=c(-1,1)); points(0,0,pch='+');
#' plot(res$eigenVal, main="Laplacian eigenvalues",pch='+'); 
#' 
#' ### Example 2: Speed and Stopping Distances of Cars
#' W <- compute.similarity.ZP(scale(iris[,-5]))
#' res <- PeronaFreemanSC(W,flagDiagZero=TRUE,verbose=TRUE)
#' plot(iris, col = res$cluster)
#' plot(res$eigenVect[,1:2], col = res$cluster, main="spectral space",
#'      xlim=c(-1,1),ylim=c(-1,1)); points(0,0,pch='+');
#' plot(res$eigenVal, main="Laplacian eigenvalues",pch='+'); 
PeronaFreemanSC <- function(W, flagDiagZero=FALSE, verbose = FALSE){
  
  #Checking the similarity matrix
  W <- checking.gram.similarityMatrix(W, flagDiagZero=flagDiagZero, verbose = verbose)
  
  #Calculation of the eigen vectors and values
  if(verbose){message("CALCULATION OF THE EIGEN VECTORS AND VALUES")}
  U <- eigen(W, symmetric=TRUE)
  if(verbose){
    message("EIGEN VALUES = ")
    print(U$values)
    message("EIGEN VECTOR = ")
    print(U$vector)
  }
  
  #Recuperation of the smallest eigenvector
  if(verbose){message("GETTING THE SMALLEST EIGENVECTOR")}
  z1 <- U$vectors[,2]
  if(verbose){cat("z1 = ",z1,"\n")}
  
  #Calculation of lambda1
  if(verbose){message("CALCULATION OF Lambda1")}
  lambda1 <- sqrt(U$value[1])
  if(verbose){cat("lambda1 = ",z1,"\n")}
  
  #Calculation of the weight vector
  if(verbose){message("CALCULATION OF THE WEIGHT VECTOR")}
  p <- as.matrix(sort(lambda1 * z1, decreasing = TRUE)) 
  if(verbose){cat("p = ",p,"\n")}
  
  #Calculation of the cluster
  if(verbose){message("CUTTING THE CLUSTER")}
  cluster <- apply(as.matrix(p), MARGIN = 1, FUN = function(x){y=1;
  if(x<0){y=2};
  out <- y
  })
  if(verbose){cat("cluster = ", cluster, "\n")}
  
  out <- list(cluster = cluster, eigenVect = U$vector, eigenVal = U$values)
}