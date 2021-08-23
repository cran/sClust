#' @title Gram similarity matrix checker
#' @author Emilie Poisson Caillault and Erwan Vincent
#' @description Function which select the number of cluster to compute thanks to a selected method
#' @param W Gram Similarity Matrix.
#' @param verbose To output the verbose in the terminal.
#' @return returns a list containing the following elements:
#' \itemize{
#'  \item{Lsym: }{a NJW laplacian matrix}
#'  \item{eigen: }{a list that contain the eigenvectors ans eigenvalues}
#'  \item{diag: }{a diagonal matrix used for the laplacian matrix}
#' }
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
#' res <- compute.laplacian.NJW(W,verbose = TRUE)
#' 
#' ### Example 2: Speed and Stopping Distances of Cars
#' W <- compute.similarity.ZP(scale(cars))
#' W <- checking.gram.similarityMatrix(W)
#' res <- compute.laplacian.NJW(W,verbose = TRUE)
compute.laplacian.NJW <- function(W, verbose= FALSE){
  
  #Calculation of the degree matrix
  if(verbose){message("CALCULATION OF THE DEGREE MATRIX")}
  degrees <- rowSums(W)
  Ds <- diag(1/sqrt(degrees))
  
  #Calculation of the Laplacian matrix
  if(verbose){message("CALCULATION OF THE LAPLACIAN MATRIX")}
  Lsym <- Ds %*% W %*% Ds
  
  #Calculation of the eigen vectors and values
  if(verbose){message("CALCULATION OF THE EIGEN VECTORS AND VALUES")}
  U <- eigen(Lsym, symmetric=TRUE)
  if(verbose){
    message("EIGEN VALUES = ")
    print(U$values)
    message("EIGEN VECTOR = ")
    print(U$vector)
  }
  
  out <- list(Lsym = Lsym, eigen = U, diag = Ds)
}