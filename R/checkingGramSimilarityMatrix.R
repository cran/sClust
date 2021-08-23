#' @title Gram similarity matrix checker
#' @author Emilie Poisson Caillault and Erwan Vincent
#' @description Function to check if a similarity matrix is Gram or not
#' @param W Gram Similarity Matrix or not.
#' @param flagDiagZero if True, Put zero on the similarity matrix W.
#' @param verbose To output the verbose in the terminal.
#' @return a Gram similarity matrix
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
#' 
#' ### Example 2: Speed and Stopping Distances of Cars
#' W <- compute.similarity.ZP(scale(cars))
#' W <- checking.gram.similarityMatrix(W)
checking.gram.similarityMatrix <- function(W, flagDiagZero = FALSE, verbose = FALSE){
  #Verification if the similarity matrix is Gram
  if(verbose){message("CHECKING THE SIMILARITY MATRIX")}
  v=eigen(W,symmetric=TRUE)
  if(length(which(v$values<(-0.1)))) {
    if(verbose){message("NON GRAM INITIAL MATRIX")}
    W <- W %*% t(W)
  }
  
  #put zero on W diagonal
  if(flagDiagZero==T){
    if(verbose){message("PUT ZERO ON SIMILARITY DIAGONAL ")}
    diag(W)=0;
  }
  if(verbose){print(W)}
  out <- W
}