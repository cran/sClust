#' @title Hierarchical Clustering
#' @author Emilie Poisson Caillault and Erwan Vincent
#' @description 
#' @param W Gram Similarity Matrix.
#' @param K number of cluster to obtain.
#' @param method method that will be used in the hierarchical clustering.
#' @param flagDiagZero if True, Put zero on the similarity matrix W.
#' @param verbose To output the verbose in the terminal.
#' @param ... Additional parameter for the hclust function.
#' @return returns a list containing the following elements:
#' \itemize{
#'  \item{cluster: }{a vector containing the cluster}
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
#' res <- HierarchicalClust(W,K=2,method="ward.D2",flagDiagZero=TRUE,verbose=TRUE)
#' plot(sameTwoDisks, col = res$cluster)
#' 
#' ### Example 2: Speed and Stopping Distances of Cars
#' W <- compute.similarity.ZP(scale(iris[,-5]))
#' res <- HierarchicalClust(W,K=2,method="ward.D2",flagDiagZero=TRUE,verbose=TRUE)
#' plot(iris, col = res$cluster)
HierarchicalClust <- function(W,K = 5, method = "ward.D2", flagDiagZero=FALSE, verbose = FALSE, ...){
  
  #Checking the similarity matrix
  W <- checking.gram.similarityMatrix(W, flagDiagZero=flagDiagZero, verbose = verbose)
  
  #Calculation of the graph
  if(verbose){message("CALCULATION OF THE GRAPH")}
  graph <- hclust(as.dist(1-W), method = method,...)
  if(verbose){plot(graph)}
  
  out <- list(cluster = cutree(graph,K))
}