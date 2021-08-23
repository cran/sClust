#' @title Fast Spectral Clustering
#' @author Emilie Poisson Caillault and Erwan Vincent
#' @description This function will sample the data before performing a classification function on the samples and then applying K nearest neighbours.
#' @param dataFrame The dataFrame.
#' @param smplPoint maximum of sample number for reduction.
#' @param stopCriteria criterion for minimizing intra-group distance and select final smplPoint.
#' @param neighbours number of points that will be selected for the similarity computation. 
#' @param similarity if True, will use the similarity matrix for the clustering function.
#' @param clustFunction the clustering function to apply on data.
#' @param ... additional arguments for the clustering function.
#' @return returns a list containing the following elements:
#' \itemize{
#'  \item{results: }{clustering results}
#'  \item{sample: }{dataframe containing the sample used}
#'  \item{quantLabels: }{quantization labels}
#'  \item{clustLabels: }{results labels}
#'  \item{kmeans: }{kmeans quantization results}
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
#' res <- fastClustering(scale(sameTwoDisks),smplPoint = 500, 
#'                       stopCriteria = 0.99, neighbours = 7, similarity = TRUE,
#'                       clustFunction = UnormalizedSC, K = 2)
#' plot(sameTwoDisks, col = as.factor(res$clustLabels))
#' 
#' ### Example 2: Speed and Stopping Distances of Cars
#' res <- fastClustering(scale(iris[,-5]),smplPoint = 500, 
#'                       stopCriteria = 0.99, neighbours = 7, similarity = TRUE,
#'                       clustFunction = spectralPAM, K = 3)
#' plot(iris, col = as.factor(res$clustLabels))
#' table(res$clustLabels,iris$Species)
fastClustering <- function(dataFrame, smplPoint, 
                           stopCriteria = 0.99, neighbours = 7, 
                           similarity = TRUE, clustFunction, ...){
  
  ## Quantization part ##
  quantization <- kmeansQuantization(dataFrame, smplPoint, stopCriteria)
  
  ## Growth part ##
  quantSample <- quantization$centers
  pointsAdd <- neighbours
  
  # if 
  if (nrow(quantSample) < trunc(smplPoint/pointsAdd)){
    count<-as.matrix(table(quantization$cluster))
    weight<-count/(min(count))
    Xadd<-min(min(count),pointsAdd)
    M<-data.frame(count=count,weight=trunc(weight*pointsAdd))
    
    nb.points<-apply(X = M, MARGIN = 1, FUN = min)
    total.points<-sum(nb.points)
    
    if(total.points>smplPoint){
      nb.points<-trunc(nb.points/total.points*smplPoint)
    }
    sample<-NULL
    for (i in 1:nrow(quantization$centers)){
      classIndex<-which(quantization$cluster==i)
      n.index<-sample(classIndex,nb.points[i])
      sample<-rbind(sample,quantSample[i,],dataFrame[n.index,])
    }
    
    # # max number of points that can be added
    # count<-as.matrix(table(quantization$cluster))
    # pointsAdd<-min(min(count),pointsAdd)
    # 
    # #Number of points to add in each cluster with their weigth
    # nbPoints <- apply(data.frame(count=count,weight=trunc(count/(min(count))*pointsAdd)),
    #                  MARGIN = 1, FUN = min)
    # totPoints <- sum(nbPoints)
    # 
    # # if there is more points to process than max points to sample
    # if(totPoints > smplPoint){
    #   nbPoints<-trunc(nbPoints / totPoints * smplPoint)
    # }
    # sample <- sapply(1:unique(nrow(quantization$centers)), function(x){
    #   # indexes of points to sample in the cluster 
    #   indexesClustPoints <- which(quantization$cluster == x)
    #   
    #   #Selecting n points by sampling
    #   indexesNPoints <- sample(indexesClustPoints,nbPoints[x])  
    #   
    #   rbind(quantSample[x,],dataFrame[indexesNPoints,])
    # })
    # 
    # sample <- matrix(unlist(sample), ncol = ncol(quantSample))
    
  # else
  }else{
    sample <- quantSample
  }
  
  # Computing clustering
  if(similarity){
    if(neighbours >= nrow(sample) && trunc(nrow(sample)/2) >= 2){
      neighbours <- trunc(nrow(sample)/2)
    }

    Wsample <- compute.similarity.ZP(sample, neighbours)
    Wsample <- checking.gram.similarityMatrix(Wsample)
    clustResults <- clustFunction(W = Wsample, ...)
  }else{
    clustResults <- clustFunction(sample,...)
  }
  
  #Recreating the labels for the complete data
  # completeLabels <- rep(0, length(quantization$cluster))
  # index = 1
  # resultsLabels <- clustResults$cluster
  # 
  # sapply(unique(quantization$cluster),function(x){
  #   allClustPoints <- which(quantization$cluster == x)
  #   completeLabels[allClustPoints] <<- resultsLabels[index]
  #   index <<- index + 1
  # })
  
  
  
  # sapply(unique(quantization$cluster),function(x){
  #   #print(completeLabels[which(quantization$cluster == x)])
  #   #print(resultsLabels[index])
  #   completeLabels[which(quantization$cluster == x)] <<- resultsLabels[index]
  #   #print(completeLabels)
  #   index <<- index + 1
  # })
  
  completeLabels <- class::knn(sample,dataFrame, cl = clustResults$cluster, k = neighbours)
  
  # output
  out <- list(results = clustResults,
              sample = sample,
              quantLabels = quantization$cluster,
              clustLabels = completeLabels,
              kmeans = quantization)
  
  
}