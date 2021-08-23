#' @title Perform a multi level clustering
#' @author Emilie Poisson Caillault and Erwan Vincent
#' @description The function, for a given dataFrame, will separate the data using the input clustering method in several levels.
#' @param dataFrame The dataFrame.
#' @param levelMax The maximum depth level.
#' @param clustFunction the clustering function to apply on data.
#' @param similarity if True, will use the similarity matrix for the clustering function.
#' @param vois number of points that will be selected for the similarity computation. 
#' @param flagDiagZero if True, Put zero on the similarity matrix W.
#' @param biparted if True, the function will not automatically choose the number of clusters to compute.
#' @param method The method that will be used. "default" to let the function choose the most suitable method. "PEV" for the Principal EigenValue method. "GAP" for the GAP method.
#' @param tolerence The tolerance allowed for the Principal EigenValue method.
#' @param threshold The threshold to select the dominant eigenvalue for the GAP method.
#' @param minPoint The minimum number of points required to compute a cluster.
#' @param verbose To output the verbose in the terminal. 
#' @param ... additional arguments for the clustering function.
#' @return returns a list containing the following elements:
#' \itemize{
#'  \item{cluster: }{vector that contain the result of the last level}
#'  \item{allLevels: }{dataframe containing the clustering results of each levels}
#'  \item{nbLevels: }{the number of computed levels}
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
#' res <- recursClust(scale(sameTwoDisks),levelMax=3, clustFunction =ShiMalikSC,
#'                    similarity = TRUE, vois = 7, flagDiagZero = FALSE,
#'                    biparted = TRUE, verbose = TRUE)
#' plot(sameTwoDisks, col = as.factor(res$cluster))
#' 
#' ### Example 2: Speed and Stopping Distances of Cars
#' res <- recursClust(scale(iris[,-5]),levelMax=4, clustFunction = spectralPAM,
#'                    similarity = TRUE, vois = 7, flagDiagZero = FALSE,
#'                    biparted = FALSE, method = "PEV", tolerence =  0.99,
#'                    threshold = 0.9, verbose = TRUE)
#' plot(iris, col = as.factor(res$cluster))
recursClust <- function(dataFrame, 
                        levelMax = 2, 
                        clustFunction, 
                        similarity = TRUE, vois = 7, flagDiagZero = FALSE,
                        biparted = FALSE, method = "default", 
                        tolerence =  0.99,threshold = 0.9,
                        minPoint = 7,
                        verbose = FALSE,
                        ...){
  stop <- FALSE
  level <- 1
  clusterToCut <- 1
  
  cl <- matrix(1, nrow = nrow(dataFrame), ncol = levelMax)
  
  
  while(!stop){
    newCluster <- c()
    cl[,level+1] <- cl[,level]
    message(level)
    sapply(clusterToCut, FUN = function(x){
      print(x)
      indices = which(cl[,level]==x)
      Xprime <- dataFrame[indices,]
      if(nrow(Xprime) > minPoint){
        #Calculating the similarity
        W <- compute.similarity.ZP(Xprime, vois=vois)
        
        # Checking if the matrix is Gram
        if(nrow(Xprime)<50){
          W <- checking.gram.similarityMatrix(W, flagDiagZero = FALSE, 
                                              verbose = verbose)
        }else{
          W <- checking.gram.similarityMatrix(W, flagDiagZero = flagDiagZero, 
                                              verbose = verbose)
        }
        
        #If not biparted clustering, compute the number of cluster to cut
        if(!biparted){
          eigenValues <- compute.laplacian.NJW(W, verbose = verbose)$eigen$values
          
          kClusters <- compute.kclust2(eigenValues, method = method, Kmax = nrow(Xprime), 
                                       tolerence =  tolerence,threshold = threshold,
                                       verbose = verbose)
          
          #If the clustering method need the similarity matrix in input
          if(similarity){
            results <- clustFunction(W = W, K = kClusters,...)
          #Else put the data in input
          }else{
            results <- clustFunction(dataFrame, K = kClusters,...)
          }
          
        #Else the clustering method as automatically K = 2
        }else{
          #If the clustering method need the similarity matrix in input
          if(similarity){
            results <- clustFunction(W = W,...)
          #Else put the data in input
          }else{
            results <- clustFunction(dataFrame,...)
          }
        }
        
        groups <- results$cluster
        print(groups)
        #Changing the cluster if necessary
        if(!is.null(groups) && length(unique(groups))>1){
          if(level == 1 ){
            cl[indices,level+1] <<- paste0(groups)
            newCluster <<- c(newCluster, unique(paste0(groups)))
          }else{
            cl[indices,level+1] <<- paste0(cl[indices,level],".",groups)
            newCluster <<- c(newCluster, unique(paste0(cl[indices,level],".",groups)))
          }
          if(verbose){
            message("CLUSTERS VECTOR")
            print(cl[,level+1])
          }
          print(newCluster)
          
        }
      }
    })
    
    clusterToCut <- newCluster
    print(clusterToCut)
    stop <- (level+1) >= levelMax
    level <- level + 1
  }
  
  out <- list(
    cluster = cl[,level],
    allLevels = cl,
    nbLevels <- level
  )
}