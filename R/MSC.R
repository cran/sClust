#' @title Multi-Level Spectral Clustering
#' @author Emilie Poisson Caillault and Erwan Vincent
#' @description The function, for a given dataFrame, will separate the data using the NJW clustering in several levels.
#' @param X The dataFrame.
#' @param levelMax The maximum depth level.
#' @param silMin The minimal silhouette allowed. Below this value, the cluster will be cut again.
#' @param vois number of points that will be selected for the similarity computation. 
#' @param flagDiagZero if True, Put zero on the similarity matrix W.
#' @param method The method that will be used. "default" to let the function choose the most suitable method. "PEV" for the Principal EigenValue method. "GAP" for the GAP method.
#' @param Kmax The maximum number of cluster which is allowed.
#' @param tolerence The tolerance allowed for the Principal EigenValue method.
#' @param threshold The threshold to select the dominant eigenvalue for the GAP method.
#' @param minPoint The minimum number of points required to compute a cluster.
#' @param verbose To output the verbose in the terminal. 
#' @return returns a list containing the following elements:
#' \itemize{
#'  \item{cluster: }{a vector containing the cluster}
#'  \item{eigenVect: }{a vector containing the eigenvectors}
#'  \item{eigenVal: }{a vector containing the eigenvalues}
#' }
#' @references Grassi, K. (2020) Definition multivariee et multi-echelle d'etats environnementaux par Machine Learning : Caracterisation de la dynamique phytoplanctonique.
#' @examples
#' ### Example 1: 2 disks of the same size
#' n<-100 ; r1<-1
#' x<-(runif(n)-0.5)*2;
#' y<-(runif(n)-0.5)*2
#' keep1<-which((x*2+y*2)<(r1*2))
#' disk1<-data.frame(x+3*r1,y)[keep1,]
#' disk2 <-data.frame(x-3*r1,y)[keep1,]
#' sameTwoDisks <- rbind(disk1,disk2)
#' res <- MSC(scale(sameTwoDisks),levelMax=5, silMin=0.7, vois=7, 
#'            flagDiagZero=TRUE, method = "default", Kmax = 20, 
#'            tolerence = 0.99,threshold = 0.7, minPoint = 7, verbose = TRUE)
#' plot(sameTwoDisks, col = as.factor(res[,ncol(res)]))
#' 
#' ### Example 2: Speed and Stopping Distances of Cars
#' res <- MSC(scale(iris[,-5]),levelMax=5, silMin=0.7, vois=7, 
#'            flagDiagZero=TRUE, method = "default", Kmax = 20, 
#'            tolerence = 0.99,threshold = 0.9, minPoint = 7, verbose = TRUE)
#' plot(iris, col = as.factor(res[,ncol(res)]))
#' table(res[,ncol(res)],iris$Species)
MSC <- function(X, levelMax, silMin=0.7,
                   vois=7,
                   flagDiagZero=FALSE,
                   method = "default", Kmax = 20, tolerence = 0.99,threshold = 0.7,
                   minPoint = 7,
                   verbose = FALSE){
  
  #Initialization
  clusterToCut <- 1 ; level <- 1 ; stop <- FALSE
  cl <- matrix(1, nrow = nrow(X), ncol = levelMax)
  Winit <- compute.similarity.ZP(X, vois=vois)
  
  #Level clustering
  while(!stop){
    silLevel <- c() ; newCluster <- c()
    cl[,level+1] <- cl[,level]
    message(level)
    #For each cluster to cut
    sapply(clusterToCut, FUN = function(x){
      indices = which(cl[,level]==x)
      Xprime <- X[indices,]
      
      if(nrow(Xprime) > minPoint){
        W <- compute.similarity.ZP(Xprime, vois=vois)
        
        if(nrow(Xprime)<50){
          W <- checking.gram.similarityMatrix(W, flagDiagZero = FALSE, 
                                              verbose = verbose)
        }else{
          W <- checking.gram.similarityMatrix(W, flagDiagZero = flagDiagZero, 
                                              verbose = verbose)
        }
        
        
        eigenValues <- compute.laplacian.NJW(W, verbose = verbose)$eigen$values
        print(eigenValues[1:10])
        
        kClusters <- compute.kclust(eigenValues, method = method, Kmax = Kmax, 
                                    tolerence =  tolerence,threshold = threshold,
                                    verbose = verbose)
        cat("K = ",kClusters,"\n")
        
        groups <- spectralPAM(W, kClusters, verbose = verbose)$cluster
        
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
          
        }
      }
    })
    
    
    clusterToCut <- c()
    
    #Calculating the silhouette
    cluster <- as.numeric(as.factor(cl[,level+1]))
    uniqueCluster <- unique(cluster)
    silLevel <- sapply(uniqueCluster, FUN = function(x){
      mean(cluster::silhouette(cluster,as.dist(1-Winit))[,"sil_width"][which(cluster == x)])
    })
    print(silLevel)
    clusterToCut <- newCluster[which(newCluster %in% unique(cl[,level+1])[which(silLevel<silMin)])]
    print(clusterToCut)
    
    stop <- ((level+1) >= levelMax) || (length(clusterToCut)==0)
    level <- level + 1
  }
  out <- cl[,1:level]
}

