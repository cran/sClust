% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/HierarchicalSC.R
\name{HierarchicalSC}
\alias{HierarchicalSC}
\title{Hierarchical Spectral Clustering}
\usage{
HierarchicalSC(
  W,
  K = 5,
  method = "ward.D2",
  flagDiagZero = FALSE,
  verbose = FALSE
)
}
\arguments{
\item{W}{Gram Similarity Matrix.}

\item{K}{number of cluster to obtain.}

\item{method}{method that will be used in the hierarchical clustering.}

\item{flagDiagZero}{if True, Put zero on the similarity matrix W.}

\item{verbose}{To output the verbose in the terminal.}
}
\value{
returns a list containing the following elements:
\itemize{
 \item{cluster: }{a vector containing the cluster}
 \item{eigenVect: }{a vector containing the eigenvectors}
 \item{eigenVal: }{a vector containing the eigenvalues}
}
}
\description{
Hierarchical Spectral Clustering
}
\examples{
### Example 1: 2 disks of the same size
n<-100 ; r1<-1
x<-(runif(n)-0.5)*2;
y<-(runif(n)-0.5)*2
keep1<-which((x*2+y*2)<(r1*2))
disk1<-data.frame(x+3*r1,y)[keep1,]
disk2 <-data.frame(x-3*r1,y)[keep1,]
sameTwoDisks <- rbind(disk1,disk2)
W <- compute.similarity.ZP(scale(sameTwoDisks))
res <- HierarchicalSC(W,K=2,method = "ward.D2",flagDiagZero=TRUE,verbose=TRUE)
plot(sameTwoDisks, col = res$cluster)
plot(res$eigenVect[,1:2], col = res$cluster, main="spectral space",
     xlim=c(-1,1),ylim=c(-1,1)); points(0,0,pch='+');
plot(res$eigenVal, main="Laplacian eigenvalues",pch='+'); 

### Example 2: Speed and Stopping Distances of Cars
W <- compute.similarity.ZP(scale(iris[,-5]))
res <- HierarchicalSC(W,K=2,method="ward.D2",flagDiagZero=TRUE,verbose=TRUE)
plot(iris, col = res$cluster)
plot(res$eigenVect[,1:2], col = res$cluster, main="spectral space",
     xlim=c(-1,1),ylim=c(-1,1)); points(0,0,pch='+');
plot(res$eigenVal, main="Laplacian eigenvalues",pch='+'); 
}
\references{
Sanchez-Garcia, R., Fernnelly, M. and al. (2014). Hierarchical Spectral Clustering of Power Grids. In IEEE Transaction on Power Systems 29.5, pages 2229-2237. ISSN : 0885-8950.
}
\author{
Emilie Poisson Caillault and Erwan Vincent
}
