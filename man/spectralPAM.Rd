% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spectralPAM.R
\name{spectralPAM}
\alias{spectralPAM}
\title{Spectral-PAM clustering}
\usage{
spectralPAM(W, K, flagDiagZero = FALSE, verbose = FALSE)
}
\arguments{
\item{W}{Gram Similarity Matrix.}

\item{K}{number of cluster to obtain.}

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
The function, for a given similarity matrix, will separate the data using a spectral space.It is based on the Jordan and Weiss algorithm. This version uses K-medoid to split the clusters.
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
res <- spectralPAM(W,K=2,flagDiagZero=TRUE,verbose=TRUE)
plot(sameTwoDisks, col = res$cluster)
plot(res$eigenVect[,1:2], col = res$cluster, main="spectral space",
     xlim=c(-1,1),ylim=c(-1,1)); points(0,0,pch='+');
plot(res$eigenVal, main="Laplacian eigenvalues",pch='+'); 
abline(h=1,lty="dashed",col="red")

### Example 2: Speed and Stopping Distances of Cars
W <- compute.similarity.ZP(scale(iris[-5]))
res <- spectralPAM(W,K=2,flagDiagZero=TRUE,verbose=TRUE)
plot(iris, col = res$cluster)
plot(res$eigenVect[,1:2], col = res$cluster, main="spectral space",
     xlim=c(-1,1),ylim=c(-1,1)); points(0,0,pch='+');
plot(res$eigenVal, main="Laplacian eigenvalues",pch='+'); 
abline(h=1,lty="dashed",col="red")
}
\author{
Emilie Poisson Caillault and Erwan Vincent
}
