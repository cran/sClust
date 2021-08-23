# Calcule matrice de similarite
#---------------------------------------------------------------------
#' @title Recherche du nb de cluster par selon le critere du gap 
#' @author Emilie Poisson Caillault v13/10/2015
#' @param val #valeur propre d'une matrice de similarite
#' @param seuil seuil
#' @param fig booleen 
#' @return Kli
#' @import grDevices
compute.nbCluster.gap<-function(val,seuil=0,fig=FALSE){
  val.accepted=val[val>seuil];
  nbVal=length(val);
  gap=abs(diff(val.accepted))
  if(fig) {dev.new(); plot(gap);}
  K=which.max(gap);
}

#---------------------------------------------------------------------
#' @title Calcule matrice de similarite gaussienn
#' @author Emilie Poisson Caillault v13/10/2015
#' @param points matrice pointsxattributs 
#' @param sigma sigma
#' @return mat
#' @import stats
compute.similarity.gaussien<-function(points, sigma){
  e <- dist(points)
  E <- as.matrix(e)
  mat <- exp(-E^2/(2*sigma^2))
  return(mat)
}

#---------------------------------------------------------------------
#' @title Recherche du voisin num id le plus proche
#' @author Emilie Poisson Caillault v13/10/2015
#' @param vdist vecteur de distance du point avec d'autres points
#' @param vois nombre de voisin a selectionner
#' @return id
search.neighboor<-function(vdist, vois){
  ind=rank(vdist,ties.method="first");
  out<-which(ind==vois)
}

#---------------------------------------------------------------------
#' @title Calcule matrice de similarite gaussienne selon Zelnik-Manor et Perona
#' @author Emilie Poisson Caillault v13/10/2015
#' @description sigma local, attention risque matrice non semi-def positive
#' @param points matrice pointsxattributs 
#' @param vois nombre de voisin qui seront selectionnes
#' @return mat
#' @import stats
compute.similarity.ZP<-function(points, vois=7){
  vois=min(vois,nrow(points))
  e <- dist(points)
  E <- as.matrix(e)
  #identification du voisin num vois pour chaque ligne (vois+1 car lui-meme insere)
  sigmaV=apply(E,1,function(x){ind <- NULL; ind<-search.neighboor(x,vois+1); out<-x[ind]+ .Machine$double.eps})
  #matrice des sigma i x sigma j
  sigma=sapply(sigmaV,function(x){x*sigmaV}) 
  mat <- exp(-E^2/sigma)
  return(mat)
}
