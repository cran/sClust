#' @title Data quantization
#' @author Emilie Poisson Caillault and Erwan Vincent
#' @description The function use kmeans algorithm to perform data quantization.
#' @param dataFrame The dataFrame.
#' @param maxData maximum of sample number for reduction.
#' @param stopCriteria criterion for minimizing intra-group distance and select final smplPoint.
#' @return kmeans result
kmeansQuantization <- function(dataFrame, maxData, stopCriteria = 0.99){
   
   maxData=min(maxData,nrow(dataFrame))
   
   begin <- 0
   end <- maxData
   middle <- round((begin + end)/2)
   
   while(abs(middle - end) != 0) {
     
     results <- kmeans(dataFrame, centers = middle, iter.max = 200, nstart = 5, algorithm = c("Hartigan-Wong"))
     if(results$ifault == 4){
       results <- kmeans(dataFrame, centers = middle, iter.max = 200, algorithm = c("MacQueen"))
     }
     
     # if variance criteria
     if ((results$betweenss/results$totss) > stopCriteria) {
       end <- middle
       middle <- round((begin + end)/2)
     }else {
       begin <- middle
       middle <- round((begin + end)/2)
     }
     # if there is no point between the middle and the end
     if (abs(middle - end) == 1) {
       end <- middle
     }
     # if there is no point between the middle and the beginning
     if (abs(middle - begin) == 1) {
       begin <- middle
     }
   }
   return(results)
 }