#' Calculate the correlation of two vectors
#' 
#' calculate the correlation of two vectors (typically: one dimension chromatography)
#' 
#' 
#' @param target Vector of one chromatography to be aligned
#' @param x Vector of another chromatography, as a reference
#' @return The similarity of the two vectors
#' @author Yizeng Liang ,Zhang Zhimin, Chen Shan
#' @keywords methods
#' @examples
#' 
#' 	data(simulate)
#' 	similarity(p1,p2)

similarity <- function(target,x){
  x <- as.vector(x)
  target <- as.vector(target)
  cx <- x-mean(x)
  ctarget <- target - mean(target)
  sim <- sum(cx*ctarget)/(sqrt(sum(cx^2))*sqrt(sum(ctarget^2)))
  return(sim)
}
