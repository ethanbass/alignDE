#' Whittaker Smoother
#'
#' penalized least squares algorithm for background fitting
#'
#' @importFrom Matrix spMatrix
#' @importFrom methods as
#' @param x raman spectrum
#' @param w binary masks (value of the mask is zero if a point belongs to peaks
#' and one otherwise)
#' @param lambda lambda is an adjustable parameter, it can be adjusted by user.
#' The larger lambda is, the smoother z will be
#' @param differences an integer indicating the order of the difference of
#' penalties
#' @return the fitted vector
#' @author Yizeng Liang ,Zhang Zhimin
#' @seealso \code{\link{widthEstimationCWT}}
#' @keywords WhittakerSmooth

# w=x.w; differences=differences

WhittakerSmooth <- function(x, w, lambda, differences=1) {
  x <- as.vector(x)
  L <- length(x)
  E <- spMatrix(L, L, i = seq(1,L), j=seq(1, L), rep(1,L))
  D <- as(diff(E, 1, differences), "dgCMatrix")
  W <- as(spMatrix(L, L, i = seq(1, L), j = seq(1,L), w), "CsparseMatrix")
  background <- solve((W+lambda*t(D)%*%D), w*x);
  return(as.vector(background))
 }
