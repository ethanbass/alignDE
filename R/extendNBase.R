#' Extend the row number of a matrix as the exponential of base N
#'
#' Extend the data as the exponential of base N by increasing row number.
#'
#' The method 'open' is padding the the matrix with the last row.
#'
#' @importFrom stats nextn
#' @param x data matrix
#' @param nLevel the level of DWT decomposition. Basically, it is equivalent to
#' changing the 'base' as base\^nLevel
#' @param base the base, 2 by default
#' @param \dots other parameters of used by \code{\link{extendLength}}
#' @return Return a extended matrix
#' @author Pan Du
#' @seealso \code{\link{extendLength}}
#' @keywords methods
#' @examples
#'
#' 	# a = matrix(rnorm(9), 3)
#' 	# extendNBase(a)	## not exposed function
#'
extendNBase <- function(x, nLevel=1, base=2, ...){
	if (!is.matrix(x)) x <- matrix(x, ncol=1)

	nR <- nrow(x)
	if (is.null(nLevel)) {
		nR1 <- nextn(nR, base)
	} else {
		nR1 <- ceiling(nR / base^nLevel) * base^nLevel
	}
	if (nR != nR1) {
		x <- extendLength(x, addLength = nR1-nR, ...)
	}

	return(x)
}

