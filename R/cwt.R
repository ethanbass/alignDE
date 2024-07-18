#' Continuous Wavelet Transform (CWT)
#'
#' CWT(Continuous Wavelet Transform) with Mexican Hat wavelet (by default) to
#' match the peaks in Mass Spectrometry spectrum
#'
#' @importFrom stats convolve
#' @param ms Mass Spectrometry spectrum (a vector of MS intensities)
#' @param scales a vector represents the scales at which to perform CWT.
#' @param wavelet The wavelet base, Mexican Hat by default. User can provide
#' wavelet Psi(x) as a form of two row matrix. The first row is the x value,
#' and the second row is Psi(x) corresponding to x.
#' @return The return is the 2-D CWT coefficient matrix, with column names as
#' the scale. Each column is the CWT coefficients at that scale.
#' @author Pan Du, Simon Lin
#' @keywords methods
#' @examples
#'
#' 	data(exampleMS)
#' 	scales <- seq(1, 64, 3)
#' 	wCoefs <- cwt(exampleMS[5000:11000], scales=scales, wavelet='mexh')
#'
#' 	## Plot the 2-D CWT coefficients as image (It may take a while!)
#' 	xTickInterval <- 1000
#' 	image(5000:11000, scales, wCoefs, col=terrain.colors(256), axes=FALSE, xlab='m/z index', ylab='CWT coefficient scale', main='CWT coefficients')
#' 	axis(1, at=seq(5000, 11000, by=xTickInterval))
#' 	axis(2, at=c(1, seq(10, 64, by=10)))
#' 	box()
#'
cwt <- function(ms, scales = 1, wavelet = 'mexh') {
	## Check for the wavelet format
	if (wavelet == 'mexh') {
		psi_xval <- seq(-8, 8, length=1024)
		psi <- (2/sqrt(3) * pi^(-0.25)) * (1 - psi_xval^2) * exp(-psi_xval^2/2)
		#plot(psi_xval, psi)
	} else if (wavelet=='haar') {
    psi_xval <- seq(0,1,length=1024)
    psi <- c(0,rep(1,511),rep(-1,511),0)
	} else if (is.matrix(wavelet)) {
		if (nrow(wavelet) == 2) {
			psi_xval <- wavelet[1,]
			psi <- wavelet[2,]
		} else if (ncol(wavelet) == 2) {
			psi_xval <- wavelet[,1]
			psi <- wavelet[,2]
		} else {
			stop('Unsupported wavelet format!')
		}
	} else {
		  stop('Unsupported wavelet!')
	  }
    oldLen <- length(ms)
	## To increase the computation effeciency of FFT, extend it as the power of 2
	## because of a strange signal length 21577 makes the FFT very slow!
	#ms <- extend.nBase(ms, nLevel=1, base=2)
	ms <- extendNBase(ms, nLevel=NULL, base=2)
	len <- length(ms)
    nbscales <- length(scales)
    wCoefs <- NULL

    psi_xval <- psi_xval - psi_xval[1]
    dxval <- psi_xval[2]
    xmax  <- psi_xval[length(psi_xval)]
    for (i in 1:length(scales)) {
		scale.i <- scales[i]
		f <- rep(0, len)
        j <- 1 + floor((0:(scale.i * xmax))/(scale.i * dxval))
        if (length(j) == 1)		j <- c(1, 1)
		lenWave <- length(j)
        f[1:lenWave] <- rev(psi[j]) - mean(psi[j])
		if (length(f) > len) stop(paste('scale', scale.i, 'is too large!'))
		wCoefs.i <- 1/sqrt(scale.i) * convolve(ms, f)
		## Shift the position with half wavelet width
		wCoefs.i <- c(wCoefs.i[(len-floor(lenWave/2) + 1) : len], wCoefs.i[1:(len-floor(lenWave/2))])
		wCoefs <- cbind(wCoefs, wCoefs.i)
    }
	if (length(scales) == 1) wCoefs <- matrix(wCoefs, ncol=1)
	colnames(wCoefs) <- scales
	wCoefs <- wCoefs[1:oldLen,,drop=FALSE]
	return(wCoefs)
}

