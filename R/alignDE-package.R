

#' Peak alignment using wavelet pattern matching and differential evolution
#' 
#' Retention time shifts badly impair qualitative or quantitative results of
#' chemometric analyses when entire chromatographic data are used. Hence,
#' chromatograms should be aligned to perform further analysis. Being inspired
#' and motivated by this purpose, a practical and handy peak alignment method
#' (alignDE) is proposed, implemented in this research for first-order
#' chromatograms, which basically consists of five steps: (1) chromatograms
#' lengths equalization using linear interpolation; (2) accurate peak pattern
#' matching by continuous wavelet transform (CWT) with the Mexican Hat and Haar
#' wavelets as its mother wavelets; (3) flexible baseline fitting utilizing
#' penalized least squares; (4) peak clustering when gap of two peaks is larger
#' than a certain threshold; (5) peak alignment using differential evolution
#' (DE) to maximize linear correlation coefficient between reference signal and
#' signal to be aligned. This method is demonstrated with simulated
#' chromatograms, chromatograms of HPLC-DAD at 202nm. It is implemented in R
#' language and available as open source software to a broad range of
#' chromatograph users (http://code.google.com/p/alignde).
#' 
#' \tabular{ll}{ Package: \tab alignDE\cr Type: \tab Package\cr Version: \tab
#' 1.0.0\cr Date: \tab 2009-10-09\cr License: \tab GPL (>= 2)\cr }
#' 
#' @aliases alignDE-package alignDE
#' @docType alignDE
#' @author yizeng liang<yizeng_liang@@263.net>, zhimin zhang
#' <zhangzhimin.csu@@gmail.com>, chen shan <chenshan.csu@@gmail.com>
#' @keywords package
NULL





#' DEoptim-methods
#' 
#' Methods for DEoptim objects.
#' 
#' 
#' @aliases DEoptim-methods plot.DEoptim summary.DEoptim
#' @param object,x An object of class~\code{DEoptim}; usually, a result of a
#' call to~\code{\link{DEoptim}}.
#' @param plot.type Should we plot the best member at each iteration, the best
#' value at each iteration or the intermediate populations?
#' @param \dots Further arguments passed to or from other methods.
#' @note Please cite the package in publications. Use
#' \code{citation("DEoptim")}.
#' @author David Ardia~\email{david.ardia@@unifr.ch}
#' @keywords methods
#' @examples
#' 
#'   ## Rosenbrock Banana function
#'   Rosenbrock <- function(x){
#'     x1 <- x[1]
#'     x2 <- x[2]
#'     100 * (x2 - x1 * x1)^2 + (1 - x1)^2
#'   }
#' 
#'   lower <- c(-10,-10)
#'   upper <- -lower
#'   r <- DEoptim(Rosenbrock, lower, upper)
#'   summary(r) ## print summary of the DEoptim object
#' 
#'   par(mfrow = c(1,2))
#'   plot(r, type = 'b') ## plot the best members
#'   plot(r, plot.type = "bestvalit", type = 'b', col = 'blue') ## plot the best values
#' 
#'   ## rerun the optimization, and store intermediate populations
#'   r <- DEoptim(Rosenbrock, lower, upper, control=list(NP=400, itermax=30, storepopfrom=1, storepopfreq=2))
#'   summary(r)
#'   
#'   par(mfrow=c(1,1))
#'   plot(r, plot.type = "storepop") ## plot intermediate populations
#' 
NULL





#' Estimate the length of the ridge
#' 
#' Estimate the length of the ridge line, which is composed of local maxima at
#' adjacent CWT scales. The ridge line is cut off at the end point, whose
#' amplitude divided by the maximum ridge amplitude is larger than the cutoff
#' amplitude ratio threshold (0.5 by default).
#' 
#' 
#' @param ridgeList a list of identified ridges
#' @param Th the cutoff amplitude ratio (the amplitude divided by the maximum
#' amplitude of the ridge) threshold of the ridge line end.
#' @return a vector of estimated ridge length
#' @author Pan Du
#' @keywords methods
NULL





#' Get the CWT coefficient values corresponding to the peak ridge
#' 
#' Get the CWT coefficient values corresponding to the peak ridge
#' 
#' 
#' @param ridgeList a list of ridge lines
#' @param wCoefs 2-D CWT coefficients
#' @param skip the CWT scale level to be skipped, by default the 0 scale level
#' (raw spectrum) is skipped.
#' @return A list of ridge values corresponding to the input ridgeList.
#' @author Pan Du
#' @keywords methods
NULL





#' Real chromatograms
#' 
#' Real chromatograms are available from Refs. [11,35]. The chromatograms were
#' analyses of extracts from fungal cultivated on Yeast Extract Sucrose agar
#' (P. cyclopium, denoted by IBT 11415 and 15670) using ultrasonic extraction
#' and HPLC, which collected at the Department of Biotechnology (IBT),
#' Technical University of Denmark. IBT11415 and IBT15670 were downloaded as an
#' MATLABTM 6 MAT-file from website mentioned in Ref. [11]
#' 
#' The eal chromatograms can be used for aligning.
#' 
#' @docType IBT
#' @format two vector \describe{ \item{list("IBT11415 ")}{vector 1 as
#' reference} \item{list("IBT15670 ")}{vector 2 to be aligned} }
#' @keywords datasets
#' @examples
#'  
#' require(alignDE)
#' data(IBT)
#' main="Chromatograms of IBT"
#' xlab ="Sample intervals"
#' ylab="mAU"
#' p1=IBT11415
#' p2=IBT15670
#' plot(p1,type='l',main=main,xlab=xlab,ylab=ylab)
#' lines(p2,lty=3)
#' legend(2400, 750, c("IBT11415", "IBT15670"),
#'        text.col = "black", lty = c(1, 3))
#' #lines(warped,col='blue')
#' 
#' 
#' scales <-seq(1, 56, 1)
#' wCoefs <- cwt(p1, scales=scales, wavelet='mexh')
#' image(1:nrow(wCoefs), scales,apply(t(wCoefs),1,rev), col=terrain.colors(256), axes=FALSE, xlab='index', ylab='CWT coefficient scale', main='CWT coefficients')
#' box()
#' localMax <- getLocalMaximumCWT(wCoefs)
#' ridgeList <- getRidge(localMax, gapTh=3, skip=2)
#' plotRidgeList(ridgeList)
#' majorPeakInfo = identifyMajorPeaks(p1, ridgeList, wCoefs, SNR.Th=3,ridgeLength=5)
#' peakWidth1=widthEstimationCWT(p1,majorPeakInfo)
#' backgr = baselineCorrectionCWT(p1,peakWidth1,lambda=100,differences=1)
#' p1c=p1-backgr
#' #plot(p1c,type='l')
#' #points(majorPeakInfo$peakIndex,p1c[majorPeakInfo$peakIndex])
#' 
#' scales <-seq(1, 56, 1)
#' wCoefs <- cwt(p2, scales=scales, wavelet='mexh')
#' image(1:nrow(wCoefs), scales,apply(t(wCoefs),1,rev), col=terrain.colors(256), axes=FALSE, xlab='index', ylab='CWT coefficient scale', main='CWT coefficients')
#' box()
#' localMax <- getLocalMaximumCWT(wCoefs)
#' ridgeList <- getRidge(localMax, gapTh=3, skip=2)
#' plotRidgeList(ridgeList)
#' majorPeakInfo = identifyMajorPeaks(p2, ridgeList, wCoefs, SNR.Th=3,ridgeLength=5)
#' peakWidth2=widthEstimationCWT(p2,majorPeakInfo)
#' backgr = baselineCorrectionCWT(p2,peakWidth2,lambda=100,differences=1)
#' p2c=p2-backgr
#' #plot(p2c,type='l')
#' #points(majorPeakInfo$peakIndex,p2c[majorPeakInfo$peakIndex])
#' x11()
#' plot(p1c,type='l',main='Peak clustered and baseline corrected chromatograms',xlab=xlab,ylab=ylab)
#' lines(p2c,lty=3)
#' points(majorPeakInfo$peakIndex,p2c[majorPeakInfo$peakIndex])
#' peakWidth2=peakClustering(peakWidth2,n=5)
#' #plotwidthEstimation(p2c,peakWidth2)
#' legend(2400, 750, c("IBT11415", "IBT15670"),
#'        text.col = "black", lty = c(1, 3))
#' 
#' ########################################################
#' 
#' cp=list(NP=200, itermax = 150, refresh = 10)
#' result=alignDE(p2c,peakWidth2,p1c,slack=100,control=cp)
#' 
#' #######################################################
#' similarity(result,p1c)
#' 
#' plot(p1c,type='l',main='Aligned and baseline corrected chromatograms',xlab=xlab,ylab=ylab)
#' #lines(p2c,col='red')
#' lines(result,lty=3)
#' legend(2400, 750, c("IBT11415", "IBT15670"),
#'        text.col = "black", lty = c(1, 3))
#' 
NULL





#' Simulated chromatograms
#' 
#' Simulated chromatogram is the sum of Gaussian peaks, sinus curve baseline
#' and random noise.
#' 
#' Simulated chromatogram is the sum of Gaussian peaks, sinus curve baseline
#' and random noise. Reference one is denoted as R, whose noise is normally
#' distributed, with variance 0.2. The one, to be aligned, is denoted as C, but
#' normally distributed noise with variance 1. The peaks of C were shifted by
#' 50 positions from peaks of R except the second peak. Both R and C are
#' created in R language, and shown in Fig. 6(a).
#' 
#' @name Simulate
#' @docType data
#' @format two vector \describe{ \item{list("p1 ")}{vector 1 as reference}
#' \item{list("p2 ")}{vector 2 to be aligned} }
#' @keywords datasets
#' @examples
#' 
#' require(alignDE) 
#' data(simulate)
#' main="Simulated chromatograms"
#' xlab ="Sample intervals"
#' ylab="mAU"
#' plot(p1,type='l',main=main,xlab=xlab,ylab=ylab)
#' lines(p2,lty=3)
#' legend(780, 120, c("R", "C"),
#'        text.col = "black", lty = c(1, 3))
#' #lines(warped,col='blue')
#' 
#' 
#' scales <-seq(1, 56, 1)
#' wCoefs <- cwt(p1, scales=scales, wavelet='mexh')
#' image(1:nrow(wCoefs), scales,apply(t(wCoefs),1,rev), col=terrain.colors(256), axes=FALSE, xlab='index', ylab='CWT coefficient scale', main='CWT coefficients')
#' box()
#' localMax <- getLocalMaximumCWT(wCoefs)
#' ridgeList <- getRidge(localMax, gapTh=3, skip=2)
#' plotRidgeList(ridgeList)
#' majorPeakInfo = identifyMajorPeaks(p1, ridgeList, wCoefs, SNR.Th=3,ridgeLength=5)
#' peakWidth1=widthEstimationCWT(p1,majorPeakInfo)
#' backgr = baselineCorrectionCWT(p1,peakWidth1,lambda=100,differences=1)
#' p1c=p1-backgr 
#' #plot(p1c,type='l')
#' #points(majorPeakInfo$peakIndex,p1c[majorPeakInfo$peakIndex])                                         
#' 
#' scales <-seq(1, 56, 1)
#' wCoefs <- cwt(p2, scales=scales, wavelet='mexh')
#' image(1:nrow(wCoefs), scales,apply(t(wCoefs),1,rev), col=terrain.colors(256), axes=FALSE, xlab='index', ylab='CWT coefficient scale', main='CWT coefficients')
#' box()
#' localMax <- getLocalMaximumCWT(wCoefs)
#' ridgeList <- getRidge(localMax, gapTh=3, skip=2)
#' plotRidgeList(ridgeList)
#' majorPeakInfo = identifyMajorPeaks(p2, ridgeList, wCoefs, SNR.Th=3,ridgeLength=5)
#' peakWidth2=widthEstimationCWT(p2,majorPeakInfo)
#' backgr = baselineCorrectionCWT(p2,peakWidth2,lambda=100,differences=1)
#' p2c=p2-backgr 
#' #plot(p2c,type='l')
#' #points(majorPeakInfo$peakIndex,p2c[majorPeakInfo$peakIndex])       
#' 
#' plot(p1c,type='l',main='Peak clustered and baseline corrected chromatograms',xlab=xlab,ylab=ylab)
#' lines(p2c,lty=3)
#' peakWidth2=peakClustering(peakWidth2,n=5)
#' plotwidthEstimation(p2c,peakWidth2)
#' legend(780, 120, c("R", "C"),
#'        text.col = "black", lty = c(1, 3))
#' 
#' ########################################################
#' 
#' cp=list(NP=50, itermax =150, refresh = 10)
#' result=alignDE(p2c,peakWidth2,p1c,100,control=cp)
#' 
#' #######################################################
#' 
#' plot(p1c,type='l',main='Aligned and baseline corrected chromatograms',xlab=xlab,ylab=ylab)
#' #lines(p2c,col='red')
#' lines(result,lty=3)
#' legend(780, 120, c("R", "A"),
#'        text.col = "black", lty = c(1, 3))
#' 
#' 
NULL



