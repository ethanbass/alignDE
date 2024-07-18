#' Plot the estimation of peak width based on the CWT
#' 
#' Plot the the estimation of peak returned by \code{\link{widthEstimationCWT}}
#' 
#' 
#' @param x raman spectrum
#' @param peakWidth returned by \code{\link{widthEstimationCWT}}
#' @return Whether the drawing is successful or not.
#' @author Yizeng Liang ,Zhang Zhimin
#' @seealso \code{\link{widthEstimationCWT}}
#' @keywords plotwidthEstimation
#' @examples
#' 
#' x=c(40*dnorm(seq(-5,5,by=0.1),sd=0.2),60*dnorm(seq(-5,5,by=0.1),sd=0.5),30*dnorm(seq(-5,5,by=0.1),sd=0.7),100*rep(0.001,100))
#' x=x+40*rnorm(length(x))*0.01
#' b=1:length(x)
#' yc=30*sin(1:length(x)/100)
#' yl=30*(1:length(x)/100)
#' xl=x+yl
#' xc=x+yc
#' x=xc  
#' scales <-seq(1, 30, 1)
#' wCoefs <- cwt(x, scales=scales, wavelet='mexh')
#' image(1:nrow(wCoefs), scales, wCoefs, col=terrain.colors(256), axes=FALSE, xlab='index', ylab='CWT coefficient scale', main='CWT coefficients')
#' box()
#' localMax <- getLocalMaximumCWT(wCoefs)
#' ridgeList <- getRidge(localMax, gapTh=3, skip=2)
#' majorPeakInfo=identifyMajorPeaks(x, ridgeList, wCoefs,SNR.Th=3,ridgeLength=5)
#' peakWidth=widthEstimationCWT(x,majorPeakInfo)  
#' plot(xc,type='l', main="pure signal with curved background and random noise",xlab=expression("Wavenumber / cm"^-1),ylab="Raman Intensity")
#' plotwidthEstimation(x,peakWidth)   
#' 
"plotwidthEstimation" <-
function(x,peakWidth) {
  x=as.vector(x)
  lmindex=1:length(x)
  peakIndex=peakWidth$peakIndex
  LR = NULL
#  plot(x,type='l')                       
  points(peakIndex,x[peakIndex])
  for (i in 1:length(peakIndex)){
    peakWidth.i=peakWidth$peakIndexLower[i]:peakWidth$peakIndexUpper[i]
    LR=c(LR,peakWidth.i[c(1,length(peakWidth.i))])  
  }
  points(LR,x[LR])
  
  return ("successful")
}
