#' Estimation of peak width based on the CWT
#'
#' Peak widths estimation by enhanced signal-to-noise ratio (SNR) derivative
#' calculation based on CWT but with the Haar wavelet as the mother wavelet
#'
#'
#' @param x Raman spectrum
#' @param majorPeakInfo returned by \code{\link{identifyMajorPeaks}}
#' @return A list cotains the peak position and the indexs of every peak.
#' @author Yizeng Liang ,Zhang Zhimin
#' @seealso \code{\link{identifyMajorPeaks}}
#' @keywords widthEstimationCWT
#' @examples
#'
#' x <- c(40*dnorm(seq(-5,5,by=0.1),sd=0.2),60*dnorm(seq(-5,5,by=0.1),sd=0.5),
#' 30*dnorm(seq(-5,5,by=0.1),sd=0.7),100*rep(0.001,100))
#' x <- x+40*rnorm(length(x))*0.01
#' b <- 1:length(x)
#' yc <- 30*sin(1:length(x)/100)
#' yl <- 30*(1:length(x)/100)
#' xl <- x+yl
#' xc <- x+yc
#' x <- xc
#' scales <- seq(1, 30, 1)
#' wCoefs <- cwt(x, scales=scales, wavelet='mexh')
#' image(1:nrow(wCoefs), scales, wCoefs, col=terrain.colors(256), axes=FALSE, xlab='index', ylab='CWT coefficient scale', main='CWT coefficients')
#' box()
#' localMax <- getLocalMaximumCWT(wCoefs)
#' ridgeList <- getRidge(localMax, gapTh=3, skip=2)
#' majorPeakInfo <- identifyMajorPeaks(x, ridgeList, wCoefs,SNR.Th=3,ridgeLength=5)
#' peakWidth <- widthEstimationCWT(x,majorPeakInfo)
#' plot(xc,type='l', main="pure signal with curved background and random noise",xlab=expression("Wavenumber / cm"^-1),ylab="Intensity")
#' plotwidthEstimation(x,peakWidth)
#'
widthEstimationCWT <- function(x, majorPeakInfo) {

  wCoefs_haar <- cwt(x, 1:max(majorPeakInfo$peakScale), wavelet='haar')
  peakIndex <- majorPeakInfo$peakIndex
  peakScale <- majorPeakInfo$peakScale[findInterval(majorPeakInfo$peakIndex,majorPeakInfo$allPeakIndex)]

  peakWidth <- list()
  peakWidth[["peakIndex"]] <- majorPeakInfo$peakIndex
  peakWidth[["peakIndexLower"]] <- majorPeakInfo$peakIndex
  peakWidth[["peakIndexUpper"]] <- majorPeakInfo$peakIndex

  for(i in 1:length(peakIndex)){
    peakIndex.i <- peakIndex[i]
    peakScale.i <- peakScale[i]
    wCoefs_haar.i <- wCoefs_haar[,peakScale.i]
    wCoefs_haar.i.abs <- abs(wCoefs_haar.i)
    localmax <- localMaximum(-wCoefs_haar.i.abs,winSize=5)
    #localmax=localmax & abs(wCoefs_haar.i)<(mean(wCoefs_haar.i.abs[localmax==1])+0.5*sd(wCoefs_haar.i.abs[localmax==1]))
    localmax <- as.numeric(localmax)
    localmax[peakIndex] <- 0
    localmax[(peakIndex.i-peakScale.i/2+1):(peakIndex.i+peakScale.i/2-1)] <- 0

    Lef <- 0
    Rig <- 0

    peakScale.i.3 <- 2*peakScale.i


    if (i == 1){
       maxIndexL <- max(c((peakIndex.i - peakScale.i.3), 1))
    } else{
       maxIndexL <- max(c((peakIndex.i - peakScale.i.3), peakIndex[i-1]))
    }

    if (i == length(peakIndex)){
       minIndexR <- min(c((peakIndex.i+peakScale.i.3), length(localmax)))
    } else{
       minIndexR <- min(c((peakIndex.i+peakScale.i.3),peakIndex[i+1]))
    }
    ignoreL <- 1:maxIndexL
    ignoreR <- minIndexR:length(localmax)
    localmax[c(ignoreL,ignoreR)] <- 0
    localmax[c(peakIndex.i, (peakIndex.i-(peakScale.i/2)):(peakIndex.i+(peakScale.i/2)))] <- 0
    bi <- which(localmax==1)

    biLeft <- bi[bi < peakIndex.i]
    useL <- maxIndexL:peakIndex.i
    minIndexLeft <- useL[tail(which(min(x[useL]) == x[useL]), 1)]

    if (length(biLeft) == 0){
      Lef <- minIndexLeft
    } else{
      minbaselineIndexLeft <- biLeft[which(min(x[biLeft]) == x[biLeft])]
      if (minIndexLeft >= (peakIndex.i-peakScale.i/2+1)){
        Lef <- minbaselineIndexLeft
      } else{
        Lef <- max(c(minIndexLeft,minbaselineIndexLeft))
      }
    }

    biRight <- bi[bi>peakIndex.i]
    useR <- peakIndex.i:minIndexR
    minIndexRight <- useR[which(min(x[useR]) == x[useR])[1]]

    if (length(biRight) == 0){
      Rig <- minIndexRight
    } else{
      minbaselineIndexRight <- biRight[which(min(x[biRight]) == x[biRight])]
      if (minIndexRight <= (peakIndex.i+peakScale.i/2 - 1)){
        Rig <- minbaselineIndexRight
      } else{
        Rig <- min(c(minIndexRight,minbaselineIndexRight))
      }
    }
    peakWidth[[paste(peakIndex.i)]] <- Lef:Rig
    peakWidth[["peakIndexLower"]][i] <- Lef
    peakWidth[["peakIndexUpper"]][i] <- Rig
  }
  return(peakWidth)
}

