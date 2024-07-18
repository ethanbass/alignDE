#' Background fitting via penalized least squares algorithm
#'
#' Intelligent background-correction algorithm is developed, which simulates
#' manual background-correction procedure of an expert. It
#' basically consists of works of three aspects: 1) the accurate peak position
#' detection of Raman spectrum by continuous wavelet transform (CWT) with the
#' Mexican Hat wavelet as the mother wavelet; 2) peak widths estimation by
#' enhanced signal-to-noise ratio (SNR) derivative calculation based on CWT but
#' with the Haar wavelet as the mother wavelet; and 3) background fitting using
#' penalized least squares with binary masks.
#'
#'
#' @param x Raman spectrum
#' @param peakWidth returned by \code{\link{widthEstimationCWT}}
#' @param threshold the user define peak shape threshold
#' @param lambda lambda is an adjustable parameter, it can be adjusted by user.
#' The larger lambda is, the smoother fitted background will be
#' @param differences an integer indicating the order of the difference of
#' penalties of Whittaker Smoother method, see \code{\link{WhittakerSmooth}}
#' @return the final background
#' @author Yizeng Liang ,Zhang Zhimin
#' @seealso \code{\link{WhittakerSmooth}}
#' @keywords widthEstimationCWT
#' @examples
#'
#' data(raman)
#' x=m[7,]
#' scales <-seq(1, 70, 1)
#' wCoefs <- cwt(x, scales=scales, wavelet='mexh')
#' image(1:nrow(wCoefs), scales, wCoefs, col=terrain.colors(256),
#' axes=FALSE, xlab='index', ylab='CWT coefficient scale', main='CWT coefficients')
#' box()
#' localMax <- getLocalMaximumCWT(wCoefs)
#' ridgeList <- getRidge(localMax, gapTh=3, skip=2)
#' plotRidgeList(ridgeList)
#' majorPeakInfo = identifyMajorPeaks(x, ridgeList, wCoefs, SNR.Th=1,ridgeLength=5)
#' peakWidth=widthEstimationCWT(x,majorPeakInfo)
#' backgr = baselineCorrectionCWT(x,peakWidth,lambda=1000,differences=1)
#' corrected=x-backgr
#' plot(xa,x,type='l',ylim=c(min(c(x,corrected)),max(c(x,corrected))),
#' main="The background-correction result of Raman Spectra",
#' xlab=expression("Wavenumber / cm"^-1),ylab="Raman Intensity/Arbitr. Units")
#' points(xa[majorPeakInfo$peakIndex],x[majorPeakInfo$peakIndex])
#' lines(xa,backgr,lty=5)
#' lines(xa,corrected)
#'
# x=p1; peakWidth=peakWidth1; lambda=100; differences=1
baselineCorrectionCWT <- function(x, peakWidth, threshold=0.5,
                                  lambda=100, differences=1) {

  # fitting an rough background with proper value of lambda
  x.w <- rep(1,length(x))
  peakIndex <- peakWidth$peakIndex
  for (i in 1:length(peakIndex)){
    x.w[peakWidth[[paste(peakIndex[i])]]]=0
  }
  backgr <- WhittakerSmooth(x, x.w, lambda, differences)

  #assigning the corresponding value of the spectrum to the part of the
  #rough background,which is larger than the spectrum. and fitting another
  #background without value lager than orignal spectrum based on the new rough background

  backgr.final <- backgr
  backgr.final[x<=backgr] <- x[x<=backgr]
  backgr.final <- WhittakerSmooth(backgr.final,rep(1,length(backgr.final)),1,differences)
  x <- backgr.final

  #refine the new rough background

  background <- NULL
  peakIndex <- peakWidth$peakIndex
  signal_baseline <- x[1:(peakWidth[[paste(peakIndex[1])]][1])-1]
  signal_baseline.w <- rep(1,length(signal_baseline))
  baseline_baseline <-  WhittakerSmooth(signal_baseline,signal_baseline.w,lambda,differences)
  background <- c(background,baseline_baseline)
  for (i in 1:(length(peakIndex)-1)){

    peakWidth.1 <- peakWidth[[paste(peakIndex[i])]]
    peakWidth.2 <- peakWidth[[paste(peakIndex[i+1])]]

    if (length(intersect(peakWidth.1, peakWidth.2))==0){
      signal_peak=x[peakWidth.1]
      signal_peak.w <- c(1,rep(0,(length(peakWidth.1)-2)),1)
      if((abs(signal_peak[length(signal_peak)] - signal_peak[1])/mean(signal_peak-min(signal_peak))) >= threshold){
        if(signal_peak[length(signal_peak)] > signal_peak[1]){
          signal_peak.w=as.numeric(!(signal_peak>signal_peak[length(signal_peak)]))
        } else{
          signal_peak.w=as.numeric(!(signal_peak>signal_peak[1]))
        }
        signal_peak.w[1] <- 1
        signal_peak.w[length(signal_peak.w)] <- 1
      }
      peak_baseline <- WhittakerSmooth(signal_peak,signal_peak.w,1,differences)
      background  <- c(background,peak_baseline)

      if(peakWidth.2[1] - peakWidth.1[length(peakWidth.1)]==1){
        # two peak just together
      } else{
        signal_baseline <- x[(peakWidth.1[length(peakWidth.1)]+1):(peakWidth.2[1]-1)]
        if(length(signal_baseline) <= 3){
          baseline_baseline <- signal_baseline
        }else{
          signal_baseline.w <- c(1,rep(1,(length(signal_baseline)-2)),1)
	        baseline_baseline <- WhittakerSmooth(signal_baseline, 
	                                             signal_baseline.w,
	                                             lambda, differences)
        }
        background <- c(background,baseline_baseline)
      }

    } else{
      if (peakWidth.1[1]>peakWidth.2[1]){
         # the last peak contain the previous peak
         peakWidth[[paste(peakIndex[i+1])]] <- 
           (peakWidth.1[1]):(peakWidth.2[length(peakWidth.2)])
      } else{
          if (length(peakWidth.1) >= length(peakWidth.2)){
            peakWidth.11 <- setdiff(peakWidth.1,intersect(peakWidth.1, peakWidth.2))
        } else{
          peakWidth.11 <- peakWidth.1
          peakWidth[[paste(peakIndex[i+1])]] <- 
            setdiff(peakWidth.2,intersect(peakWidth.1, peakWidth.2))
        }

        if(length(peakWidth.11) <= 2){
          peakWidth[[paste(peakIndex[i+1])]] <- 
            (peakWidth.1[1]):(peakWidth.2[length(peakWidth.2)])
        } else{
          signal_peak <- x[peakWidth.11]
          signal_peak.w <- c(1,rep(0,(length(peakWidth.11)-2)),1)
          if((abs(signal_peak[length(signal_peak)]-signal_peak[1])/mean(signal_peak-min(signal_peak)))>=threshold){
            if (signal_peak[length(signal_peak)] > signal_peak[1]){
              signal_peak.w <- as.numeric(!(signal_peak > signal_peak[length(signal_peak)]))
            } else{
              signal_peak.w <- as.numeric(!(signal_peak>signal_peak[1]))
            }
            signal_peak.w[1] <- 1
            signal_peak.w[length(signal_peak.w)] <- 1
          }
          peak_baseline <- WhittakerSmooth(signal_peak,signal_peak.w,1,differences)
          background  <- c(background,peak_baseline)
        }
      }
    }
  }

  peakWidth.end <- peakWidth[[paste(peakIndex[length(peakIndex)])]]
  signal_peak <- x[peakWidth.end]
  signal_peak.w <- c(1,rep(0,(length(peakWidth.end)-2)),1)
  if((abs(signal_peak[length(signal_peak)]-signal_peak[1])/mean(signal_peak-min(signal_peak)))>=threshold){
    if(signal_peak[length(signal_peak)] > signal_peak[1]){
      signal_peak.w <- as.numeric(!(signal_peak>signal_peak[length(signal_peak)]))
    }else{
      signal_peak.w <- as.numeric(!(signal_peak>signal_peak[1]))
    }
  }
  peak_baseline <- WhittakerSmooth(signal_peak,signal_peak.w,1,differences)

  signal_baseline <- x[(peakWidth.end[length(peakWidth.end)]+1):length(x)]
  signal_baseline.w <- rep(1,length(signal_baseline))
  baseline_baseline  <- WhittakerSmooth(signal_baseline,signal_baseline.w,lambda,differences)
  background  <- c(background,peak_baseline)
  background  <- c(background,baseline_baseline)

  return(background)
}
