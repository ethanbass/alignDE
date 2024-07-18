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
#' @importFrom stats approx
#' @param x chromatogram to be aligned
#' @param peakWidth returned by \code{\link{widthEstimationCWT}}
#' @param target reference chromatograms
#' @param slack shifts the peak position can be adjusted
#' @param n Gap for peak clustering
#' @param control List datatype for DEoptim
#' @param verbose Whether to print progress to console.
#' @return the aligned chromatograms
#' @author Yizeng Liang ,Zhang Zhimin, Chen Shan
#' @keywords alignDE
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
#' majorPeakInfo <- identifyMajorPeaks(p2, ridgeList, wCoefs, SNR.Th=3,ridgeLength=5)
#' peakWidth2 <- widthEstimationCWT(p2,majorPeakInfo)
#' backgr <- baselineCorrectionCWT(p2,peakWidth2,lambda=100,differences=1)
#' p2c <- p2-backgr
#' #plot(p2c,type='l')
#' #points(majorPeakInfo$peakIndex,p2c[majorPeakInfo$peakIndex])
#'
#' plot(p1c,type='l',main='Peak clustered and baseline corrected chromatograms',xlab=xlab,ylab=ylab)
#' lines(p2c,lty=3)
#' peakWidth2 <- peakClustering(peakWidth2,n=5)
#' plotwidthEstimation(p2c,peakWidth2)
#' legend(780, 120, c("R", "C"),
#'        text.col = "black", lty = c(1, 3))
#'
#' ########################################################
#'
#' cp <- list(NP=50, itermax =150, refresh = 10)
#' result <- alignDE(p2c,peakWidth2,p1c,100,control=cp)
#'
#' #######################################################
#'
#' plot(p1c,type='l',main='Aligned and baseline corrected chromatograms',xlab=xlab,ylab=ylab)
#' #lines(p2c,col='red')
#' lines(result,lty=3)
#' legend(780, 120, c("R", "A"),
#'        text.col = "black", lty = c(1, 3))
#' @export

alignDE <- function(x, peakWidth, target, slack = 70, n = 4,
                    verbose = getOption("verbose"),
                    control = list(termax = 100, CR = 0.7, refresh = 1,
                                 verbose = getOption("verbose"))){

#  control=list(NP=50, itermax =100,CR = 0.7, refresh = 10)
#  slack=100
#  x=p2c
#  target=p1c
#  peakWidth=peakWidth21

  nPeakNum <- length(peakWidth$peakIndex)
  result <- list()
  shift <- c()
  for(i in 1:floor(nPeakNum/n))
  {
    lower=rep(-slack, n)
    upper=rep(slack, n)
    peakWidths <- list()
    peakWidths$peakIndex <- peakWidth$peakIndex[(1+((i-1)*n)):(n+(i-1)*n)]
    peakWidths$peakIndexLower <- peakWidth$peakIndexLower[(1+((i-1)*n)):(n+(i-1)*n)]
    peakWidths$peakIndexUpper <- peakWidth$peakIndexUpper[(1+((i-1)*n)):(n+(i-1)*n)]

    if (i == 1){
      startindex <- 1
    } else{
      startindex <- peakWidth$peakIndexUpper[(1+((i-1)*n))-1]+shift[length(shift)]+1
    }
    endindex <- length(x)
    result[[i]] <- DEoptim(objfun,lower = lower,upper = upper,
                           control = control, peakWidth=peakWidths,
                        x=x, target=target, startindex=startindex,
                        endindex=endindex, verbose=verbose)
    shift <- c(shift,result[[i]]$optim$bestmem)
  }


  nLeftPeaks <- nPeakNum%%n
  i <- i+1
  if (nLeftPeaks >= 1){
    lower <- rep(-slack,nLeftPeaks)
    upper <- rep(slack,nLeftPeaks)
    peakWidths <- list()
    peakWidths$peakIndex <- peakWidth$peakIndex[(nPeakNum-nLeftPeaks+1):nPeakNum]
    peakWidths$peakIndexLower <- peakWidth$peakIndexLower[(nPeakNum-nLeftPeaks+1):nPeakNum]
    peakWidths$peakIndexUpper <- peakWidth$peakIndexUpper[(nPeakNum-nLeftPeaks+1):nPeakNum]
    if (i == 1){
      startindex <- 1
    } else{
      startindex <- peakWidth$peakIndexUpper[(1 + ((i-1)*n)) - 1] + 
        shift[length(shift)]+1
    }
    endindex <- length(x)
    result[[i+1]] <- DEoptim(objfun, lower, upper, control, peakWidth=peakWidths,
                          x=x,target=target,startindex=startindex,
                          endindex=endindex, verbose = verbose)
    shift <- c(shift,result[[i+1]]$optim$bestmem)
  }
#######################################################################################

  x <- as.vector(x)
  shift <- ceiling(shift)
  peakIndex <- peakWidth$peakIndex
  peak <- peakWidth$peakIndexLower[1]:peakWidth$peakIndexUpper[1]
  xs=c()
    baseline_warp <- approx(x[1:(peak[1]-1)],n=peak[1]+shift[1]-1)
  xs <- c(xs,baseline_warp$y)
  for(i in 1:(nPeakNum-1)){
    peakCurrent <- peakWidth$peakIndexLower[i]:peakWidth$peakIndexUpper[i]
    peakNext <- peakWidth$peakIndexLower[i+1]:peakWidth$peakIndexUpper[i+1]
    baselineStart <- peakCurrent[length(peakCurrent)]+1
    baselineEnd <- peakNext[1]-1

    xs <- c(xs,x[peakCurrent])
    baseline_warp <- approx(x[baselineStart:baselineEnd],
                         n=((baselineEnd+shift[i+1])-(baselineStart+shift[i])+1) )
    xs <- c(xs,baseline_warp$y)
  }
  peak <- peakWidth$peakIndexLower[nPeakNum]:peakWidth$peakIndexUpper[nPeakNum]
  xs <- c(xs, x[peak])
  baseline_warp <- approx(x[(peak[length(peak)]+1):length(x)],
                          n=(length(x)-(peak[length(peak)]+shift[nPeakNum]+1)+1))
  xs <- c(xs,baseline_warp$y)
###############################################################################################

  return(xs)
}
