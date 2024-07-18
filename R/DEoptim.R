##'.onLoad' <- function(lib, pkg){
##  cat("\nDEoptim package")
##  cat("\nDifferential Evolution algorithm")
##  cat("\nAuthor and maintainer : David Ardia <david.ardia@unifr.ch>\n")
##}

## Differential Evolution Optimization
## David Ardia -- 20081203

#' Differential Evolution Optimization
#'
#' Performs evolutionary optimization via the Differential Evolution algorithm.
#'
#' \code{DEoptim} performs minimization of \code{FUN}.
#'
#' The \code{control} argument is a list that can supply any of the following
#' components:
#'
#' \describe{ \item{list("VTR")}{The value to be reached. The optimization
#' process will stop if either the maximum number of iterations \code{itermax}
#' is reached or the best parameter vector \code{bestmem} has found a value
#' \code{FUN(bestmem) <= VTR}. Default to \code{-Inf}.}
#' \item{list("itermax")}{The maximum iteration (population generation)
#' allowed.  Default is \code{200}.} \item{list("NP")}{Number of population
#' members. Default to \code{50}.} \item{list("F")}{Stepsize from interval
#' [0,2]. Default to \code{0.8}.} \item{list("CR")}{Crossover probability from
#' interval [0,1]. Default to \code{0.5}.} \item{list("initial")}{An initial
#' population used as a starting population in the optimization procedure.
#' Maybe useful to speed up the convergence. Defaults to \code{NULL}.}
#' \item{list("storepopfrom")}{From which population should the following
#' intermediate populations be stored in memory. Default to \code{itermax+1},
#' i.e., no intermediate population is stored.} \item{list("storepopfreq")}{The
#' frequency of populations' storage. Default to \code{1}, i.e. every
#' intermediate population is memorized.} \item{list("strategy")}{Defines the
#' binomial DE-strategy used in the optimization procedure: \describe{
#' \item{list("1")}{best/1} \item{list("2")}{rand/1}
#' \item{list("3")}{rand-to-best/1} \item{list("4")}{best/2}
#' \item{list("5")}{rand/2} }
#'
#' By default \code{strategy} is \code{2}. See references below for details.}
#'
#' \item{list("refresh")}{The frequency of reports. Default to every \code{10}
#' iterations.} \item{list("digits")}{The number of digits to print when
#' printing numeric values at each iteration. Default to 4.} }
#'
#' @importFrom stats runif
#' @param FUN A function to be minimized, with first argument the vector of
#' parameters over which minimization is to take place. It should return a
#' scalar result. \code{NA} and \code{NaN} values are not allowed.
#' @param lower,upper Bounds on the variables.
#' @param control A list of control parameters. See *Details*.
#' @param ... Further arguments to be passed to \code{FUN}.
#' @return A list of lists of the class \code{DEoptim}.\cr
#'
#' list \code{optim} contains the followings:\cr \code{bestmem}: the best set
#' of parameters found.\cr \code{bestval}: the value of \code{FUN}
#' corresponding to \code{bestmem}.\cr \code{nfeval}: number of function
#' evaluations.\cr \code{iter}: number of procedure iterations.\cr
#'
#' list \code{member} contains the followings:\cr \code{lower}: the lower
#' boundary.\cr \code{upper}: the upper boundary.\cr \code{bestvalit}: the best
#' value of \code{FUN} at each iteration.\cr \code{bestmemit}: the best member
#' at each iteration.\cr \code{pop}: the population generated at the last
#' iteration.\cr \code{storepop}: a list containing the intermediate
#' populations.\cr
#' @note \code{DEoptim} is a \R-vectorized variant of the Differential
#' Evolution algorithm initialy developed by Rainer Storn
#' \email{storn@icsi.berkeley.edu}, International Computer Science Institute
#' (ICSI), 1947 Center Street, Suite 600, Berkeley, CA 94704.
#'
#' If you experience misconvergence in the optimization process you usually
#' have to increase the value for \code{NP}, but often you only have to adjust
#' \code{F} to be a little lower or higher than \code{0.8}. If you increase
#' \code{NP} and simultaneously lower \code{F} a little, convergence is more
#' likely to occur but generally takes longer, i.e. \code{DEoptim} is getting
#' more robust (there is always a convergence speed/robustness tradeoff).
#'
#' \code{DEoptim} is much more sensitive to the choice of \code{F} than it is
#' to the choice of \code{CR}. \code{CR} is more like a fine tuning element.
#' High values of \code{CR} like \code{CR=1} give faster convergence if
#' convergence occurs. Sometimes, however, you have to go down as much as
#' \code{CR=0} to make \code{DEoptim} robust enough for a particular problem.
#'
#' The \R-adaptation \code{DEoptim} has properties which differ from the
#' original DE version: \describe{ \item{1.}{The random selection of vectors is
#' performed by shuffling the population array. Hence a certain vector cannot
#' be chosen twice in the same term of the perturbation expression.}
#' \item{2.}{Due to the vectorized expressions \code{DEoptim} executes fairly
#' fast.} \item{3.}{The parameters are constrained within boundaries.}
#' \item{4.}{An initial population can be given as a starting point for the
#' DE-optimization. This may speed up the convergence if the optimization
#' procedure has to be run many times for sligthly different data sets.} }
#'
#' To perform a maximization (instead of minimization) of a given function,
#' simply define a new function which is the opposite of the function to
#' maximize and apply \code{DEoptim} to it.
#'
#' To integrate additional constraints on the parameters \code{x} of
#' \code{FUN(x)}, for instance \code{x[1] + x[2]^2 < 2}, integrate the
#' constraint within the function to optimize, for instance: \preformatted{ FUN
#' <- function(x){ if (x[1] + x[2]^2 < 2){ r <- Inf else{ ...  } return(r) } }
#'
#' Note that \code{DEoptim} stops if any \code{NA} or \code{NaN} value is
#' obtained. You have to redefine your function to handle these values (for
#' instance, set \code{NA} to \code{Inf} in your objective function).
#'
#' Please cite the package in publications. Use \code{citation("DEoptim")}.
#' @author David Ardia \email{david.ardia@@unifr.ch} for the \R-port; Rainer
#' Storn \email{storn@@icsi.berkeley.edu} for the Differential Evolution
#' algorithm.
#' @seealso \code{\link{DEoptim-methods}} for methods on \code{DEoptim} object;
#' \code{\link{optim}} or \code{\link{constrOptim}} for constrained
#' optimization.
#' @references Differential Evolution homepage :
#'
#' \url{http://www.icsi.berkeley.edu/~storn/code.html}
#'
#' Some useful books:
#'
#' Price, K.V., Storn, R.M., Lampinen J.A. (2005).  \emph{Differential
#' Evolution - A Practical Approach to Global Optimization}.  Springer-Verlag.
#' ISBN 3540209506.
#'
#' Nocedal, J. and Wright, S.J. (1999).  \emph{Numerical Optimization}.
#' Springer-Verlag. ISBN 0387987932.
#' @keywords nonlinear optimize
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
#'   DEoptim(Rosenbrock, lower, upper)
#'   DEoptim(Rosenbrock, lower, upper,
#'     control = list(NP = 100, refresh = 1))
#'   DEoptim(Rosenbrock, lower, upper,
#'     control = list(NP = 50, itermax = 200, F = 1.5,
#'     CR = 0.2, refresh = 1))
#'   DEoptim(Rosenbrock, lower, upper,
#'     control = list(NP = 80, itermax = 400, F = 1.2,
#'     CR = 0.7, refresh = 1))
#'
#'   ## 'Wild' function, global minimum at about -15.81515
#'   Wild <- function(x)
#'     10 * sin(0.3*x) * sin(1.3*x^2) +
#'        0.00001 * x^4 + 0.2 * x + 80
#'   plot(Wild, -50, 50, n = 1000,
#'     main = "DEoptim minimizing 'Wild function'")
#'   DEoptim(Wild, lower = -50, upper = 50,
#'       control = list(NP = 50, refresh = 1))
#'   DEoptim(Wild, lower = -50, upper = 50,
#'       control = list(NP = 50, refresh = 1, digits = 8))
#'

DEoptim <- function(FUN, lower, upper, control = list(),
                    verbose = getOption("verbose"),...) {
  if (missing(FUN))
    stop("'FUN' is missing")
  FUN <- match.fun(FUN)

  if (missing(lower) || missing(upper))
    stop("'lower' or 'upper' is missing")
  if (length(lower) != length(upper))
    stop("'lower' and 'upper' are not of same length")
  if (!is.vector(lower))
    lower <- as.vector(lower)
  if (!is.vector(upper))
    upper <- as.vector(upper)
  if (any(lower > upper))
    stop("'lower' > 'upper'")
  if (any(lower == "Inf"))
    warning("you set a component of 'lower' to 'Inf'. May imply 'NaN' results")
  if (any(lower == "-Inf"))
    warning("you set a component of 'lower' to '-Inf'. May imply 'NaN' results")
  if (any(upper == "Inf"))
    warning("you set a component of 'upper' to 'Inf'. May imply 'NaN' results")
  if (any(upper == "-Inf"))
    warning("you set a component of 'upper' to '-Inf'. May imply 'NaN' results")

  ## Sub-functions
  fn.zeros <- function(nr, nc)
    matrix(rep.int(0, nr * nc), nrow = nr)

  fn.checkBoundaries <- function(x, lower, upper) {
    r <- apply(rbind(lower, x), 2, max)
    apply(rbind(upper, r), 2, min)
  }

  d <- length(lower)
  con <- list(VTR = -Inf, itermax = 200,
              initial = NULL,
              storepopfrom = NULL, storepopfreq = 1,
              NP = 50, F = 0.8, CR = 0.5, strategy = 2,
              refresh = 10, digits = 4)
  con[names(control)] <- control

  if (con$itermax <= 0) {
    warning("'itermax' <= 0; set to default value 200\n", immediate. = TRUE)
    con$itermax <- 200
  }
  if (con$NP < 1) {
    warning("'NP' < 1; set to default value 50\n", immediate. = TRUE)
    con$NP <- 50
  }
  NP <- con$NP
  if (con$F < 0 | con$F > 2) {
    warning("'F' not in [0,2]; set to default value 0.8\n", immediate. = TRUE)
    con$F <- 0.8
  }
  if (con$CR < 0 | con$CR > 1) {
    warning("'CR' not in [0,1]; set to default value 0.5\n", immediate. = TRUE)
    con$CR <- 0.5
  }
  if (con$strategy < 1 | con$strategy > 5) {
    warning("'strategy' not in {1,...,5}; set to default value 2\n", immediate. = TRUE)
    con$strategy <- 2
  }
  con$refresh <- floor(con$refresh)
  if (con$refresh > con$itermax)
    con$refresh <- 1

  if (is.null(con$initial)) {
    ## Initialize population and some arrays
    pop <- matrix(rep.int(lower, NP), nrow = NP, byrow = TRUE) +
      matrix(runif(NP * d), nrow = NP) *
        matrix(rep.int(upper - lower, NP), nrow = NP, byrow = TRUE)
  }
  else{
    warning("'initial' population is set by the user\n", immediate. = TRUE)
    if (!is.matrix(con$initial)){
      warning("'initial' must be a matrix; set it to a matrix\n", immediate. = TRUE)
      pop <- matrix(con$initial, nrow = NP, ncol = d)
    }
    else{
      warning("'NP' determined by the number of rows of the 'initial' population\n", immediate = TRUE)
      NP <- nrow(con$initial)
      pop <- con$initial
      if (d != ncol(pop))
        warning ("modify the length of 'lower' and 'upper' to match the dimension of 'initial'\n", immediate = TRUE)
    }
  }

  if (is.null(con$storepopfrom)) {
    con$storepopfrom <- con$itermax+1
  }

  con$storepopfreq <- floor(con$storepopfreq)
  if (con$storepopfreq > con$itermax)
    con$storepopfreq <- 1
  storepopiter <- 1
  storepop <- list()

  ## initialization
  popold <- fn.zeros(NP,d) ## toggle population
  val <- rep.int(0,NP) ## create and reset the "cost array"
  bestmem <- bestmemit <- rep.int(0,d) ## best population member ever and iteration

  ## Evaluate the best member after initialization
  nfeval <- NP ## number of function evaluations
  val <- apply(pop, 1, FUN, ...)
  if (any(is.nan(val)))
    stop ("your function returns 'NaN'; modify it or change 'lower' or 'upper' boundaries")
  if (any(is.na(val)))
    stop ("your function returns 'NA'; modify it or change 'lower' or 'upper' boundaries")

  bestval <- bestvalit <- min(val)
  ibest <- match(bestvalit, val)
  bestmem <- pop[ibest,]
  bestmemit <- matrix(bestmem, nrow = 1)

  ## DE - optimization
  ##
  ## popold is the population which has to compete. It is
  ## static through one iteration. pop is the newly emerging population.
  pm1 <- pm2 <- pm3 <- pm4 <- pm5 <- fn.zeros(NP,d) ## initialize population matrix 1 - 5
  bm <- ui <- mui <- mpo <- fn.zeros(NP,d)
  rot <- seq(from = 0, by = 1, to = (NP-1))## rotating index array (size NP)
  rotd <- seq(from = 0, by = 1, to = (d-1)) ## rotating index array (size d)
  rt <- fn.zeros(NP,NP) ## another rotating index array
  rtd <- fn.zeros(d,d) ## rotating index array for exponential crossover
  a1 <- a2 <- a3 <- a4 <- a5 <- fn.zeros(NP,NP) ## index array 1 - 5
  ind <- fn.zeros(4,4)

  iter <- 1
  while (iter <= con$itermax & bestval >= con$VTR){
    popold <- pop ## save old population

    ind <- sample(1:4) ## index pointer array

    a1 <- sample(1:NP) ## shuffle locations and rotate vectors
    rt <- (rot + ind[1]) %% NP
    a2 <- a1[rt + 1]
    rt <- (rot + ind[2]) %% NP
    a3 <- a2[rt + 1]
    rt <- (rot + ind[3]) %% NP
    a4 <- a3[rt + 1]
    rt <- (rot + ind[4]) %% NP
    a5 <- a4[rt + 1]

    pm1 <- popold[a1,] ## shuffled populations 1 - 5
    pm2 <- popold[a2,]
    pm3 <- popold[a3,]
    pm4 <- popold[a4,]
    pm5 <- popold[a5,]

    bm <- matrix(rep.int(bestmemit[iter,], NP), nrow = NP, byrow = TRUE) ## population filled with
    ## the best member of the last iteration

    mui <- matrix(runif(NP * d), nrow = NP) < con$CR ## all random numbers < CR are 1, 0 otherwise
    mpo <- mui < 0.5

    if (con$strategy == 1) { ## best / 1
      ui <- bm + con$F * (pm1 - pm2) ## differential variation
      ui <- popold * mpo + ui * mui ## crossover
    }
    else if (con$strategy == 2) { ## rand / 1
      ui <- pm3 + con$F * (pm1 - pm2) ## differential variation
      ui <- popold * mpo + ui * mui ## crossover
    }
    else if (con$strategy == 3) { ## rand-to-best / 1
      ui <- popold + con$F * (bm - popold) + con$F * (pm1 - pm2) ## differential variation
      ui <- popold * mpo + ui * mui ## crossover
    }
    else if (con$strategy == 4) { ## best / 2
      ui <- bm + con$F * (pm1 - pm2 + pm3 - pm4) ## differential variation
      ui <- popold * mpo + ui * mui ## crossover
    }
    else { ## rand / 2
      ui <- pm5 + con$F * (pm1 - pm2 + pm3 - pm4) ## differential variation
      ui <- popold * mpo + ui * mui ## crossover
    }

    for (i in 1:NP)
      ui[i,] <- fn.checkBoundaries(ui[i,], lower, upper) ## check whether
    ## the components are within the boundaries

    nfeval <- nfeval + NP
    tempval <- apply(ui, 1, FUN, ...) ## check cost of competitor
    if (any(is.nan(tempval)))
      stop ("'your function returns 'NaN'; modify it or change 'lower' or 'upper' boundaries")
    if (any(is.na(tempval)))
      stop ("your function returns 'NA'; modify it or change 'lower' or 'upper' boundaries")
    ichange <- tempval <= val
    val[ichange] <- tempval[ichange]
    pop[ichange,] <- ui[ichange,]
    bestval <- min(val)
    bestvalit <- c(bestvalit, bestval)
    ibest <- match(bestval, val)
    bestmem <- pop[ibest,]
    bestmemit <- rbind(bestmemit, bestmem)

    ## keeppop
    if (iter >= con$storepopfrom & iter %% con$storepopfreq == 0){
      storepop[[storepopiter]] <- pop
      storepopiter <- storepopiter + 1
    }

    ## refresh output
    if (verbose && con$refresh > 0 & iter %% con$refresh == 0) {
      message("iteration: ", iter,
          "best member: " , round(bestmem, con$digits),
          "best value: ", round(bestval, con$digits), "\n")
    }
    iter <- iter + 1
  }

  if (!is.null(names(lower)))
    nam <- names(lower)
  else if (!is.null(names(upper)) & is.null(names(lower)))
    nam <- names(upper)
  else
    nam <- paste("par", 1:length(lower), sep = "")

  names(lower) <- names(upper) <- names(bestmem) <- nam
  dimnames(bestmemit) <- list(1:iter, nam)
  r <- list(optim = list(
              bestmem = bestmem,
              bestval = bestval,
              nfeval = nfeval,
              iter = iter-1),
            member = list(
              lower = lower,
              upper = upper,
              bestvalit = bestvalit,
              bestmemit = bestmemit,
              pop = pop,
              storepop = storepop))

  attr(r, "class") <- "DEoptim"
  return(r)
}

summary.DEoptim <- function(object, ...){
  digits <- max(5, getOption('digits') - 2)
  z <- object$optim

  cat("\n***** summary of DEoptim object *****",
      "\nbest member   : ", round(z$bestmem, digits),
      "\nbest value    : ", round(z$bestval, digits),
      "\nafter         : ", round(z$iter), "iterations",
      "\nFUN evaluated : ", round(z$nfeval), "times",
      "\n*************************************\n")

  invisible(z)
}

#' @importFrom graphics abline
#' @importFrom graphics points
#' @importFrom graphics par
#' @importFrom grDevices gray

plot.DEoptim <- function(x, plot.type = c("bestmemit","bestvalit","storepop"), ...){
  plot.type <- plot.type[1]
  z <- x$member
  n <- length(z$bestvalit)
  npar <- length(z$lower)
  nam <- names(z$lower)
  k <- length(z$storepop)
  if (plot.type == "bestmemit") {
    if (npar == 1) {
      plot(1:n, z$bestmemit,
           xlab = "iteration", ylab = "value", main = nam, ...)
      abline(h = c(z$lower, z$upper), col = 'red')
    }
    else if (npar == 2) {
      plot(z$bestmemit[,1], z$bestmemit[,2],
           xlab = nam[1], ylab = nam[2], ...)
      abline(h = c(z$lower[1], z$upper[1]), col = 'red')
      abline(v = c(z$lower[2], z$upper[2]), col = 'red')
    }
    else{
      par(mfrow = c(npar,1))
      for (i in 1:npar){
        plot(1:n, z$bestmemit[,i],
             xlab = "iteration", ylab = "value", main = nam[i], ...)
        abline(h = c(z$lower[i], z$upper[i]), col = 'red')
      }
    }
  }
  else if (plot.type == "bestvalit") {
    plot(1:n, z$bestvalit,
         xlab = "iteration", ylab = "function value",
         main = "convergence plot", ...)
  }
  else if (plot.type == "storepop" & k>0) {
    the.col <- gray(k:1/k)
    if (npar == 1) {
      for (i in 1:k) {
        plot(rep(i/k, nrow(z$storepop[[i]])),
             z$storepop[[i]], xlim = c(0,1), ylim = c(z$lower, z$upper),
             las=1, pch=20, col=the.col[i], xlab = nam[1], ylab = "", main = "")
        par(new=TRUE)
      }
    }
    else {
      if (npar>2)
        par(mfrow=c(npar-1,npar-1))
      for (i in 1:k) {
        for (j in 1:(npar-1)) {
          for (l in (j+1):npar) {
            par(mfg=c(j,l-1))
            plot(z$storepop[[i]][,c(j,l)], col=the.col[i], las=1, pch=20,
                 xlim = c(z$lower[j],z$upper[j]), ylim = c(z$lower[l],z$upper[l]),
                 xlab = nam[j], ylab = nam[l], main = "", ...)
            par(new = TRUE)
          }
        }
      }
    }
    par(new = FALSE)
  }
  else {
    warning("'plot.type' does not correspond to any type", immediate. = TRUE)
  }
}
