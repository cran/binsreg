######################################################################################
#'@title Data-driven IMSE-Optimal Partitioning/Binning Selection for Binscatter
#'@description \code{binsregselect} implements data-driven procedures for selecting the number of bins for binscatter
#'             estimation. The selected number is optimal in minimizing integrated mean squared error (IMSE).
#'@param  y outcome variable. A vector.
#'@param  x independent variable of interest. A vector.
#'@param  w control variables. A matrix or a vector.
#'@param  deriv  derivative order of the regression function for estimation, testing and plotting.
#'               The default is \code{deriv=0}, which corresponds to the function itself.
#'@param  bins a vector. \code{bins=c(p,s)} set a piecewise polynomial of degree \code{p} with \code{s} smoothness constraints
#'             for data-driven (IMSE-optimal) selection of the partitioning/binning scheme. The default is
#'             \code{bins=c(0, 0)}, which corresponds to piecewise constant (canonical binscatter).
#'@param  binspos position of binning knots.  The default is \code{binspos="qs"}, which corresponds to quantile-spaced
#'                binning (canonical binscatter).  The other options is \code{"es"} for evenly-spaced binning.
#'@param  binsmethod method for data-driven selection of the number of bins. The default is \code{binsmethod="dpi"},
#'                   which corresponds to the IMSE-optimal direct plug-in rule.  The other option is: \code{"rot"}
#'                   for rule of thumb implementation.
#'@param  nbinsrot initial number of bins value used to construct the DPI number of bins selector.
#'                 If not specified, the data-driven ROT selector is used instead.
#'@param  simsgrid number of evaluation points of an evenly-spaced grid within each bin used for evaluation of
#'                 the supremum (or infimum) operation needed to construct confidence bands and hypothesis testing
#'                 procedures. The default is \code{simsgrid=20}, which corresponds to 20 evenly-spaced
#'                 evaluation points within each bin for approximating the supremum (or infimum) operator.
#'@param  savegrid If true, a data frame produced containing grid.
#'@param  vce procedure to compute the variance-covariance matrix estimator. Options are
#'           \itemize{
#'           \item \code{"const"} homoskedastic variance estimator.
#'           \item \code{"HC0"} heteroskedasticity-robust plug-in residuals variance estimator
#'                              without weights.
#'           \item \code{"HC1"} heteroskedasticity-robust plug-in residuals variance estimator
#'                              with hc1 weights. Default.
#'           \item \code{"HC2"} heteroskedasticity-robust plug-in residuals variance estimator
#'                              with hc2 weights.
#'           \item \code{"HC3"} heteroskedasticity-robust plug-in residuals variance estimator
#'                              with hc3 weights.
#'           }
#'@param  useeffn effective sample size to be used when computing the (IMSE-optimal) number of bins. This option
#'                is useful for extrapolating the optimal number of bins to larger (or smaller) datasets than
#'                the one used to compute it.
#'@param  cluster cluster ID. Used for compute cluster-robust standard errors.
#'@param  dfcheck adjustments for minimum effective sample size checks, which take into account number of unique
#'                values of \code{x} (i.e., number of mass points), number of clusters, and degrees of freedom of
#'                the different stat models considered. The default is \code{dfcheck=c(20, 30)}.
#'                See \href{https://arxiv.org/abs/1902.09615}{Cattaneo, Crump, Farrell and Feng (2019b)} for more details.
#'@param  masspoints how mass points in \code{x} are handled. Available options:
#'                   \itemize{
#'                   \item \code{"on"} all mass point and degrees of freedom checks are implemented. Default.
#'                   \item \code{"noadjust"} mass point checks and the corresponding effective sample size adjustments are omitted.
#'                   \item \code{"nolocalcheck"} within-bin mass point and degrees of freedom checks are omitted.
#'                   \item \code{"off"} "noadjust" and "nolocalcheck" are set simultaneously.
#'                   \item \code{"veryfew"} forces the function to proceed as if \code{x} has only a few number of mass points (i.e., distinct values).
#'                                          In other words, forces the function to proceed as if the mass point and degrees of freedom checks were failed.
#'                   }
#'@param  norotnorm if true, a uniform density rather than normal density used for ROT selection.
#'@param  numdist  number of distinct for selection. Used to speed up computation.
#'@param  numclust number of clusters for selection. Used to speed up computation.
#'@param  weights an optional vector of weights to be used in the fitting process. Should be \code{NULL} or
#'                a numeric vector. For more details, see \code{\link{lm}}.
#'@param  subset optional rule specifying a subset of observations to be used.
#'@return \item{\code{nbinsrot.poly}}{ROT number of bins, unregularized.}
#'        \item{\code{nbinsrot.regul}}{ROT number of bins, regularized.}
#'        \item{\code{nbinsrot.uknot}}{ROT number of bins, unique knots.}
#'        \item{\code{nbinsdpi}}{DPI number of bins.}
#'        \item{\code{nbinsdpi.uknot}}{DPI number of bins, unique knots.}
#'        \item{\code{opt}}{ A list containing options passed to the function, as well as total sample size \code{n},
#'                           number of distinct values \code{Ndist} in \code{x}, and number of clusters \code{Nclust}.}
#'        \item{\code{data.grid}}{A data frame containing grid.}
#'@author
#' Matias D. Cattaneo, University of Michigan, Ann Arbor, MI. \email{cattaneo@umich.edu}.
#'
#' Richard K. Crump, Federal Reserve Bank of New York, New York, NY. \email{richard.crump@ny.frb.org}.
#'
#' Max H. Farrell, University of Chicago, Chicago, IL. \email{max.farrell@chicagobooth.edu}.
#'
#' Yingjie Feng (maintainer), University of Michigan, Ann Arbor, MI. \email{yjfeng@umich.edu}.
#'
#'@references
#' Cattaneo, M. D., R. K. Crump, M. H. Farrell, and Y. Feng. 2019a: \href{https://arxiv.org/abs/1902.09608}{On Binscatter}. Working Paper.
#'
#' Cattaneo, M. D., R. K. Crump, M. H. Farrell, and Y. Feng. 2019b: \href{https://arxiv.org/abs/1902.09615}{Binscatter Regressions}. Working Paper.
#'
#'@seealso \code{\link{binsreg}}, \code{\link{binsregtest}}.
#'
#'@examples
#'  x <- runif(500); y <- sin(x)+rnorm(500)
#'  est <- binsregselect(y,x)
#'  summary(est)
#'@export

binsregselect <- function(y, x, w=NULL, deriv=0,
                          bins=c(0,0), binspos="qs", binsmethod="dpi", nbinsrot=NULL,
                          simsgrid=20, savegrid=F,
                          vce="HC1", useeffn=NULL, cluster=NULL,
                          dfcheck=c(20,30), masspoints="on", weights=NULL, subset=NULL,
                          norotnorm=F, numdist=NULL, numclust=NULL) {

  # param for internal use
  rot.lb <- 1
  qrot <- 2

  ####################
  ### prepare data ###
  ####################
  if (is.data.frame(y)) y <- y[,1]
  xname <- NULL
  if (is.data.frame(x)) {
    xname <- colnames(x)[1]; x <- x[,1]
  }
  if (!is.null(w))      w <- as.matrix(w)

  # substract subset
  if (!is.null(subset)) {
    y <- y[subset]
    x <- x[subset]
    w <- w[subset, , drop = F]
  }

  na.ok <- complete.cases(x) & complete.cases(y)
  if (!is.null(w)) na.ok <- na.ok & complete.cases(w)

  y <- y[na.ok]
  x <- x[na.ok]
  w <- w[na.ok, , drop = F]
  xmin <- min(x); xmax <- max(x)

  #############################error checking
  exit <- 0
  if (length(bins)==2) if (bins[1]<bins[2]) {
    print("p<s not allowed.")
    exit <- 1
  }
  if (binsmethod!="dpi" & binsmethod!="rot") {
    print("bin selection method incorrectly specified.")
    exit <- 1
  }
  if (binspos!="es" & binspos!="qs") {
    print("binspos incorrectly specified.")
    exit <- 1
  }
  if (exit>0) stop()

  #####################################
  rot.fewobs <- dpi.fewobs <- F
  localcheck <- massadj <- T
  if (masspoints=="on") {
    localcheck <- T; massadj <- T
  } else if (masspoints=="off") {
    localcheck <- F; massadj <- F
  } else if (masspoints=="noadjust") {
    localcheck <- T; massadj <- F
  } else if (masspoints=="nolocalcheck") {
    localcheck <- F; massadj <- T
  } else if (masspoints=="veryfew") {
    rot.fewobs <- dpi.fewobs <- T
  }

  # effective size
  eN <- N <- length(x)
  Ndist <- NA
  if (massadj) {
    if (!is.null(numdist)) {
      Ndist <- numdist
    } else {
      Ndist <- length(unique(x))
    }
    eN <- min(eN, Ndist)
  }
  Nclust <- NA
  if (!is.null(cluster)) {
    if (!is.null(numclust)) {
      Nclust <- numclust
    } else {
      Nclust <- length(unique(cluster))
    }
    eN <- min(eN, Nclust)
  }

  # Prepare params
  p <- bins[1]; s <- bins[2]

  if (binspos == "es") {
    es <- T
  } else {
    es <- F
  }

  # Store options
  if (es)  {
    position <- "Evenly-spaced"
  } else {
    position <- "Quantile-spaced"
  }

  if (binsmethod == "dpi") {
    selectmethod <- "IMSE direct plug-in"
  } else {
    selectmethod <- "IMSE rule-of-thumb"
  }

  # Run rot selection
  J.rot.regul <- J.rot.poly <- NA
  if (!is.null(nbinsrot)) J.rot.regul <- nbinsrot
  if (is.na(J.rot.regul) & !rot.fewobs) {
    # checking
    if (eN <= dfcheck[1]+p+1+qrot) {
      rot.fewobs <- T
      warning("too small effective sample size for bin selection.")
    }

    if (!rot.fewobs) {
      J.rot.poly <- binsregselect.rot(y, x, w, p, s, deriv, es=es, eN=eN, qrot=qrot, norotnorm=norotnorm, weights=weights)
    }
    J.rot.regul <- max(J.rot.poly, ceiling((2*(p+1-deriv)/(1+2*deriv)*rot.lb*eN)^(1/(2*p+3))))
  }
  # repeated knots?
  J.rot.uniq <- J.rot.regul
  if (!es & !is.na(J.rot.regul)) {
     J.rot.uniq <- length(unique(genKnot.qs(x, J.rot.regul)[-1]))
  }

  # Run dpi selection
  J.dpi <- NA
  if (binsmethod == "dpi" & !dpi.fewobs) {
    # check if dpi can be implemented
    if (!is.na(J.rot.uniq)) {
      if ((p-s+1)*(J.rot.uniq-1)+p+2+dfcheck[2]>=eN) {
         dpi.fewobs <- T
         warning("too small effective sample size for DPI selection.")
      }

      # check empty bins
      if (localcheck) {
         uniqmin <- binsreg.checklocalmass(x, J.rot.regul, es, knot=NULL) # mimic STATA
         if (uniqmin < p+2) {
           dpi.fewobs <- T
           warning("some bins have too few distinct values of x for DPI selection.")
         }
      }
    } else {
      dpi.fewobs <- T
    }

    if (!dpi.fewobs) {
      if (massadj) {
        if (is.null(cluster)) {
          cluster <- x
          # note: clustered at mass point level
        } else {
          if (Nclust > Ndist) {
            cluster <- x
            warning("# mass points < # clusters. Clustered at mass point level.")
          }
        }
      }
      J.dpi <- binsregselect.dpi(y, x, w, p, s, deriv, es=es, vce=vce, cluster=cluster, nbinsrot=J.rot.uniq, weights=weights)
    }
  }
  J.dpi.uniq <- J.dpi

  if (!is.null(useeffn)) {
    scaling <- (useeffn/eN)^(1/(2*p+2+1))
    if (!is.na(J.rot.poly))  J.rot.poly  <- J.rot.poly * scaling
    if (!is.na(J.rot.regul)) J.rot.regul <- J.rot.regul * scaling
    if (!is.na(J.rot.uniq))  J.rot.uniq  <- J.rot.uniq * scaling
    if (!is.na(J.dpi))       J.dpi       <- J.dpi * scaling
    if (!is.na(J.dpi.uniq))  J.dpi.uniq  <- J.dpi.uniq * scaling
  }

  # Generate a knot vector
  if (binsmethod == "rot") {
    Jselect <- J.rot.uniq
  } else {
    Jselect <- J.dpi
  }
  knot <- NA; data.grid <- NA
  if (!is.na(Jselect) & is.null(useeffn)) {
    if (es) {
      knot <- genKnot.es(xmin, xmax, Jselect)
    } else {
      knot <- genKnot.qs(x, Jselect)
    }
    knot <- c(knot[1], unique(knot[-1]))
    Jselect <- length(knot)-1
    if (binsmethod=="dpi") {
      J.dpi.uniq <- Jselect
    }

    # a grid dataset
    if (savegrid) {
      grid <- binsreg.grid(knot=knot, ngrid=simsgrid, addmore=T)
      data.grid <- cbind(grid$eval, grid$bin, grid$isknot)
      if (!is.null(w)) {
         data.grid <- cbind(data.grid, matrix(0, nrow(data.grid), ncol(w)))
      }
      data.grid <- data.frame(data.grid)
      if (is.null(xname)) {
         colnames(data.grid) <- c("x", "binreg_bin", "binsreg_isknot", colnames(w))
      } else {
         colnames(data.grid) <- c(xname, "binreg_bin", "binsreg_isknot", colnames(w))
      }
    }
  }

  ######################
  #######output#########
  ######################
  out <- list(nbinsrot.poly=J.rot.poly, nbinsrot.regul=J.rot.regul, nbinsrot.uknot=J.rot.uniq,
              nbinsdpi=J.dpi, nbinsdpi.uknot=J.dpi.uniq,
              opt = list(bins.p=p, bins.s=s, deriv=deriv,
                         binspos=position, binsmethod=selectmethod,
                         n=N, Ndist=Ndist, Nclust=Nclust),
              knot=knot, data.grid=data.grid)
  out$call   <- match.call()
  class(out) <- "CCFFbinsregselect"
  return(out)
}

##########################################################################
#' Internal function.
#'
#' @param x Class \code{CCFFbinsregselect} objects.
#'
#' @keywords internal
#' @export
#'
print.CCFFbinsregselect <- function(x, ...) {
  cat("Call: binsregselect\n\n")

  cat(paste("Sample size (n)                    =  ", x$opt$n,          "\n", sep=""))
  cat(paste("# of distinct values (Ndist)       =  ", x$opt$Ndist,      "\n", sep=""))
  cat(paste("# of clusters (Nclust)             =  ", x$opt$Nclust,     "\n", sep=""))
  cat(paste("Derivative (deriv)                 =  ", x$opt$deriv,      "\n", sep=""))
  cat(paste("Bin selection:", "\n"))
  cat(paste("  Method (binsmethod)              =  ", x$opt$binsmethod, "\n", sep=""))
  cat(paste("  Placement (binspos)              =  ", x$opt$binspos,    "\n", sep=""))
  cat(paste("  degree (p)                       =  ", x$opt$bins.p,     "\n", sep=""))
  cat(paste("  smooth (s)                       =  ", x$opt$bins.s,     "\n", sep=""))
  cat(paste("  # of bins (ROT-POLY)             =  ", sprintf("%1.0f",x$nbinsrot.poly),  "\n", sep=""))
  cat(paste("  # of bins (ROT-REGUL)            =  ", sprintf("%1.0f",x$nbinsrot.regul), "\n", sep=""))
  cat(paste("  # of bins (ROT-UKNOT)            =  ", sprintf("%1.0f",x$nbinsrot.uknot), "\n", sep=""))
  if (x$opt$binsmethod=="IMSE direct plug-in") {
  cat(paste("  # of bins (DPI)                  =  ", sprintf("%1.0f",x$nbinsdpi),       "\n", sep=""))
  cat(paste("  # of bins (DPI-UKNOT)            =  ", sprintf("%1.0f",x$nbinsdpi.uknot), "\n", sep=""))
  }
  cat("\n")
}

################################################################################
#' Internal function.
#'
#' @param object Class \code{CCFFbinsregselect} objects.
#'
#' @keywords internal
#' @export
summary.CCFFbinsregselect <- function(object, ...) {
  x <- object
  args <- list(...)

  cat("Call: binsregselect\n\n")

  cat(paste("Sample size (n)                    =  ", x$opt$n,          "\n", sep=""))
  cat(paste("# of distinct values (Ndist)       =  ", x$opt$Ndist,      "\n", sep=""))
  cat(paste("# of clusters (Nclust)             =  ", x$opt$Nclust,     "\n", sep=""))
  cat(paste("Derivative (deriv)                 =  ", x$opt$deriv,      "\n", sep=""))
  cat(paste("Bin selection:", "\n"))
  cat(paste("  Method (binsmethod)              =  ", x$opt$binsmethod, "\n", sep=""))
  cat(paste("  Placement (binspos)              =  ", x$opt$binspos,    "\n", sep=""))
  cat(paste("  degree (p)                       =  ", x$opt$bins.p,     "\n", sep=""))
  cat(paste("  smooth (s)                       =  ", x$opt$bins.s,     "\n", sep=""))
  cat(paste("  # of bins (ROT-POLY)             =  ", sprintf("%1.0f",x$nbinsrot.poly),  "\n", sep=""))
  cat(paste("  # of bins (ROT-REGUL)            =  ", sprintf("%1.0f",x$nbinsrot.regul), "\n", sep=""))
  cat(paste("  # of bins (ROT-UKNOT)            =  ", sprintf("%1.0f",x$nbinsrot.uknot), "\n", sep=""))
  if (x$opt$binsmethod=="IMSE direct plug-in") {
  cat(paste("  # of bins (DPI)                  =  ", sprintf("%1.0f",x$nbinsdpi),       "\n", sep=""))
  cat(paste("  # of bins (DPI-UKNOT)            =  ", sprintf("%1.0f",x$nbinsdpi.uknot), "\n", sep=""))
  }
  cat("\n")

  cat(paste(rep("=", 15 + 15 + 15), collapse="")); cat("\n")
  cat(format("method",     width=15, justify="right"))
  cat(format("# of bins",  width=15, justify="right"))
  cat(format("df",         width=15, justify="right"))
  cat("\n")

  cat(paste(rep("-", 15 + 15 + 15), collapse="")); cat("\n")
  cat(format("ROT-POLY",  width= 15, justify="right"))
  cat(format(sprintf("%3.0f", x$nbinsrot.poly), width=15 , justify="right"))
  ROT.poly.df <- NA
  if (!is.na(x$nbinsrot.poly)) ROT.poly.df <- x$opt$bins.p+1+(x$nbinsrot.poly-1)*(x$opt$bins.p-x$opt$bins.s+1)
  cat(format(sprintf("%3.0f", ROT.poly.df), width=15, justify="right"))
  cat("\n")

  cat(format("ROT-REGUL",  width= 15, justify="right"))
  cat(format(sprintf("%3.0f", x$nbinsrot.regul), width=15 , justify="right"))
  ROT.regul.df <- NA
  if (!is.na(x$nbinsrot.regul)) ROT.regul.df <- x$opt$bins.p+1+(x$nbinsrot.regul-1)*(x$opt$bins.p-x$opt$bins.s+1)
  cat(format(sprintf("%3.0f", ROT.regul.df), width=15, justify="right"))
  cat("\n")

  cat(format("ROT-UKNOT",  width= 15, justify="right"))
  cat(format(sprintf("%3.0f", x$nbinsrot.uknot), width=15 , justify="right"))
  ROT.uknot.df <- NA
  if (!is.na(x$nbinsrot.uknot)) ROT.uknot.df <- x$opt$bins.p+1+(x$nbinsrot.uknot-1)*(x$opt$bins.p-x$opt$bins.s+1)
  cat(format(sprintf("%3.0f", ROT.uknot.df), width=15, justify="right"))
  cat("\n")

  if (x$opt$binsmethod=="IMSE direct plug-in") {
    cat(format("DPI",  width=15, justify="right"))
    cat(format(sprintf("%3.0f", x$nbinsdpi), width=15 , justify="right"))
    DPI.df <- NA
    if (!is.na(x$nbinsdpi)) DPI.df <- x$opt$bins.p+1+(x$nbinsdpi-1)*(x$opt$bins.p-x$opt$bins.s+1)
    cat(format(sprintf("%3.0f", DPI.df), width=15, justify="right"))
    cat("\n")

    cat(format("DPI-UKNOT",  width=15, justify="right"))
    cat(format(sprintf("%3.0f", x$nbinsdpi.uknot), width=15 , justify="right"))
    DPI.uknot.df <- NA
    if (!is.na(x$nbinsdpi.uknot)) DPI.uknot.df <- x$opt$bins.p+1+(x$nbinsdpi.uknot-1)*(x$opt$bins.p-x$opt$bins.s+1)
    cat(format(sprintf("%3.0f", DPI.uknot.df), width=15, justify="right"))
    cat("\n")
  }
  cat(paste(rep("-", 15 + 15 + 15), collapse=""))
  cat("\n")
}
