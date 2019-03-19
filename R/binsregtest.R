################################################################################################################
#'@title  Data-driven Nonparametric Shape Restriction and Parametric Model Specification Testing using Binscatter
#'@description \code{binsregtest} implements binscatter-based hypothesis testing procedures for parametric functional
#'             forms of and nonparametric shape restrictions on the regression function estimators, following the results
#'             in \href{https://arxiv.org/abs/1902.09608}{Cattaneo, Crump, Farrell and Feng (2019a)}.
#'             If the binning scheme is not set by the user,
#'             the companion function \code{\link{binsregselect}} is used to implement binscatter in a
#'             data-driven (optimal) way and inference procedures are based on robust bias correction.
#'             Binned scatter plots can be constructed using the companion function \code{\link{binsreg}}.
#'@param  y outcome variable. A vector.
#'@param  x independent variable of interest. A vector.
#'@param  w control variables. A matrix or a vector.
#'@param  deriv  derivative order of the regression function for estimation, testing and plotting.
#'               The default is \code{deriv=0}, which corresponds to the function itself.
#'@param  testmodel a vector. \code{testmodel=c(p,s)} sets a piecewise polynomial of degree \code{p} with \code{s}
#'                  smoothness constraints for parametric model specification testing.  The default is
#'                  \code{testmodel=c(3,3)}, which corresponds to a cubic B-spline estimate of the regression
#'                  function of interest for testing against the fitting from a parametric model specification.
#'@param  testmodelparfit a data frame or matrix which contains the evaluation grid and fitted values of the model(s) to be
#'                        tested against.  The column contains a series of evaluation points
#'                        at which the binscatter model and the parametric model of interest are compared with
#'                        each other.  Each parametric model is represented by other columns, which must
#'                        contain the fitted values at the corresponding evaluation points.
#'@param  testmodelpoly degree of a global polynomial model to be tested against.
#'@param  testshape a vector. \code{testshape=c(p,s)} sets a piecewise polynomial of degree \code{p} with \code{s}
#'                  smoothness constraints for nonparametric shape restriction testing. The default is
#'                  \code{testshape=c(3,3)}, which corresponds to a cubic B-spline estimate of the regression
#'                  function of interest for one-sided or two-sided testing.
#'@param  testshapel a vector of null boundary values for hypothesis testing. Each number \code{a} in the vector
#'                   corresponds to one boundary of a one-sided hypothesis test to the left of the form
#'                   \code{H0: sup_x mu(x)<=a}.
#'@param  testshaper a vector of null boundary values for hypothesis testing. Each number \code{a} in the vector
#'                   corresponds to one boundary of a one-sided hypothesis test to the right of the form
#'                   \code{H0: inf_x mu(x)>=a}.
#'@param  testshape2 a vector of null boundary values for hypothesis testing. Each number \code{a} in the vector
#'                   corresponds to one boundary of a two-sided hypothesis test ofthe form
#'                   \code{H0: sup_x |mu(x)-a|=0}.
#'@param  bins Degree and smoothness for bin selection.
#'@param  nbins number of bins for partitioning/binning of \code{x}.  If not specified, the number of bins is
#'              selected via the companion function \code{binsregselect} in a data-driven, optimal way whenever possible.
#'@param  binspos position of binning knots. The default is \code{binspos="qs"}, which corresponds to quantile-spaced
#'                binning (canonical binscatter).  The other options are \code{"es"} for evenly-spaced binning, or
#'                a vector for manual specification of the positions of inner knots (which must be within the range of
#'                \code{x}).
#'@param  binsmethod method for data-driven selection of the number of bins. The default is \code{binsmethod="dpi"},
#'                   which corresponds to the IMSE-optimal direct plug-in rule.  The other option is: \code{"rot"}
#'                   for rule of thumb implementation.
#'@param  nbinsrot initial number of bins value used to construct the DPI number of bins selector.
#'                 If not specified, the data-driven ROT selector is used instead.
#'@param  nsims number of random draws for constructing confidence bands and hypothesis testing. The default is
#'              \code{nsims=500}, which corresponds to 500 draws from a standard Gaussian random vector of size
#'              \code{[(p+1)*J - (J-1)*s]}.
#'@param  simsgrid number of evaluation points of an evenly-spaced grid within each bin used for evaluation of
#'                 the supremum (or infimum) operation needed to construct confidence bands and hypothesis testing
#'                 procedures. The default is \code{simsgrid=20}, which corresponds to 20 evenly-spaced
#'                 evaluation points within each bin for approximating the supremum (or infimum) operator.
#'@param  simsseed  seed for simulation.
#'@param  vce Procedure to compute the variance-covariance matrix estimator. Options are
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
#'@param  cluster cluster ID. Used for compute cluster-robust standard errors.
#'@param  dfcheck adjustments for minimum effective sample size checks, which take into account number of unique
#'                values of \code{x} (i.e., number of mass points), number of clusters, and degrees of freedom of
#'                the different stat models considered. The default is \code{dfcheck=c(20, 30)}.
#'                See \href{https://sites.google.com/site/nppackages/binsreg/Cattaneo-Crump-Farrell-Feng_2019_Stata.pdf}{Cattaneo, Crump, Farrell and Feng (2019b)} for more details.
#'@param  masspoints how mass points in \code{x} are handled. Available options:
#'                   \itemize{
#'                   \item \code{"on"} all mass point and degrees of freedom checks are implemented. Default.
#'                   \item \code{"noadjust"} mass point checks and the corresponding effective sample size adjustments are omitted.
#'                   \item \code{"nolocalcheck"} within-bin mass point and degrees of freedom checks are omitted.
#'                   \item \code{"off"} "noadjust" and "nolocalcheck" are set simultaneously.
#'                   \item \code{"veryfew"} forces the function to proceed as if \code{x} has only a few number of mass points (i.e., distinct values).
#'                                          In other words, forces the function to proceed as if the mass point and degrees of freedom checks were failed.
#'                   }
#'@param  weights an optional vector of weights to be used in the fitting process. Should be \code{NULL} or
#'                a numeric vector. For more details, see \code{\link{lm}}.
#'@param  subset Optional rule specifying a subset of observations to be used.
#'@param  numdist  Number of distinct for selection. Used to speed up computation.
#'@param  numclust Number of clusters for selection. Used to speed up computation.
#'@return \item{\code{testshapeL}}{Results for \code{testshapel}, including: \code{testvalL}, null boundary values;
#'                                 \code{stat.shapeL}, test statistics; and \code{pval.shapeL}, p-value.}
#'        \item{\code{testshapeR}}{Results for \code{testshaper}, including: \code{testvalR}, null boundary values;
#'                                 \code{stat.shapeR}, test statistics; and \code{pval.shapeR}, p-value.}
#'        \item{\code{testshape2}}{Results for \code{testshape2}, including: \code{testval2}, null boundary values;
#'                                 \code{stat.shape2}, test statistics; and \code{pval.shape2}, p-value.}
#'        \item{\code{testpoly}}{Results for \code{testmodelpoly}, including: \code{testpoly}, the degree of global polynomial;
#'                               \code{stat.poly}, test statistic; \code{pval.poly}, p-value.}
#'        \item{\code{testmodel}}{Results for \code{testmodelparfit}, including: \code{stat.model}, test statistics;
#'                                \code{pval.model}, p-values.}
#'        \item{\code{opt}}{ A list containing options passed to the function, as well as total sample size \code{n},
#'                           number of distinct values \code{Ndist} in \code{x}, number of clusters \code{Nclust}, and
#'                           number of bins \code{nbins}.}
#'
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
#'@seealso \code{\link{binsreg}}, \code{\link{binsregselect}}.
#'
#'@examples
#'  x <- runif(500); y <- sin(x)+rnorm(500)
#'  est <- binsregtest(y,x, testmodelpoly=1)
#'  summary(est)
#'@export

binsregtest <- function(y, x, w=NULL, deriv=0,
                        testmodel=c(3,3), testmodelparfit=NULL, testmodelpoly=NULL,
                        testshape=c(3,3), testshapel=NULL,
                        testshaper=NULL, testshape2=NULL,
                        bins=c(0,0), nbins=NULL, binspos="qs",
                        binsmethod="dpi", nbinsrot=NULL,
                        nsims=500, simsgrid=20, simsseed=666,
                        vce="HC1", cluster=NULL,
                        dfcheck=c(20,30), masspoints="on", weights=NULL, subset=NULL,
                        numdist=NULL, numclust=NULL) {

  # param for internal use
  qrot <- 2

  ####################
  ### prepare data ###
  ####################
  if (is.data.frame(y)) y <- y[,1]
  if (is.data.frame(x)) x <- x[,1]
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


  ##################################error checking
  exit <- 0
  if (!is.character(binspos)) {
    if (min(binspos)<=xmin|max(binspos)>=xmax) {
      print("knots out of allowed range")
      exit <- 1
    }
  } else {
    if (binspos!="es"&binspos!="qs") {
      print("binspos incorrectly specified.")
      exit <- 1
    }
  }
  if (length(bins)==2) if (bins[1]<bins[2]) {
    print("p<s not allowed.")
    exit <- 1
  }
  if (length(testshape)==2) if (testshape[1]<testshape[2]) {
    print("p<s not allowed.")
    exit <- 1
  }
  if (length(testmodel)==2) if (testmodel[1]<testmodel[2]) {
    print("p<s not allowed.")
    exit <- 1
  }
  if (testshape[1]<deriv|testmodel[1]<deriv) {
    print("p<deriv not allowed.")
    exit <- 1
  }
  if (binsmethod!="dpi" & binsmethod!="rot") {
    print("bin selection method incorrectly specified.")
    exit <- 1
  }
  if (masspoints == "veryfew") {
    print("veryfew not allowed for testing.")
    exit <- 1
  }
  if (exit>0) stop()

  #################################################
  localcheck <- massadj <- T
  if (masspoints=="on") {
    localcheck <- T; massadj <- T
  } else if (masspoints=="off") {
    localcheck <- F; massadj <- F
  } else if (masspoints=="noadjust") {
    localcheck <- T; massadj <- F
  } else if (masspoints=="nolocalcheck") {
    localcheck <- F; massadj <- T
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

  # prepare params
  bins.p <- bins[1]; bins.s <- bins[2]
  tsha.p <- testshape[1]; tsha.s <- testshape[2]
  tmod.p <- testmodel[1]; tmod.s <- testmodel[2]
  nL <- length(testshapel); nR <- length(testshaper); nT <- length(testshape2)
  ntestshape <- nL + nR + nT

  if (binsmethod=="dpi") {
    selectmethod <- "IMSE direct plug-in"
  } else {
    selectmethod <- "IMSE rule-of-thumb"
  }

  knot <- NULL
  if (!is.character(binspos)) {
    nbins <- length(binspos)+1
    knot <- c(xmin, binspos, xmax)
    position <- "User-specified"
    es <- F
    selectmethod <- "User-specified"
  } else {
    if (binspos == "es") {
      es <- T
      position <- "Evenly-spaced"
    } else {
      es <- F
      position <- "Quantile-spaced"
    }
  }

  #### bin selection if needed ########
  if (is.null(nbins)) {
    # check if rot can be implemented
    if (is.null(nbinsrot)) {
      if (eN <= dfcheck[1]+bins.p+1+qrot) {
        print("too small effective sample size for bin selection.")
        stop()
      }
    }
    if (is.na(Ndist))  {
      Ndist.sel <- NULL
    } else {
      Ndist.sel <- Ndist
    }
    if (is.na(Nclust)) {
      Nclust.sel <- NULL
    } else {
      Nclust.sel <- Nclust
    }
    binselect <- binsregselect(y, x, w, deriv=deriv,
                              bins=bins, binspos=binspos,
                              binsmethod=binsmethod, nbinsrot=nbinsrot,
                              vce=vce, cluster=cluster,
                              dfcheck=dfcheck, masspoints=masspoints, weights=weights,
                              numdist=Ndist.sel, numclust=Nclust.sel)

    if (is.na(binselect$nbinsrot.regul)) {
      print("bin selection fails.")
      stop()
    }
    if (binsmethod == "rot") {
      nbins <- binselect$nbinsrot.regul
    } else if (binsmethod == "dpi") {
      nbins <- binselect$nbinsdpi
      if (is.na(nbins)) {
        warning("DPI selection fails. ROT choice used.")
        nbins <- binselect$nbinsrot.regul
      }
    }
  }

  # check eff size for testing
  tsha.fewobs <- F
  if (ntestshape!=0) if ((nbins-1)*(tsha.p-tsha.s+1)+tsha.p+1+dfcheck[2]>=eN) {
    tsha.fewobs <- T
    warning("too small eff. sample size for testing shape.")
  }

  tmod.fewobs <- F
  if (!is.null(testmodelparfit)|!is.null(testmodelpoly)) if ((nbins-1)*(tmod.p-tmod.s+1)+tmod.p+1+dfcheck[2]>=eN) {
    tmod.fewobs <- T
    warning("too small eff. sample size for testing models.")
  }

  # Generate knot
  if (is.null(knot)) {
    if (es) {
      knot <- genKnot.es(xmin, xmax, nbins)
    } else {
      knot <- genKnot.qs(x, nbins)
    }
  }
  # check local mass points
  knot <- c(knot[1], unique(knot[-1]))
  if (nbins!=length(knot)-1) {
    warning("repeated knots. Some bins dropped.")
    nbins <- length(knot)-1
  }
  if (localcheck) {
    uniqmin <- binsreg.checklocalmass(x, nbins, es, knot=knot) # mimic STATA
    if (ntestshape != 0) {
      if (uniqmin < tsha.p+1) {
        tsha.fewobs <- T
        warning("some bins have too few distinct values of x for testing shape.")
      }
    }
    if (!is.null(testmodelparfit)|!is.null(testmodelpoly)) {
      if (uniqmin < tmod.p+1) {
        tmod.fewobs <- T
        warning("some bins have too few distinct values of x for testing models.")
      }
    }
  }

  # seed
  if (!is.null(simsseed)) set.seed(simsseed)
  # prepare grid for uniform inference
  x.grid <- binsreg.grid(knot, simsgrid)$eval

  #####################################
  ####### Shape restriction test ######
  #####################################
  stat.shapeL <- pval.shapeL <- NA
  stat.shapeR <- pval.shapeR <- NA
  stat.shape2 <- pval.shape2 <- NA
  if (nL>0) {
    stat.shapeL <- matrix(NA, nL, 2)
    pval.shapeL <- c()
  }
  if (nR>0) {
    stat.shapeR <- matrix(NA, nR, 2)
    pval.shapeR <- c()
  }
  if (nT>0) {
    stat.shape2 <- matrix(NA, nT, 2)
    pval.shape2 <- c()
  }

  if (ntestshape != 0 & !tsha.fewobs) {
    B    <- binsreg.spdes(eval=x, p=tsha.p, s=tsha.s, knot=knot, deriv=0)
    k    <- ncol(B)
    P    <- cbind(B, w)
    model  <- lm(y ~ P-1, weights=weights)
    beta <- model$coeff[1:k]
    basis.sha <- binsreg.spdes(eval=x.grid, p=tsha.p, s=tsha.s, knot=knot, deriv=deriv)
    pred.sha  <- binsreg.pred(basis.sha, model, type = "all", vce=vce, cluster=cluster)
    fit.sha <- pred.sha$fit
    se.sha  <- pred.sha$se

    pos   <- !is.na(beta)
    k.new <- sum(pos)
    vcv.sha <- binsreg.vcov(model, type=vce, cluster=cluster)[1:k.new, 1:k.new]

    for (j in 1:ntestshape) {
      if (j <= nL) {
        stat.shapeL[j,2] <- 1
        stat.shapeL[j,1] <- max((fit.sha - testshapel[j]) / se.sha)
      } else if (j <= nL+nR) {
        stat.shapeR[j-nL,2] <- 2
        stat.shapeR[j-nL,1] <- min((fit.sha - testshaper[j-nL]) / se.sha)
      } else {
        stat.shape2[j-nL-nR,2] <- 3
        stat.shape2[j-nL-nR,1] <- max(abs((fit.sha - testshape2[j-nL-nR]) / se.sha))
      }
    }
    stat.shape <- NULL
    if (nL > 0) {
       stat.shape <- rbind(stat.shape, stat.shapeL)
    }
    if (nR > 0) {
       stat.shape <- rbind(stat.shape, stat.shapeR)
    }
    if (nT > 0) {
       stat.shape <- rbind(stat.shape, stat.shape2)
    }

    # Compute p-val
    Sigma.root <- lssqrtm(vcv.sha)
    num        <- basis.sha[,pos,drop=F] %*% Sigma.root
    pval.shape <- binsreg.pval(num, se.sha, nsims, tstat=stat.shape, side=NULL)$pval
    if (nL!=0) {
       stat.shapeL <- stat.shapeL[,1]
       pval.shapeL <- pval.shape[1:nL]
    }
    if (nR!=0) {
       stat.shapeR <- stat.shapeR[,1]
       pval.shapeR <- pval.shape[(nL+1):(nL+nR)]
    }
    if (nT!=0) {
       stat.shape2 <- stat.shape2[,1]
       pval.shape2 <- pval.shape[(nL+nR+1):ntestshape]
    }
  }

  ############################
  #### Specification Test ####
  ############################
  stat.poly <- pval.poly <- NA
  stat.mod <- pval.mod <- NA
  if (!is.null(testmodelpoly)) {
    stat.poly <- matrix(NA, 1,2)
    pval.poly <- c()
  }
  if (!is.null(testmodelparfit)) {
    stat.mod <- matrix(NA,ncol(testmodelparfit)-1,2)
    pval.mod <- c()
  }

  if ((!is.null(testmodelparfit)|!is.null(testmodelpoly))&!tmod.fewobs) {
    if (tmod.p==tsha.p & tmod.s==tsha.s & exists("vcv.sha")) {
      exist.mod <- T
      #beta.mod <- beta.sha
      vcv.mod  <- vcv.sha
      fit.mod  <- fit.sha
      se.mod   <- se.sha
      basis.mod <- basis.sha
      pos <- pos
      k.new <- k.new
      model <- model
    } else {
      exist.mod <- F
      B    <- binsreg.spdes(eval=x, p=tmod.p, s=tmod.s, knot=knot, deriv=0)
      k    <- ncol(B)
      P    <- cbind(B, w)
      model  <- lm(y ~ P-1, weights=weights)
      beta <- model$coeff[1:k]
      pos <- !is.na(beta)
      # beta.mod <- beta[pos]
      k.new <- sum(pos)
      vcv.mod <- binsreg.vcov(model, type=vce, cluster=cluster)[1:k.new, 1:k.new]
    }

    ######################
    # Test poly reg
    if (!is.null(testmodelpoly)) {
      if (!exist.mod) {
         basis.mod <- binsreg.spdes(eval=x.grid, p=tmod.p, s=tmod.s, knot=knot, deriv=deriv)
         pred.mod  <- binsreg.pred(basis.mod, model, type="all", vce=vce, cluster=cluster)
         fit.mod <- pred.mod$fit
         se.mod  <- pred.mod$se
      }

      # Run a poly reg
      x.p <- matrix(NA, N, testmodelpoly+1)
      for (j in 1:(testmodelpoly+1))  x.p[,j] <- x^(j-1)
      P.poly <- cbind(x.p, w)
      model.poly <- lm(y~P.poly-1, weights=weights)
      beta.poly <- model.poly$coefficients
      poly.fit <- 0
      for (j in deriv:testmodelpoly) {
          poly.fit <- poly.fit + x.grid^(j-deriv)*beta.poly[j+1]*factorial(j)/factorial(j-deriv)
      }

      stat.poly[1,2] <- 3
      stat.poly[1,1] <- max(abs((fit.mod - poly.fit) / se.mod))

      Sigma.root   <- lssqrtm(vcv.mod)
      num          <- basis.mod[,pos,drop=F] %*% Sigma.root
      pval.poly    <- binsreg.pval(num, se.mod, nsims, tstat=stat.poly, side=NULL)$pval
      stat.poly <- stat.poly[,1]
    }

    ##################################
    # Test model in another data frame
    if (!is.null(testmodelparfit)) {
       x.grid <- testmodelparfit[,1]
       basis.mod <- binsreg.spdes(eval=x.grid, p=tmod.p, s=tmod.s, knot=knot, deriv=deriv)
       pred.mod  <- binsreg.pred(basis.mod, model, type="all", vce=vce, cluster=cluster)
       fit.mod <- pred.mod$fit
       se.mod  <- pred.mod$se

       for (j in 2:ncol(testmodelparfit)) {
         stat.mod[j-1,2] <- 3
         stat.mod[j-1,1] <- max(abs((fit.mod - testmodelparfit[,j]) / se.mod))
       }
       Sigma.root   <- lssqrtm(vcv.mod)
       num          <- basis.mod[,pos,drop=F] %*% Sigma.root
       pval.mod  <- binsreg.pval(num, se.mod, nsims, tstat=stat.mod, side=NULL)$pval
       stat.mod <- stat.mod[,1]
    }
  }

  ######################
  #######output#########
  ######################
  out <- list(testshapeL=list(testvalL=testshapel, stat.shapeL=stat.shapeL, pval.shapeL=pval.shapeL),
              testshapeR=list(testvalR=testshaper, stat.shapeR=stat.shapeR, pval.shapeR=pval.shapeR),
              testshape2=list(testval2=testshape2, stat.shape2=stat.shape2, pval.shape2=pval.shape2),
              testpoly=list(testpoly=testmodelpoly, stat.poly=stat.poly, pval.poly=pval.poly),
              testmodel=list(stat.model=stat.mod, pval.model=pval.mod),
              opt = list(bins.p=bins.p, bins.s=bins.s, deriv=deriv,
                         testshape=testshape, testmodel=testmodel,
                         binspos=position, binsmethod=selectmethod,
                         n=N, Ndist=Ndist, Nclust=Nclust,
                         nbins=nbins, knot=knot))
  out$call   <- match.call()
  class(out) <- "CCFFbinsregtest"
  return(out)

}

##########################################################################
#' Internal function.
#'
#' @param x Class \code{CCFFbinsregtest} objects.
#'
#' @keywords internal
#' @export
#'

print.CCFFbinsregtest <- function(x, ...) {
  cat("Call: binsregtest\n\n")

  cat(paste("Sample size (n)                    =  ", x$opt$n,         "\n", sep=""))
  cat(paste("# of distinct values (Ndist)       =  ", x$opt$Ndist,     "\n", sep=""))
  cat(paste("# of clusters (Nclust)             =  ", x$opt$Nclust,    "\n", sep=""))
  cat(paste("Derivative (deriv)                 =  ", x$opt$deriv,     "\n", sep=""))
  cat(paste("Bin selection:", "\n"))
  cat(paste("  Method (binsmethod)              =  ", x$opt$binsmethod,"\n", sep=""))
  cat(paste("  Placement (binspos)              =  ", x$opt$binspos,   "\n", sep=""))
  cat(paste("  degree (p)                       =  ", x$opt$bins.p,     "\n", sep=""))
  cat(paste("  smooth (s)                       =  ", x$opt$bins.s,     "\n", sep=""))
  cat(paste("  # of bins (nbins)                =  ", x$opt$nbins,     "\n", sep=""))
  cat("\n")

}

################################################################################
#' Internal function.
#'
#' @param object Class \code{CCFFbinsregtest} objects.
#'
#' @keywords internal
#' @export
summary.CCFFbinsregtest <- function(object, ...) {
  x <- object
  args <- list(...)

  cat("Call: binsregtest\n\n")

  cat(paste("Sample size (n)                    =  ", x$opt$n,         "\n", sep=""))
  cat(paste("# of distinct values (Ndist)       =  ", x$opt$Ndist,     "\n", sep=""))
  cat(paste("# of clusters (Nclust)             =  ", x$opt$Nclust,    "\n", sep=""))
  cat(paste("Derivative (deriv)                 =  ", x$opt$deriv,     "\n", sep=""))
  cat(paste("Bin selection:", "\n"))
  cat(paste("  Method (binsmethod)              =  ", x$opt$binsmethod,"\n", sep=""))
  cat(paste("  Placement (binspos)              =  ", x$opt$binspos,   "\n", sep=""))
  cat(paste("  degree (p)                       =  ", x$opt$bins.p,     "\n", sep=""))
  cat(paste("  smooth (s)                       =  ", x$opt$bins.s,     "\n", sep=""))
  cat(paste("  # of bins (nbins)                =  ", x$opt$nbins,     "\n", sep=""))
  cat("\n")

  if (!is.null(x$testshapeL$testvalL)|!is.null(x$testshapeR$testvalR)|!is.null(x$testshape2$testval2)) {
    cat(paste("Shape Restriction Tests:")); cat("\n")
    cat(paste("degree (p) = ", x$opt$testshape[1], sep="")); cat(paste("   "))
    cat(paste("smooth (s) = ", x$opt$testshape[2], sep="")); cat("\n")

    if (!is.null(x$testshapeL$testvalL)) {
      cat(paste(rep("=", 15 + 12 + 12), collapse="")); cat("\n")
      cat(format("H0: sup mu <=", width = 15, justify = "left"))
      cat(format("sup T",         width = 12, justify = "right"))
      cat(format("p value",       width = 12, justify = "right"))
      cat("\n")
      cat(paste(rep("-", 15 + 12 + 12), collapse="")); cat("\n")
      cat("\n")
      for (i in 1:length(x$testshapeL$testvalL)) {
        cat(format(paste(x$testshapeL$testvalL[i]),    width = 15, justify = "centre"))
        cat(format(sprintf("%3.3f", x$testshapeL$stat.shapeL[i]), width = 12, justify = "right"))
        cat(format(sprintf("%3.3f", x$testshapeL$pval.shapeL[i]), width = 12, justify = "right"))
        cat("\n")
      }
      cat(paste(rep("-", 15 + 12 + 12), collapse=""))
      cat("\n")
      cat("\n")
    }

    if (!is.null(x$testshapeR$testvalR)) {
      cat(paste(rep("=", 15 + 12 + 12), collapse="")); cat("\n")
      cat(format("H0: inf mu >=", width = 15, justify = "left"))
      cat(format("inf T",         width = 12, justify = "right"))
      cat(format("p value",       width = 12, justify = "right"))
      cat("\n")
      cat(paste(rep("-", 15 + 12 + 12), collapse="")); cat("\n")
      cat("\n")
      for (i in 1:length(x$testshapeR$testvalR)) {
        cat(format(paste(x$testshapeR$testvalR[i]),               width = 15, justify = "centre"))
        cat(format(sprintf("%3.3f", x$testshapeR$stat.shapeR[i]), width = 12, justify = "right"))
        cat(format(sprintf("%3.3f", x$testshapeR$pval.shapeR[i]), width = 12, justify = "right"))
        cat("\n")
      }
      cat(paste(rep("-", 15 + 12 + 12), collapse=""))
      cat("\n")
      cat("\n")
    }

    if (!is.null(x$testshape2$testval2)) {
      cat(paste(rep("=", 15 + 12 + 12), collapse="")); cat("\n")
      cat(format("H0: mu =",      width = 15, justify = "left"))
      cat(format("sup |T|",       width = 12, justify = "right"))
      cat(format("p value",       width = 12, justify = "right"))
      cat("\n")
      cat(paste(rep("-", 15 + 12 + 12), collapse="")); cat("\n")
      cat("\n")
      for (i in 1:length(x$testshape2$testval2)) {
        cat(format(paste(x$testshape2$testval2[i]),               width = 15, justify = "centre"))
        cat(format(sprintf("%3.3f", x$testshape2$stat.shape2[i]), width = 12, justify = "right"))
        cat(format(sprintf("%3.3f", x$testshape2$pval.shape2[i]), width = 12, justify = "right"))
        cat("\n")
      }
      cat(paste(rep("-", 15 + 12 + 12), collapse=""))
      cat("\n")
      cat("\n")
    }
  }

  if (!is.null(x$testpoly$testpoly)|!all(is.na(x$testmodel$stat.model))) {
    cat(paste("Model Specification Tests:")); cat("\n")
    cat(paste("degree (p) = ", x$opt$testmodel[1], sep="")); cat(paste("   "))
    cat(paste("smooth (s) = ", x$opt$testmodel[2], sep="")); cat("\n")

    if (!is.null(x$testpoly$testpoly)) {
      cat(paste(rep("=", 15 + 12 + 12), collapse="")); cat("\n")
      cat(format("H0: mu =",      width = 15, justify = "left"))
      cat(format("sup |T|",       width = 12, justify = "right"))
      cat(format("p value",       width = 12, justify = "right"))
      cat("\n")
      cat(paste(rep("-", 15 + 12 + 12), collapse="")); cat("\n")
      cat("\n")
      cat(format(paste("poly. p=", x$testpoly$testpoly, sep=""),  width = 15, justify = "left"))
      cat(format(sprintf("%3.3f", x$testpoly$stat.poly), width = 12, justify = "right"))
      cat(format(sprintf("%3.3f", x$testpoly$pval.poly), width = 12, justify = "right"))
      cat("\n")
      cat(paste(rep("-", 15 + 12 + 12), collapse=""))
      cat("\n")
    }

    if (!all(is.na(x$testmodel$stat.model))) {
      cat(paste(rep("=", 15 + 12 + 12), collapse="")); cat("\n")
      cat(format("H0: mu =",      width = 15, justify = "left"))
      cat(format("sup |T|",       width = 12, justify = "right"))
      cat(format("p value",       width = 12, justify = "right"))
      cat("\n")
      cat(paste(rep("-", 15 + 12 + 12), collapse="")); cat("\n")
      cat("\n")
      for (i in 1:length(x$testmodel$stat.model)) {
        cat(format(paste("model ", i, sep=""),   width = 15, justify = "centre"))
        cat(format(sprintf("%3.3f", x$testmodel$stat.model[i]),   width = 12, justify = "right"))
        cat(format(sprintf("%3.3f", x$testmodel$pval.model[i]),   width = 12, justify = "right"))
        cat("\n")
      }
      cat(paste(rep("-", 15 + 12 + 12), collapse=""))
      cat("\n")
    }
  }

}
