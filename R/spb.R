#' Separable-plus-banded Covariance Estimation
#'
#' Estimates components of \code{cov(X) = A1 x A2 + B} where \code{B} is banded
#' by \code{d}. If bandwidth \code{d} is not provided, it is estimated via
#' cross-validation. \cr
#' CV can be either fit-based, i.e tailored to minimize mean squared estimation
#' error, or prediction-based, i.e. tailored to minimize the mean squared
#' prediction error. The fit-based CV is faster. \cr
#' \code{B} is considered stationary by default, this can be
#' changed by toggling \code{stationary=F}.
#'
#' @param X data array of size \code{N} x \code{K1} x \code{K2}
#' @param d the bandwidth, an integer between \code{0} and \code{min(K1,K2)}
#' @param stationary a flag whether the banded part should be stationary,
#'                   defaults to \code{TRUE}
#' @param mu matrix of size \code{K1} x \code{K2}, the mean, defaults to the
#'           empirical mean
#' @param predict flag whether prediction-based cross-validation should be used,
#'                the default is fit-based cross-validation
#' @param maxd integer, maximum bandwidth considered by CV, defaults to
#'             \code{min(K1,K2)/10}
#' @param mind minimum bandwidth considered by CV, defaults to zero
#' @param Folds number of folds for cross-validation, defaults to 10
#' @param perc vector of 2, if \code{predict=T}, what proportions (in time and
#'             space) of every observation should be predicted. When
#'             \code{perc}>=1, specific number of grid lines is predicted, e.g.
#'             \code{perc}=1 corresponds to one-step-ahead prediction (in both
#'             space and time)
#'
#' @return a list of objects, depending on \code{stationary} and \code{d}: \cr
#' * \code{A1} - matrix of size \code{K1} x \code{K1}, the temporal kernel
#' * \code{A2} - matrix of size \code{K2} x \code{K2}, the spatial kernel
#' * \code{B} - matrix of size \code{d} x \code{d}, i.e. the symbol of
#'   the banded part, if \code{stationary=T}. If \code{stationary=F}, then
#'   the value depends on \code{d}: for \code{d=1} it is a matrix of size
#'   \code{K1} x \code{K2} representing location-wise variances, while for
#'   \code{d>1} a structure produced by [banded_nonstat_est] is returned
#'
#' if cross-validation was conducted, the list contains also:
#'
#' * \code{cv} - the CV objective values: a larger value is better for fit-based
#'   CV and a smaller value is better for the prediction-based CV (when
#'   \code{predict=T})
#' * \code{d} - the chosen bandwidth
#' @export
#'
#' @examples
#' N <- 100
#' K1 <- K2 <- 7
#' X <- array(rnorm(N*K1*K2),c(N,K1,K2))
#' A1 <- array(0,c(K1,K1))
#' A1[,1] <- 1
#' A1[,2] <- -3:3
#' A1 <- A1 %*% t(A1)
#' A2 <- A1 <- mat_root(A1)
#' for(n in 1:N){
#'   X[n,,] <- A1 %*% X[n,,] %*% A2 + matrix(rnorm(K1*K2),K1)
#' }
#' # we have separable-plus-banded model, where B is stationary and d=1
#' Chat <- spb(X)
#' Chat$d # estimated bandwidth by fit-based CV
#' Chat <- spb(X,stationary=FALSE) # when we do not use stationarity
#' Chat$d # estimated bandwidth by fit-based CV without stationarity
spb <- function(X, d=NULL, stationary=TRUE, mu=NULL, predict=FALSE,
                maxd=NULL, mind=NULL, Folds=NULL, perc=c(1/4,1/4)){
  ### checking inputs
  # data set X
  if(length(dim(X)) != 3){
    stop(paste0("The data set has to form a 3-dimensional array."))
  }
  if(sum(I(dim(X) > 1)) < 3){
    stop(paste0("The data set has to form a 3-dimensional array with all
                dimensions being non-trivial."))
  }
  if(sum(!is.finite(X)) > 0){
    stop(paste0("Infinite or missing values in the data set."))
  }
  if(sum(!is.numeric(X)) > 0){
    stop(paste0("Non-numeric values in the data set."))
  }
  if(sum(I(apply(X, c(2,3), stats::var) < 10*.Machine$double.eps)) > 0){
    warning(paste0("Not enough variability in the data."))
  }
  N <- dim(X)[1]
  Kmin <-  min(dim(X)[-1])
  # bandwidth d
  if(length(d) > 0){
    if(d > Kmin){
      stop(paste0("Provided bandwidth is too large."))
    }
    if(d < 1) d <- floor(d*Kmin)
  }
  # range of bandwidths mind, maxd for CV
  if(length(mind) + length(maxd) == 1){
    stop(paste0("Please provide either both mind and max, or neither."))
  }
  if(length(mind) + length(maxd) == 2){
    if(mind > Kmin){
      stop(paste0("Minimum bandwidth is too large."))
    }
    if(mind < 1) mind <- floor(mind*Kmin)
    if(maxd > Kmin) maxd <- Kmin - 1
    if(maxd < 1) maxd <- floor(maxd*Kmin)
    if(maxd <= mind){
      d <- min(mind,maxd)
      warning(paste0("Since maxd <= mind, taking d=min(mind,maxd)."))
    }
  }
  # mean mu
  if(length(mu) > 0){
    if(dim(mu)[1] != dim(X)[2] || dim(mu)[2] != dim(X)[3]){
      warning(paste0("Dimensions of the mean mu do not correspond to the
                     data size. Discarding the provided mu."))
      mu <- NULL
    }
  }
  # prediction flag predict
  predict <- as.logical(predict)
  # number of splits for the cross-validation Folds
  if(length(Folds) > 0){
    Folds <- as.integer(Folds)
    Folds <- min(max(2,Folds),floor(N/2))
  }
  # percentages for prediction
  if(length(perc) > 1) perc <- perc[1:2] else perc <- c(perc,perc)
  if(perc[1] < 0 || perc[2] < 0){
    warning(paste0("Percentages for prediction need to be positive,
                   collapsing to the default."))
    perc <- c(1/4,1/4)
  }
  if(perc[1] >= Kmin || perc[2] >= Kmin){
    stop(paste0("Cannot understand the provided perc argument."))
  }
  if(perc[1] >= 1) perc[1] <- perc[1]/dim(X)[2]
  if(perc[2] >= 1) perc[2] <- perc[2]/dim(X)[3]
  ### set unused variables to their default
  if(length(mu)==0){ # use empirical mean to center the data
    mu <- apply(X, c(2,3), mean)
  }
  if(length(mind)==0) mind <- 0
  if(length(Folds)==0) Folds <- min(10, floor(N/2))
  if(length(d) == 0){
    if(N < 2*Folds || Folds==1){
      stop(paste0("Not enough surfaces to do cross-validation, choose the
                  bandwidth manually."))
    }
  }
  # cross-validation (CV)
  X <- sweep(X, c(2,3), mu)
  CV <- FALSE
  if(length(d)==0){
    CV <- TRUE
    if(predict){ # prediction-based CV
      if(stationary){
        if(length(maxd)==0) maxd <- min(Kmin-1,max(2,floor(Kmin/10)))
        cvscores <- cvband_pred(X, Folds, maxd, mind, perc)
        names(cvscores) <- mind:maxd
        d <- min(localMaxima(-cvscores)) + mind - 1
      }else{
        stop(paste0("Prediction-based cross-validation only available when
                    `stationary=TRUE'."))
      }
    }else{ # fit-based CV
      if(stationary){
        if(length(maxd)==0) maxd <- min(Kmin-1,max(2,floor(Kmin/10)))
        cvscores <- cvband(X, Folds, maxd, mind)
        cvscores <- cvscores[1,]/cvscores[2,]
        names(cvscores) <- mind:maxd
        d <- min(localMaxima(cvscores)) + mind - 1
      }else{
        if(length(maxd)==0) maxd <- min(2,Kmin-1)
        cvscores <- cvband_nonstat(X, Folds, maxd, mind)
        cvscores <- cvscores[1,]/cvscores[2,]
        names(cvscores) <- mind:maxd
        d <- min(localMaxima(cvscores)) + mind - 1
      }
    }
  }
  ### estimation
  Est <- spt(X, d)
  B <- as.matrix(0)
  if(d > 0){
    if(stationary){
      B <- as.matrix(Ta(X, Est$A1, Est$A2, d))
    }else{
      if(d == 1){
        B <- diagonal_band(X, Est$A1, Est$A2)
      }else{
        B <- banded_nonstat_est(X, d, Est$A1, Est$A2)
      }
    }
  }
  if(CV){
    return(list(mu=mu, A1=Est$A1, A2=Est$A2, B=B, cv=cvscores, d=d))
  }else{
    return(list(mu=mu, A1=Est$A1, A2=Est$A2, B=B))
  }
}
