#' Separable Component Decomposition of the Empirical Covariance
#'
#' Performs the power iteration method to find the separable expansion for
#' \code{C = cov(X)}, i.e. to find the decomposition \code{C = sigma1 A1} x
#' \code{B1} + \code{sigma2 A2} x \code{B2} + ... + \code{sigmaR AR} x
#' \code{BR}. When the degree-of-separability \code{R} is not provided,
#' cross-validation (CV) is performed.
#'
#' @param X data array of size \code{N} x \code{K1} x \code{K2}
#' @param R integer, the degree-of-separability
#' @param B array of size \code{R} x \code{K2} x \code{K2} specifying the
#'          initial values for \code{B1},...,\code{BR}. Defaults to partially
#'          traced initializations. If \code{R} not provided by the user,
#'          \code{B} should be array of size \code{K2} x \code{K2}, and this
#'          single starting point is used for every degree-of-separability.
#' @param mu matrix of size \code{K1} x \code{K2}, the mean, defaults to the
#'           empirical mean
#' @param predict flag whether prediction-based cross-validation should be used,
#'                the default is fit-based cross-validation
#' @param maxiter integer, maximum number of power iterations, defaults to 10
#' @param maxR integer, maximum value of the degree-of-separability \code{R}
#'             to be probed via CV (used only when \code{R} not provided
#'             directly)
#' @param Folds number of folds for cross-validation, defaults to 10
#' @param perc vector of 2, if \code{predict=T}, what proportions (in time and
#'             space) of every observation should be predicted. When
#'             \code{perc}>=1, specific number of grid lines is predicted, e.g.
#'             \code{perc}=1 corresponds to one-step-ahead prediction (in both
#'             space and time)
#'
#' @return list of the following objects:list of 3:
#' * \code{A} - array of size \code{R} x \code{K1} x \code{K1}, the first
#'   \code{R} temporal kernels;
#' * \code{B} - array of size \code{R} x \code{K2} x \code{K2}, the first
#'   \code{R} spatial kernels;
#' * \code{sigma} - vector of length \code{R}, the separable component scores
#'
#' if cross-validation was conducted, the list contains also:
#'
#' * \code{cv} - the CV objective values, a smaller value is better
#' * \code{R} - the chosen bandwidth
#'
#' @export
#'
#' @examples
#' N <- 10
#' K1 <- K2 <- 10
#' x <- 1:K1/(K1+1)
#' e0 <- rep(1,K1)
#' e1 <- sin(2*pi*x)
#' e2 <- cos(2*pi*x)
#' X <- array(0,c(N,K1,K2))
#' set.seed(123)
#' for(n in 1:N){
#'   X[n,,] <- 2*rnorm(1)*outer(e0,e0) + rnorm(1)*outer(e1,e1) +
#'     rnorm(1)*outer(e2,e2)
#' }
#' SCD <- scd(X)
#' SCD$R # the degree-of-separability chosen by fit-based CV
#' SCD <- scd(X,predict=TRUE)
#' SCD$R # the degree-of-separability leading to best predictions
#' SCD <- scd(X,predict=TRUE,perc=1)
#' SCD$R # the degree-of-separability leading to best one-step-ahead (both in
#'       # space and time) predictions
scd <- function(X, R=NULL, B=NULL, mu=NULL, predict=FALSE, maxiter=10,
                maxR=NULL, Folds=NULL, perc=c(1/4,1/4)){
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
    warning(paste0("Not enough variability in some locations."))
  }
  N <- dim(X)[1]
  Kmin <-  min(dim(X)[-1])
  # degree-of-separability R
  if(length(R) > 0){
    R <- as.integer(R)
    if(!is.integer(R) || R < 1){
      stop(paste0("The degree-of-separability R has to be a positive integer."))
    }
    if(R > N || R > Kmin)
    {
      stop(paste0("The chose degree-of-separability R is too high.
                     It needs to be R <= min(dim(X))."))
    }
  }
  # maximal degree for cross-validation maxR
  if(length(maxR)>0){
    maxR <- as.integer(maxR)
    if(!is.integer(maxR) || maxR < 1){
      stop(paste0("Them maximal degree-of-separability maxR has to be
                  a positive integer."))
    }
    if(maxR > N || maxR > Kmin)
    {
      stop(paste0("The chosen maximal degree-of-separability maxR is too high.
                  It needs to be maxR <= min(dim(X))"))
    }
  }
  # starting point - below, only when R is provided and no CV utilized
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
  # maximum no. of iterations maxiter
  maxiter <- as.integer(maxiter)
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
  #
  if(length(maxR)==0) maxR <- 3
  if(length(Folds)==0) Folds <- min(10, floor(N/2))
  if(N < 2*Folds || Folds==1){
    stop(paste0("Not enough surfaces to do cross-validation, choose the
                bandwidth manually."))
  }
  ### cross-validation (CV)
  CV <- FALSE
  X <- sweep(X, c(2,3), mu)
  if(length(R)==0){
    CV <- TRUE
    if(predict){
      cvscores <- cvR_pred(X,Folds,maxR,maxiter,perc)
    }else{
      cvscores <- cvR(X,Folds,maxR,maxiter)
    }
    R <- min(localMaxima(-cvscores))
  }
  # check out the starting point B, if provided by the user
  if(length(B) > 0){
    if(length(dim(B)) != 3){
      if(length(dim(B)) == 2){
        Bb <- B
        B <- array(0,c(R,dim(B)))
        for(r in 1:R) B[r,,] <- Bb
      }else{
        warning(paste0("Cannot use the provided starting point B."))
        B <- NULL
      }
    }
    if(dim(B)[2] != dim(X)[3] || dim(B)[3] != dim(X)[3]){
      warning(paste0("Dimensions of the starting point B do not correspond to
                     the spatial covariance kernel of X. Discarding
                     the provided B."))
      B <- NULL
    }
  }
  ### Estimation
  if(CV){
    Res <- scd_est(X,R,maxiter)
    return(list(A=Res$A, B=Res$B, sigma=Res$sigma,
                cv=cvscores, R=R))
  }else return(scd_est(X,R,maxiter,B))
}
