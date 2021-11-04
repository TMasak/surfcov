#' Goodness-of-fit testing for separable-plus-banded model based on empirical
#' bootstrap.
#'
#' @param X data array of size \code{N} x \code{K1} x \code{K2}
#' @param d the bandwidth, an integer between \code{0} and \code{min(K1,K2)}
#' @param stationary a flag whether the banded part should be stationary,
#'                   defaults to \code{TRUE}
#' @param I number of temporal eigenfunctions used to assess the fit
#' @param J number of spatial eigenfunctions used to assess the fit
#' @param noB number of bootstrap samples, defaults to 1000
#' @param mu matrix of size \code{K1} x \code{K2}, the mean, defaults to the
#'           empirical mean
#' @param verbose flag whether the progress should be printed
#'
#' @return the p-value (a number - if small, it suggests the model is not valid)
#' @export
#'
#' @examples
#' set.seed(123)
#' N <- 100
#' K1 <- K2 <- 7
#' X <- array(rnorm(N*K1*K2),c(N,K1,K2))
#' A1 <- array(0,c(K1,K1))
#' A1[,1] <- 1
#' A1[,2] <- -3:3
#' A1 <- A1 %*% t(A1)
#' A2 <- A1 <- mat_root(A1)
#' for(n in 1:N){
#'   X[n,,] <- A1 %*% X[n,,] %*% A2 + 2*matrix(rnorm(K1*K2),K1)
#' }
#' # the truth is separable+stationary (with d=1), i.e. it is also
#' # separable+diagonal at the same time
#' # set e.g. noB=1e4 in code below to get
#' test_spb(X,d=0,stationary=TRUE,noB=10)
#' # with d=0 the model is separable, which is rejected here
#' test_spb(X,d=1,stationary=TRUE,noB=10)
#' # separable+stationary model is not rejected
#' test_spb(X,d=1,stationary=FALSE,noB=10)
#' # separable+diagonal model is also not rejected
test_spb <- function(X,d=1,stationary=TRUE,I=1,J=1,noB=1e3,mu=NULL,
                     verbose=TRUE){
  ### checking inputs
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
  if(length(d) > 0){
    if(d > Kmin){
      stop(paste0("Provided bandwidth is too large."))
    }
    if(d < 1) d <- floor(d*Kmin)
  }
  if(length(mu) > 0){
    if(dim(mu)[1] != dim(X)[2] || dim(mu)[2] != dim(X)[3]){
      warning(paste0("Dimensions of the mean mu do not correspond to the
                     data size. Discarding the provided mu."))
      mu <- NULL
    }
  }
  if(d != 1 && stationary==FALSE){
    warning(paste0("Non-stationary banded part currently available only for
                   d=1. Setting d=1."))
    d <- 1
  }
  ### bootstrap test
  T_U <- dis2sep(X,d,mu,I,J,stationary)
  stat <- sum(T_U^2)
  stat_boot <- rep(0,noB)
  for(b in 1:noB){
    cat(paste0("Bootstrapped sample no.: ",b,"/",noB),"\r")
    Xb <- X[sample.int(N, N, replace = TRUE),,]
    T_Ub <- dis2sep(Xb,d,mu,I,J,stationary)
    stat_boot[b] <- sum((T_U - T_Ub)^2)
  }
  cat("\n")
  return(mean(stat_boot > stat))
}

#' Projected distances to separability, i.e. the bootstrapped statistic.
#'
#' @param X data array of size \code{N} x \code{K1} x \code{K2}
#' @param d the bandwidth, an integer between \code{0} and \code{min(K1,K2)}
#' @param mu matrix of size \code{K1} x \code{K2}, the mean, defaults to the
#'           empirical mean
#' @param I number of temporal eigenfunctions used to assess the fit
#' @param J number of spatial eigenfunctions used to assess the fit
#' @param stationary a flag whether the banded part should be stationary,
#'                   defaults to \code{TRUE}
#'
#' @return array of projected distances to separability
#' @export
#'
#' @examples
dis2sep <- function(X,d,mu,I,J,stationary){
  N <- dim(X)[1]
  if(length(mu)==0) mu <- apply(X, c(2,3), mean)
  X <- sweep(X, c(2,3), mu)
  Fit <- spb(X,d=d,stationary=stationary)
  if(sum(!is.finite(Fit$A1)) + sum(!is.finite(Fit$A2)) +
     sum(!is.finite(Fit$B)) > 1){
    stop(paste0("Not enough variability in a bootstrap sample. Potential reason:
                a part of the dataset is only a single repeated observation."))
  }

  SVD <- svd(Fit$A1)
  Es <- SVD$u[,1:I,drop=F]
  lambda <- SVD$d[1:I]
  SVD <- svd(Fit$A2)
  Fs <- SVD$u[,1:J,drop=F]
  gamma <- SVD$d[1:J]

  T_ij <- array(0,c(I,J))
  for(i in 1:I){
    for(j in 1:J){
      pom <- 0
      for(n in 1:N){
        pom <- pom + ( t(Es[,i]) %*% X[n,,] %*% Fs[,j])^2
      }
      pom <- pom/N
      T_ij[i,j] <- pom - lambda[i]*gamma[j]
      pom <- Es[,i] %*% t(Fs[,j])
      if(!stationary){
        # reduction of the test statistic by what variability is explained by
        # the banded part in case of separable+noise model (i.e. diagonal B)
        T_ij[i,j] <- T_ij[i,j] - sum(Fit$B*pom^2)
      }else{
        # same as above in case of separable+stationary model
        eigvals <- stats::fft(to_book_format(Fit$B,dim(X)[2],dim(X)[3]))
        RHS <- BXfast(Re(eigvals), pom)
        T_ij[i,j] <- T_ij[i,j] - sum(pom*RHS) # reduction by banded part
      }
    }
  }
  return(sqrt(N) * T_ij)
}
