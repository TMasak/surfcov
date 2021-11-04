#' Estimate the Separable Component Decomposition
#'
#' Calculates the R-separable estimator of the covariance of X via the
#' generalized power iteration, starting from B.
#'
#' @param X data set, array of size \code{N} x \code{K1} x \code{K2}
#' @param R integer, the degree of separability
#' @param maxiter maximum number of iteration
#' @param B a spatial kernel used as a starting point
#'
#' @return list of 3:
#'         \code{A} - array of size \code{R} x \code{K1} x \code{K1},
#'                    the first \code{R} temporal kernels;
#'         \code{B} - array of size \code{R} x \code{K2} x \code{K2},
#'                    the first \code{R} spatial kernels;
#'         \code{sigma} - vector of length \code{R},
#'                        the separable component scores
#' @export
#'
#' @examples
#' X <- array(runif(20*3*4),c(20,3,4))
#' scd_est(X,R=2)
scd_est <- function(X,R=1,maxiter=10,B=NULL){
  N <- dim(X)[1]
  K1 <- dim(X)[2]
  K2 <- dim(X)[3]
  if(length(B) == 0){
    B <- array(0, c(R,K2,K2))
    for(r in 1:R){
      B[r,,] <- diag(K2)
    }
  }
  A <- array(0,c(R,K1,K1))
  sigma <- rep(0,R)
  for(r in 1:R){
    iter <- 0
    while(iter < maxiter){
      A[r,,] <- T1(X,A,B,sigma,r,N,K1,K2)
      A[r,,] <- A[r,,]/norm(A[r,,],type="F")
      B[r,,] <- T2(X,A,B,sigma,r,N,K1,K2)
      sigma[r] <- norm(B[r,,],type="F")
      B[r,,] <- B[r,,]/sigma[r]
      iter <- iter +1
    }
  }
  return(list(A=A,B=B,sigma=sigma))
}

#' Partial Inner Product w.r.t. the First Argument
#'
#' Calculates T_1(X x X, B_r), i.e. the partial inner product of the
#' covariance of X with the current weight B_r
#' resulting in a temporal kernel (proxy for A_r), which is then
#' orthogonalized w.r.t A_1,...,A_{r-1}.
#'
#' @param X data set, array of size \code{N} x \code{K1} x \code{K2}
#' @param A temporal kernels, array of size (\code{>= r}) x \code{K1} x
#'          \code{K1}
#' @param B spatial kernels, array of size (\code{>= r}) x \code{K2} x \code{K2}
#' @param sigma separable component scores, vector of size \code{>= r}
#' @param r current sought separable rank
#' @param N sample size
#' @param K1 temporal grid size
#' @param K2 spatial grid size
#'
#' @return \code{r}-th temporal kernel A_r, array of size \code{K1} x \code{K1}
#' @export
#'
#' @examples
#' X <- array(rnorm(20*3*4),c(20,3,4))
#' A <- array(runif(2*3^2),c(2,3,3))
#' B <- array(runif(2*4^2),c(2,4,4))
#' T1(X, A, B, c(1,1), 2, 20, 3, 4)
T1 <- function (X,A,B,sigma,r,N,K1,K2){
  Res <- array(0,c(K1,K1))
  for(n in 1:N){
    Res <- Res + X[n,,] %*% B[r,,] %*% t(X[n,,])
  }
  Res <- Res/N
  if(r > 1){
    for(j in 1:(r-1)){
      Res <- Res - sigma[j] * sum(B[j,,]*B[r,,]) * A[j,,]
    }
  }
  return(Res)
}

#' Partial Inner Product w.r.t. the Second Argument
#'
#' Calculates T_2(X x X, A_r), i.e. the partial inner product of the
#' covariance of X with the current weight A_r
#' resulting in a temporal kernel (proxy for B_r), which is then
#' orthogonalized w.r.t B_1,...,B_{r-1}.
#'
#' @param X data set, array of size \code{N} x \code{K1} x \code{K2}
#' @param A temporal kernels, array of size (\code{>= r}) x \code{K1} x
#'          \code{K1}
#' @param B spatial kernels, array of size (\code{>= r}) x \code{K2} x \code{K2}
#' @param sigma separable component scores, vector of size \code{>= r}
#' @param r current sought separable rank
#' @param N sample size
#' @param K1 temporal grid size
#' @param K2 spatial grid size
#'
#' @return \code{r}-th spatial kernel B_r, array of size \code{K1} x \code{K1}
#' @export
#'
#' @examples
#' X <- array(rnorm(20*3*4),c(20,3,4))
#' A <- array(runif(2*3^2),c(2,3,3))
#' B <- array(runif(2*4^2),c(2,4,4))
#' T2(X, A, B, c(1,1), 2, 20, 3, 4)
T2 <- function (X,A,B,sigma,r,N,K1,K2){
  Res <- array(0,c(K2,K2))
  for(n in 1:N){
    Res <- Res + t(X[n,,]) %*% A[r,,] %*% X[n,,]
  }
  Res <- Res/N
  if(r>1){
    for(j in 1:(r-1)){
      Res <- Res - sigma[j] * sum(A[j,,]*A[r,,]) * B[j,,]
    }
  }
  return(Res)
}


