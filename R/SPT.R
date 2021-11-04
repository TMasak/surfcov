#' Shifted Partial Tracing
#'
#' Estimates the separable part of the separable-plus-banded model via Shifted
#' Partial Tracing. Choosing \code{d}=0 defaults to standard (non-shifted) partial
#' tracing estimator of the separable model.
#'
#' @param X data array of size \code{N} x \code{K1} x \code{K2}
#' @param d bandwidth, integer \code{>= 0}
#' @param standardize whether the results should be standardized,
#'                    defaults to True, option False is only for comparison
#'                    purposes with the partial inner product
#' @return list of two:
#'        \code{A1} - matrix of size \code{K1} x \code{K1}, the temporal kernel;
#'         \code{A2} - matrix of size \code{K2} x \code{K2}, the spatial kernel;
#' @export
#'
#' @examples
#' X <- array(runif(20*3*4),c(20,3,4))
#' spt(X,0)
spt <- function(X, d=0, standardize=TRUE){
  N <- dim(X)[1]
  K1 <- dim(X)[2]
  K2 <- dim(X)[3]
  W1 <- matrix(0,K1,K1) # weighting matrix, will have ones on the `d'-off-
  indx <- 1:(K1-d)  # diagonals (both of them, making the results symmetric)
  W1[cbind(indx+d,indx)] <- W1[cbind(indx,indx+d)] <- 1
  W2 <- matrix(0,K2,K2)
  indx <- 1:(K2-d)
  W2[cbind(indx+d,indx)] <- W2[cbind(indx,indx+d)] <- 1

  A1 <- T1(X, NULL, ext2arr(W2), NULL, 1, N, K1, K2)
  A2 <- T2(X, ext2arr(W1), NULL, NULL, 1, N, K1, K2)

  A1 <- A1*sign(sum(diag(A1)))
  A2 <- A2*sign(sum(diag(A2))) # flip signs if shifted traces negative
  if(standardize){
    A1 <- A1/sqrt(abs(sum(A1*W1)))
    A2 <- A2/sqrt(abs(sum(A2*W2)))
  }
  return(list(A1 = A1, A2 = A2))
}

#' Estimate the Entry-wise Variance
#'
#' Estimates the variance of heteroscedastic white noise after the separable
#' part of the model has already been estimated.
#'
#' @param X data array of size \code{N} x \code{K1} x \code{K2}
#' @param A1 temporal kernel of the separable part
#' @param A2 spatial kernel of the separable part
#'
#' @return matrix of size K1 x K2 providing the entry-wise variance
#' @export
#'
#' @examples
#' X <- array(runif(20*3*4),c(20,3,4))
#' A1 <- matrix(runif(3^2),3)
#' A2 <- matrix(runif(4^2),4)
#' diagonal_band(X,A1,A2)
diagonal_band <- function(X,A1,A2){
  K1 <- dim(A1)[1]
  K2 <- dim(A2)[1]
  B_var <- array(0,c(K1,K2))
  for(k1 in 1:K1){
    for(k2 in 1:K2){
      B_var[k1,k2] <- stats::var(X[,k1,k2]) - A1[k1,k1]*A2[k2,k2]
    }
  }
  return(B_var)
}

#' Toeplitz Averaging
#'
#' @param X data array of size \code{N} x \code{K1} x \code{K2}
#' @param A1,A2 the separable part of the covariance
#' @param d bandwidth
#'
#' @return the stationary band - matrix of size \code{d} x \code{d}
#' @export
#'
#' @examples
#' X <- array(rnorm(5*10*7),c(5,10,7))
#' Ta(X,d=2)
Ta <- function(X, A1=NULL, A2=NULL, d){
  N <- dim(X)[1]
  K1 <- dim(X)[2]
  K2 <- dim(X)[3]
  band <- array(0, dim = c(K1, K2))
  for(n in 1:N){
    IFTFT <- stats::fft(abs(stats::fft(X[n,,]))^2, inverse = TRUE)/K1/K2
    band <- band + Re(IFTFT)
  }
  for(t in 0:(d-1)){
    for(s in 0:(d-1)){
      band[t+1,s+1] <- band[t+1,s+1]/K1/K2
    }
  }
  if(is.null(A1)){
    return(band[1:d,1:d]/N)
  }else{
    band2 <- array(0, dim = c(d, d))
    for(t in 0:(d-1)){
      ind <- 1:(K1-t)
      alpha <- sum(A1[cbind(ind+t,ind)])
      for(s in 0:(d-1)){
        ind <- 1:(K2-s)
        beta <- sum(A2[cbind(ind+s,ind)])
        band2[t+1,s+1] <- alpha*beta/K1/K2
      }
    }
    return(band[1:d,1:d]/N - band2)
  }
}
