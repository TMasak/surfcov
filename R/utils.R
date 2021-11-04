#' Find Local Maxima of a Vector
#'
#' @param x a vector
#'
#' @return maximum
#' @export
#'
#' @examples
#' localMaxima(runif(9))
localMaxima <- function(x) {
  y <- diff(c(-.Machine$integer.max, x)) > 0L
  rle(y)$lengths
  y <- cumsum(rle(y)$lengths)
  y <- y[seq.int(1L, length(y), 2L)]
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  return(y)
}

#' Epanechnikov Kernel
#'
#' Returns f(x) where f is the Epanechnikov kernel supported on (-1,1).
#'
#' @param x double
#'
#' @return value of epanechnikov kernel at \code{x}
#' @export
#'
#' @examples epanechnik(0)
epanechnik <- function(x){
  y <- 3/4*(1-x^2)*I(-1 <= x)*I(x <= 1)
  return(y)
}

#' Covariance of the Wiener Process
#'
#' @param K grid size
#'
#' @return the covariance, i.e. matrix of size \code{K} x \code{K}
#' @export
#'
#' @examples
#' brownian_cov(20)
brownian_cov <- function(K){
  B <- array(0,c(K,K))
  for(i in 1:K){
    for(j in 1:K){
      B[i,j] <- min(i,j)
    }
  }
  return(B)
}

#' Extend to Array
#'
#' Extends a matrix M into array of size 1 x \code{dim(M)[1]} x
#' \code{dim(M)[2]}.
#'
#' @param M a numeric matrix
#'
#' @return array of size \code{1} x \code{dim(M)[1]} x \code{dim(M)[2]}
#' @export
#'
#' @examples ext2arr(matrix(rnorm(9),3))
ext2arr <- function(M){
  A <- array(0,c(1,dim(M)))
  A[1,,] <- M
  return(A)
}

#' Transform the symbol of B to format suitable for 2D FFT
#'
#' @param band the symbol of B, matrix of size \code{d} x \code{d}, where
#'   \code{d} is an integer
#' @param K1 temporal grid size, integer \code{>= 2}
#' @param K2 spatial grid size, integer \code{>= 2}
#'
#' @return matrix of size \code{2*K1} x \code{2*K2}
#' @export
#'
#' @examples to_book_format(matrix(runif(4),2), 3, 5)
to_book_format <- function(band, K1, K2){
  Bnew <- array(0, c(2*K1,2*K2))
  delta <- min(dim(band))
  Bnew[1:delta,1:delta] <- band[1:delta,1:delta]
  if(delta > 1){
    Bnew[(2*K1):(K1+2),1:K2] <- Bnew[2:K1,1:K2]
    Bnew[1:K1,(2*K2):(K2+2)] <- Bnew[1:K1,2:K2]
    Bnew[(2*K1):(K1+2),(2*K2):(K2+2)] <- Bnew[2:K1,2:K2]
  }
  return(Bnew)
}

#' Fast Multiplication for a Two-level Toeplitz Matrix
#'
#' We wish to calculate the tensor-matrix product Y = B X fast, where B is a
#' stationary tensor, i.e. vec(Y) = mat(B) vec(X) is a product of a Two-level
#' Toeplitz matrix with a vector, which can be calculated fast using 2D FFT.
#'
#' @param eigvals eigenvalues of a Two-level Toeplitz Matrix
#' @param x vectorization of matrix observation
#'
#' @return vector vec(Y), which can be matrized to obtain the product Y
#' @export
#'
#' @examples
#' band <- matrix(4:1,2)
#' eigvals <- stats::fft(to_book_format(band,3,3))
#' BXfast(Re(eigvals), runif(9))
BXfast <- function(eigvals,x){
  n2 <- dim(eigvals)[1]
  m2 <-dim(eigvals)[2]
  n <- n2/2
  m <- m2/2
  X <- matrix(x, ncol=m)
  evals <- array(0,c(n2,m2))
  evals[1:n,1:m] <- X
  y <- stats::fft(evals)
  y <- eigvals*y
  y <- stats::fft(y, inverse = T)/m2/n2
  y <- y[1:n,1:m]
  return(c(Re(y)))
}

#' Matrix Square-root
#'
#' For a symmetric matrix \code{X} calculates the matrix \code{A} such that
#' \code{X} = \code{A}^2.
#'
#' @param X a symmetric matrix
#'
#' @return a symmetric matrix \code{A}
#' @export
#'
#' @examples
#' X <- matrix(runif(9),3)
#' X <- X %*% t(X)
#' mat_root(X)
mat_root <- function(X){
  EIG <- eigen(X)
  lambda <- EIG$values
  V <- EIG$vectors
  lambda <- (lambda + abs(lambda))/2
  return(V %*% diag(sqrt(lambda)) %*% t(V))
}
