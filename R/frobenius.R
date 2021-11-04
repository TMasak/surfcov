#' Frobenius Norm of an Array
#'
#' @param X array of arbitrary dimension and size
#'
#' @return a number, the Frobenius norm of \code{X}
frobenius <- function(X){
  return(sqrt(sum(X^2)))
}

#' Frobenius Norm of a Difference of Two Separable-plus-banded Covariances
#'
#' @param Chat separable-plus-banded model \code{A1} x \code{A2} + \code{B},
#'             i.e. a list of 3 elements:
#'             matrix \code{A1}, matrix \code{A2}, the symbol of \code{B}
#' @param C as \code{Chat} above
#' @param relative flag whether the error relative to \code{C} should be computed,
#'                 defaults to \code{TRUE}
#' @return a number, the Frobenius norm
#' @export
#'
#' @examples
#' X <- array(runif(20*4*5), c(20,4,5))
#' Y <- array(runif(20*4*5), c(20,4,5))
#' Chat <- spb(X,2)
#' C <- spb(Y,2)
#' frobenius_spb(Chat,C)
frobenius_spb <- function(Chat, C=NULL, relative=T){
  K1 <- dim(Chat$A1)[1]
  K2 <- dim(Chat$A2)[1]
  FROB <- 0
  FROB_down <- 0
  delta1 <- dim(Chat$B)[1]
  if(is.null(C)){
    A1 <- matrix(0, K1, K1)
    A2 <- matrix(0, K2, K2)
    B <- matrix(0, delta1, delta1)
    C <- list(A1=A1, A2=A2, B=B)
  }
  delta2 <- dim(C$B)[1]
  Band1 <- matrix(0, K1, K2)
  Band2 <- matrix(0, K1, K2)
  Band1[1:delta1,1:delta1] <- Chat$B
  Band2[1:delta2,1:delta2] <- C$B
  for(i in 0:(K1-1)){
    ind <- 1:(K1-i)
    a1 <- Chat$A1[cbind(ind+i,ind)]
    a2 <- C$A1[cbind(ind+i,ind)]
    for(j in 0:(K2-1)){
      ind <- 1:(K2-j)
      b1 <- Chat$A2[cbind(ind+j,ind)]
      b2 <- C$A2[cbind(ind+j,ind)]
      beta <- Band1[i+1,j+1] - Band2[i+1,j+1]
      times <- 2
      if(i==0 && j==0) times <- 1
      if(i>0 && j>0) times <- 4
      # we can only use fast calculation if we have enough elements for
      # Gram-Schmidt orthogonalization
      if(i <= K1-3 && j <= K2-3){
        FROB <- FROB + times*frobenius_piece(a1,b1,a2,b2,beta)
        FROB_down <- FROB_down +
                     times*frobenius_piece(a2,b2,NULL,NULL,Band2[i+1,j+1])
      } else {
        FROB <- FROB + times*frobenius_piece_slow(a1,b1,a2,b2,beta)
        FROB_down <- FROB_down +
                     times*frobenius_piece_slow(a2,b2,NULL,NULL,Band2[i+1,j+1])
      }
    }
  }
  if(relative) FROB <- sqrt(FROB/FROB_down)
  else FROB <- sqrt(FROB)
  return(FROB)
}

#' Frobenius Norm of a Separable-plus-banded Covariance with Non-stationary
#' Banded Part
#'
#' @param A1 matrix, the temporal kernel
#' @param A2 matrix, the spatial kernel
#' @param B 3-dimensional array specifying the banded part in the form of an
#'          output of \code{\link{banded_nonstat_est}}
#'
#' @return the value of the Frobenius norm
#' @export
#'
#' @examples
#' N <- 20
#' K1 <- 10
#' K2 <- 10
#' delta <- 2
#' X <- array(runif(N*K1*K2),c(N,K1,K2))
#' Est <- spt(X,delta)
#' B <- banded_nonstat_est(X,delta,Est$A1,Est$A2)
#' frobenius_spb_nonstat(Est$A1,Est$A2,B)

frobenius_spb_nonstat <- function(A1, A2, B){
  K1 <- dim(A1)[1]
  K2 <- dim(A2)[1]
  FROB <- 0
  d <- dim(B)[1]
  for(i in 0:(K1-1)){
    ind <- 1:(K1-i)
    a1 <- rep(0,K1)
    a1[ind] <- A1[cbind(ind+i,ind)]
    for(j in 0:(K2-1)){
      ind <- 1:(K2-j)
      a2 <- rep(0,K2)
      a2[ind] <- A2[cbind(ind+j,ind)]
      if(i < d && j < d){
        bb <- B[1+i,1+j,]
        FROB <- FROB + frobenius(outer(a1,a2) + bb)^2
        if(i==0 && j==0) DIAG <- frobenius(outer(a1,a2) + bb)^2
      }else{
        FROB <- FROB + frobenius(a1)^2 * frobenius(a2)^2
      }
    }
  }
  for(i in 1:(K1-1)){
    ind <- 1:(K1-i)
    a1 <- rep(0,K1)
    a1[ind] <- A1[cbind(ind+i,ind)]
    for(j in (-K2+1):(-1)){
      ind <- 1:(K2+j)
      a2 <- rep(0,K2)
      a2[ind] <- A2[cbind(ind-j,ind)]
      a2 <- rev(a2)
      if(i < d && -j < d){
        bb <- B[1+i,2*d+j,]
        FROB <- FROB + frobenius(outer(a1,a2) + bb)^2
      }else{
        FROB <- FROB + frobenius(a1)^2 * frobenius(a2)^2
      }
    }
  }
  # everything but the diagonal is there twice
  FROB <- 2*FROB - DIAG
  return(sqrt(FROB))
}

#' Frobenius Norm of a rank-3 Matrix (Vanilla)
#'
#' Calculates the Frobenius norm of a single sub-matrix of the difference of the
#' Van Loan's permutations of a separable-plus-banded covariance, which must be
#' at most rank 3, with left factors \code{a1}, \code{b1} and \code{1}, and
#' right factors \code{a2}, \code{b2}, 1. This is a vanilla version used for
#' matrices with \code{ncol <= 3} or \code{nrow <= 3}, because the Gram-Schmidt
#' orthogonalization cannot be used for such a small matrix.
#'
#' @param a1 left singular vector of the minuend
#' @param b1 right singular vector of the minuend
#' @param a2 left singular vector of the subtrahend
#' @param b2 left singular vector of the minuend
#' @param beta corresponding to the constant matrix \code{1} otimes \code{1}
#'
#' @return a number, the Frobenius norm
#'
#' @examples
#' a1 <- runif(5)
#' b1 <- runif(10)
#' a2 <- runif(5)
#' b2 <- runif(10)
#' beta <- runif(1)
#' frobenius_piece(a1,b1,a2,b2,beta)
frobenius_piece_slow <- function(a1,b1,a2,b2,beta){
  if(length(a2)==0) return(frobenius(a1 %*% t(b1) + beta)^2)
  else return(frobenius(a1 %*% t(b1) - a2 %*% t(b2) + beta)^2)
}

#' Frobenius Norm of a rank-3 Matrix
#'
#' Calculates the Frobenius norm of a single sub-matrix of the difference of the
#' Van Loan's permutations of a separable-plus-banded covariance, which must be
#' at most rank 3, with left factors \code{a1}, \code{b1} and \code{1}, and
#' right factors \code{a2}, \code{b2}, \code{1}. Gram-Schmidt orthogonalization
#' can be used to calculate the Frobenius norm of such a matrix of size \code{K}
#' x \code{K} in only 3*\code{K} operations, where \code{K} is the dimension of
#' the matrix.
#'
#' @param a1 left singular vector of the minuend
#' @param b1 right singular vector of the minuend
#' @param a2 left singular vector of the subtrahend
#' @param b2 left singular vector of the minuend
#' @param beta corresponding to the constant matrix 1 otimes 1
#'
#' @return a number, the Frobenius norm
#' @export
#'
#' @examples
#' a1 <- runif(5)
#' b1 <- runif(10)
#' a2 <- runif(5)
#' b2 <- runif(10)
#' beta <- runif(1)
#' frobenius_piece(a1,b1,a2,b2,beta)
frobenius_piece <- function(a1,b1,a2,b2,beta){
  K1 <- length(a1)
  K2 <- length(b1)
  if(is.null(a2)){
    a2 <- rep(0,K1)
    b2 <- rep(0,K2)
  }
  M <- array(0,c(K1,3))
  M[,1] <- a1
  M[,2] <- a2
  M[,3] <- sqrt(abs(beta))
  pom <- qr(M)
  S1 <- pom$qr[1:3,1:3]
  Ind <- lower.tri(S1)
  S1[Ind] <- 0

  M <- array(0,c(K2,3))
  M[,1] <- b1
  M[,2] <- -b2
  M[,3] <- sign(beta)*sqrt(abs(beta))
  pom <- qr(M)
  S2 <- pom$qr[1:3,1:3]
  Ind <- lower.tri(S2)
  S2[Ind] <- 0
  return(frobenius( S1 %*% t(S2) )^2)
}
