#' Preconditioned Conjugate Gradient for R-separable Covariance
#'
#' PCG method for solving \code{A1 X B1 + ... + AR X BR + eps X = Y}
#' using the preconditioner \code{(A1 x B1 + eps I)^(-1)}
#'
#' @param Est list of 3 giving the left-hand side, produced by
#'            \code{\link{scd_est}}
#' @param Y the right-hand side matrix
#' @param eps ridge regularization constant
#' @param Ig the initial guess, random by default
#' @param tol the tolerance used for stopping
#'
#' @return list of 2: \code{U} - the solution X as a matrix;
#'                    \code{iter} - the number of iterations
#' @export
#'
#' @examples
#' X <- array(runif(20*3*4),c(20,3,4))
#' X_new <- array(runif(3*4),c(3,4))
#' Est <- scd_est(X,R=2)
#' Y <- apply_lhs(Est, X_new)
#' X_hat <- pcg(Est, Y, 1e-5)
#' norm(X_new - X_hat$U, type="F")/norm(X_new, type="F")
pcg <- function(Est,Y,eps,Ig=NULL,tol=1e-10){
  K1 <- dim(Y)[1]
  K2 <- dim(Y)[2]
  R <- length(Est$sigma)
  Est_new <- list(sigma=rep(0,R+1),A=array(0,c(R+1,K1,K1)),
                  B=array(0,c(R+1,K2,K2)))
  Est_new$sigma[1:R] <- Est$sigma
  Est_new$A[1:R,,] <- Est$A
  Est_new$B[1:R,,] <- Est$B
  Est_new$sigma[R+1] <- eps
  Est_new$A[R+1,,] <- diag(K1)
  Est_new$B[R+1,,] <- diag(K2)
  mmax = 2000; # the maximum number of iterations
  if(length(Ig)==0) Ig <- matrix(stats::rnorm(K1*K2),ncol=K2)
  U <- Ig
  R <- Y-apply_lhs(Est_new,U)
  E <- rep(0,mmax)
  E[1] <- sqrt(sum(R^2))
  if(!is.finite(E[1])){
    E[1] <- 1
    iter <- mmax
    U <- array(0,c(K1,K2))
  }
  iter = 1
  t1 = 1.0
  D = array(0,c(K1,K2))
  EIG <- eigen(Est_new$B[1,,])
  UU <- EIG$vectors
  alpha1 <- EIG$values
  EIG <- eigen(Est_new$A[1,,])
  VV <- EIG$vectors
  alpha2 <- EIG$values
  HH <- ( (alpha2 %*% t(alpha1))*Est_new$sigma[1] + eps )^(-1)
  while(iter < mmax && E[iter]/E[1] > tol){
    # analytic solution to the Stein's equation
    Z <- VV %*% ( HH * (t(VV) %*% R %*% UU) ) %*% t(UU)
    t1old = t1
    t1 = sum(Z*R)
    beta = t1/t1old
    D = Z+beta*D
    S = apply_lhs(Est_new,D)
    suma = sum(D*S)
    tau = t1/suma
    U = U+tau*D
    R = R-tau*S
    iter = iter+1
    E[iter] = sqrt(sum(R^2))
    if(!is.finite(E[iter])){
      E[iter] <- 0
      U <- array(0,c(K1,K2))
    }
  }
  return(list(U=U,iter=iter))
}

#' Fast Application of R-separable Covariance
#'
#' Calculates 'Est(Y)' = \code{sigma1 A1 Y B1 + ... + sigmaR AR Y BR}.
#'
#' @param Est list of 3, as produced by \code{\link{scd_est}}:
#'           \code{A} - array of size \code{R} x \code{K1} x \code{K1},
#'                      the first \code{R} temporal kernels;
#'           \code{B} - array of size \code{R} x \code{K2} x \code{K2},
#'                      the first \code{R} spatial kernels;
#'           \code{sigma} - vector of length \code{R},
#'                          the separable component scores
#' @param Y matrix argument
#'
#' @return `Est(Y)'
#' @export
#'
#' @examples
#' X <- array(runif(20*3*4),c(20,3,4))
#' Y <- array(runif(3*4),c(3,4))
#' Est <- scd_est(X,R=2)
#' apply_lhs(Est, Y)
apply_lhs <- function(Est,Y){
  R <- length(Est$sigma)
  K1 <- dim(Est$A)[2]
  K2 <- dim(Est$B)[2]
  Z <- array(0,c(K1,K2))
  for(r in 1:R){
    Z <- Z + Est$sigma[r] * ( Est$A[r,,] %*% Y %*% Est$B[r,,] )
  }
  return(Z)
}
