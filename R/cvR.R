#' Cross-validation for Degree-of-Separability Selection
#'
#' Calculates the objective value of the cross-validation criterion to select
#' the degree of separability R=1,...,maxR best fitting the data.
#'
#' @param X data set, array of size \code{N} x \code{K1} x \code{K2}
#' @param Folds number of folds
#' @param maxR maximum degree-of-separability considered
#'             (defaults to log-sample size)
#' @param maxiter maximum number of iterations for \code{scd_est()}
#'
#' @return CV objective
#' @export
#'
#' @examples
#' X <- array(runif(20*3*4),c(20,3,4))
#' cvscores <- cvR(X,10,5)
#' min(localMaxima(-cvscores)) # best R is chosen as the smallest local minimum
#'                             # of the CV objective
cvR <- function(X,Folds=10,maxR=7,maxiter=10){
  N <- dim(X)[1]
  K1 <- dim(X)[2]
  K2 <- dim(X)[3]
  Nnew <- (N %/% Folds)*Folds
  Ind <- matrix(sample(1:Nnew), nrow=Folds)
  ERR <- array(0,c(Folds,maxR))
  Norms <- rep(0,maxR)
  for(fold in 1:Folds){
    cat(paste0("Fold: ",fold,"/",Folds),"\r")
    Xtrain <- X[-Ind[fold,],,]
    Xtest <- X[Ind[fold,],,]
    Ntest <- dim(Xtest)[1]

    Res <- scd_est(Xtrain,maxR,maxiter)
    for(r in 1:maxR){
      Norms[r] <- sum(Res$sigma[1:r]^2)
      akt_err <- rep(0,Ntest)
      akt_Res <- Res
      akt_Res$A <- Res$A[1:r,,,drop=F]
      akt_Res$B <- Res$B[1:r,,,drop=F]
      akt_Res$sigma <- Res$sigma[1:r]
      for(n in 1:Ntest){
        akt_err[n] <- sum(apply_lhs(akt_Res,Xtest[n,,])*Xtest[n,,])
      }
      ERR[fold,r] <- mean(akt_err)
    }
  }
  ERR <- colMeans(ERR)
  # CVscores <- Norms - 2*ERR
  CVscores <- matrix(c(ERR,Norms),nrow=2,byrow=T)
  return(CVscores)
}

#' Cross-validation for Degree-of-separability Selection via Prediction
#'
#' Examines the prediction performance of R-separable model for
#' the degree of separability R=1,...,maxR via a cross-validated scheme:
#' for every surface in the current hold-out sample, part of the surface
#' (a fixed number of rows and columns given by \code{perc}) is predicted based
#' on the remainder, and relative prediction error is calculated. Two random
#' splits of every observation are predicted.
#'
#' @param X data set, array of size \code{N} x \code{K1} x \code{K2}
#' @param Folds number of folds
#' @param maxR maximum degree-of-separability considered
#' @param maxiter maximum number of iterations for scd_est()
#' @param perc number in (0,1) giving the portion of rows and columns held out
#'             for prediction, defaults to 1/2, in that case 3/4 of every
#'             held-out surface is predicted based on the remaining 1/4
#'
#' @return vector of prediction performances for different values of R
#' @export
#'
#' @examples
#' X <- array(runif(20*3*4),c(20,3,4))
#' cvscores <- cvR_pred(X,10,5)
#' min(localMaxima(-cvscores)) # best R is chosen as the smallest local minimum
#'                             # of the CV objective
cvR_pred <- function(X,Folds=10,maxR=7,maxiter=10,perc=NULL){
  N <- dim(X)[1]
  K1 <- dim(X)[2]
  K2 <- dim(X)[3]
  if(is.null(perc)) perc <- c(1/2,1/2)
  Nnew <- (N %/% Folds)*Folds
  Ind <- matrix(sample(1:Nnew), nrow=Folds)
  ERR <- array(0,c(Folds,maxR))
  for(fold in 1:Folds){
    cat(paste0("Fold: ",fold,"/",Folds),"\r")
    Xtrain <- X[-Ind[fold,],,]
    Xtest <- X[Ind[fold,],,]
    Ntest <- dim(Xtest)[1]

    Res <- scd_est(Xtrain,maxR,maxiter)
    EIG <- eigvals(Res)
    if(EIG$min < 0) eps <- -EIG$min + 1e-5 else eps <- 1e-5
    for(r in 1:maxR){
      akt_err <- array(0,c(2,Ntest))
      akt_Res <- Res
      akt_Res$A <- Res$A[1:r,,,drop=F]
      akt_Res$B <- Res$B[1:r,,,drop=F]
      akt_Res$sigma <- Res$sigma[1:r]
      for(n in 1:Ntest){
        I <- sample(1:K1,max(floor(K1*perc),1))
        J <- sample(1:K2,max(floor(K2*perc),1))
        Xhat <- blup(Xtest[n,,],akt_Res,I,J,eps)
        akt_err[1,n] <- frobenius(Xtest[n,,] - Xhat)/frobenius(Xtest[n,,])
        # second split of the same observation
        I <- sample(1:K1,max(floor(K1*perc),1))
        J <- sample(1:K2,max(floor(K2*perc),1))
        Xhat <- blup(Xtest[n,,],akt_Res,I,J,eps)
        akt_err[2,n] <- frobenius(Xtest[n,,] - Xhat)/frobenius(Xtest[n,,])
      }
      ERR[fold,r] <- mean(akt_err)
    }
  }
  ERR <- colMeans(ERR)
  return(ERR)
}

#' Best Linear Unbiased Prediction in R-separable Model
#'
#' @param Xnew new observation to be predicted, given as a matrix
#' @param Res list of 3, R-separable covariance estimated by \code{scd_est()}
#' @param I indices of unobserved rows
#' @param J indices of unobserved columns
#' @param eps ridge regularization constant
#'
#' @return prediction o \code{Xnew} given as a matrix
#' @export
#'
#' @examples
#' X <- array(runif(20*3*4),c(20,3,4))
#' Res <- scd_est(X,2)
#' Xnew <- array(runif(3*4), c(3,4))
#' blup(Xnew, Res, 3, 4, 1e-5)
blup <- function(Xnew,Res,I,J,eps){
  K1 <- dim(Xnew)[1]
  K2 <- dim(Xnew)[2]
  C22 <- Res
  C22$A <- Res$A[,-I,-I,drop=F]
  C22$B <- Res$B[,-J,-J,drop=F]
  Xobs <- Xnew[-I,-J,drop=F]
  Z <- pcg(C22,Xobs,eps)
  Zfull <- array(0,c(K1,K2))
  Zfull[-I,-J] <- Z$U
  Xhat <- apply_lhs(Res,Zfull)
  Xhat[-I,-J] <- Xobs
  return(Xhat)
}

#' Maximum and Minimum Eigenvalue of an R-separable Covariance
#'
#' @param Est R-separable covariance given as a list of 3:
#'           \code{A} - array of size R x K1 x K1, the first R temporal kernels;
#'            \code{B} - array of size R x K2 x K2, the first R spatial kernels;
#'            \code{sigma} - vector of length R, the separable component scores
#' @param tol stopping criterion
#' @param maxiter maximum number of iterations
#'
#' @return list of 2: \code{min} - the minimum eigenvalue;
#'                    \code{max} - the maximum eigenvalue
#' @export
#'
#' @examples
#' X <- array(runif(20*3*4),c(20,3,4))
#' Res <- scd_est(X,2)
#' eigvals(Res)
eigvals <- function(Est,tol=1e-10,maxiter=10000){
  # calculates the largest and the smallest eigenvalue of the estimated covariance
  K1 <- dim(Est$A)[2]
  K2 <- dim(Est$B)[2]
  V <- matrix(stats::runif(K1*K2),ncol=K2) # initial guess
  V <- V/frobenius(V)
  eps <- 1
  iter <- 0
  # TODO can V explode and the while condition fail?
  while(eps > tol && iter < maxiter){
    iter <- iter+1
    V_old <- V
    V <- apply_lhs(Est,V)
    V <- V/frobenius(V)
    eps <- frobenius(V - V_old)
  }
  lambda_max <- sum(V*apply_lhs(Est,V))

  V <- matrix(stats::runif(K1*K2),ncol=K2) # initial guess
  V <- V/frobenius(V)
  eps <- 1
  # TODO can V explode and the while condition fail?
  while(eps > tol && iter < maxiter){
    iter <- iter+1
    V_old <- V
    V <- apply_lhs(Est,V)
    V <- lambda_max*V_old - V
    V <- V/frobenius(V)
    eps <- frobenius(V - V_old)
  }
  lambda_min <- sum(V*apply_lhs(Est,V))
  return(list(min=lambda_min, max=lambda_max))
}
