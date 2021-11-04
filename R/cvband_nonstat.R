#' Cross-validation for Bandwidth Selection in Separable-plus-banded Model
#' without Stationarity
#'
#' Similar to \code{cvband}, but the banded part of the model is not assumed to
#' be stationary here. Because of that this function is computationally more
#' demanding compared to \code{cvband}. The purpose of this function
#' is to check whether there might be heteroscedastic noise in the data,
#' corresponding to a separable-plus-diagonal model, i.e. to \code{d=1}. For
#' that, one needs to look past the separable-plus-diagonal, i.e. search for
#' \code{d=0,1,..,maxd}. However, relatively small \code{maxd} should be enough.
#'
#' @param X data set, array of size \code{N} x \code{K1} x \code{K2}
#' @param Folds number of folds
#' @param maxd maximum bandwidth value to check
#' @param mind minimum bandwidth value to check
#'
#' @return matrix with two rows: first row gives the fit, the second gives the
#'         norm, columns correspond to different bandwidths
#' @export
#'
#' @examples
#' N <- 20
#' K1 <- 5
#' K2 <- 7
#' set.seed(517)
#' X <- array(0,c(N,K1,K2))
#' A <- matrix(rnorm(K1^2),K1)
#' B <- matrix(rnorm(K2^2),K2)
#' A <- A %*% t(A)
#' B <- B %*% t(B)
#' A <- mat_root(A)
#' B <- mat_root(B)
#' for(n in 1:N){
#'   X[n,,] <- A %*% X[n,,] %*% B + 2*matrix(rnorm(K1*K2),K1)
#' }
#' cvscores <- cvband_nonstat(X,10,4,0) # two rows: fit and norm
#' cvscores[1,]/cvscores[2,] # cross-validation objective: fit/norm
#' min(localMaxima(cvscores[1,]/cvscores[2,])) # choose bandwidth as the
#'                                             # smallest local maximum
cvband_nonstat <- function(X,Folds=10,maxd=2,mind=0){
  N <- dim(X)[1]
  K1 <- dim(X)[2]
  K2 <- dim(X)[3]
  Nnew <- (N %/% Folds)*Folds
  Ind <- matrix(sample(1:Nnew), nrow=Folds)
  Fit <- array(0,c(Folds,maxd-mind+1))
  Norms <- rep(0,maxd-mind+1)
  for(delta in mind:maxd){
    Est <- spt(X,delta)
    if(delta > 1){
      B <- banded_nonstat_est(X,delta,Est$A1,Est$A2)
      Norms[delta-mind+1] <- frobenius_spb_nonstat(Est$A1,Est$A2,B)
    }else if(delta==1){
      B <- diagonal_band(X,Est$A1,Est$A2)
      pom <- frobenius(Est$A1)^2*frobenius(Est$A2)^2 -
             frobenius(diag(Est$A1))^2*frobenius(diag(Est$A2))^2
      Norms[delta-mind+1] <- sqrt(pom +
                             frobenius(outer(diag(Est$A1),diag(Est$A2)) + B)^2)
    }else{ # delta==0
      Norms[delta-mind+1] <- frobenius(Est$A1)*frobenius(Est$A2)
    }
  }
  for(fold in 1:Folds){
    cat(paste0("Fold: ",fold,"/",Folds),"\r")
    Xtrain <- X[-Ind[fold,],,]
    Xtest <- X[Ind[fold,],,]
    Ntest <- dim(Xtest)[1]
    for(delta in mind:maxd){
      akt_err <- rep(0,Ntest)
      Est <- spt(Xtrain,delta)
      if(delta > 1){
        B <- banded_nonstat_est(Xtrain,delta,Est$A1,Est$A2)
        for(n in 1:Ntest){
          pom <- Est$A1 %*% Xtest[n,,] %*% Est$A2
          akt_err[n] <- sum(pom*Xtest[n,,]) + banded_nonstat_ip(Xtest[n,,],B)
        }
      }else if(delta==1){
        for(n in 1:Ntest){
          B <- diagonal_band(Xtrain,Est$A1,Est$A2)
          pom <- Est$A1 %*% Xtest[n,,] %*% Est$A2
          akt_err[n] <- sum(pom*Xtest[n,,]) + sum(Xtest[n,,]*B*Xtest[n,,])
        }
      }else{ # delta==0
        for(n in 1:Ntest){
          pom <- Est$A1 %*% Xtest[n,,] %*% Est$A2
          akt_err[n] <- sum(pom*Xtest[n,,])
        }
      }
      Fit[fold,delta-mind+1] <- mean(akt_err)
    }
  }
  Fit <- colMeans(Fit)

  RES <- array(0,c(2,maxd-mind+1))
  RES[1,] <- Fit
  RES[2,] <- Norms
  return(RES)
}

#' Inner Product of an Observation with itself through a Non-stationary Banded
#' Covariace
#'
#' @param X observation, a \code{K1} x \code{K2} matrix
#' @param B banded non-stationary covariance in the form of output of
#'          \code{banded_nonstat_est()}
#'
#' @return the inner product, i.e. a number
#' @export
#'
#' @examples
#' X <- array(runif(20*3*4), c(20,3,4))
#' B <- banded_nonstat_est(X,2)
#' banded_nonstat_ip(X[1,,],B)
banded_nonstat_ip <- function(X,B){
  d <- dim(B)[1]
  Tot <- 0
  for(i in 0:(d-1)){
    for(j in 0:(d-1)){
      Xshift <- shift_obs(X,i,j)
      Tot <- Tot + sum(X*Xshift*B[1+i,1+j,])
    }
    for(j in (-d+1):(-1)){
      Xshift <- shift_obs(X,i,j)
      Tot <- Tot + sum(X*Xshift*B[i+1,2*d+j,])
    }
  }
  2*Tot-sum(X^2*B[1,1,])
}

#' Estimate a Non-stationary Banded Covariance
#'
#' Calculates the banded version of the empirical covariance for the purpose
#' of cross-validating bandwidth without assuming stationarity. This can be done
#' simply by 1) calculating the empirical covariance and 2) setting entries
#' outside of a given band to zero, which is however computationally wasteful,
#' because many calculations from step 1 are thrown away in step 2.
#' Instead, here we calculate directly only those entries of the empirical
#' covariance, which will be kept non-zero.
#'
#' @param X data set, array of size \code{N} x \code{K1} x \code{K2}
#' @param d integer >1, the discrete bandwidth
#' @param A1 matrix of size \code{K1} x \code{K1}, the temporal kernel, defaults
#'           to zero
#' @param A2 matrix of size \code{K2} x \code{K2}, the spatial kernel, defaults
#'           to zero
#'
#' @return array B of size \code{d} x \code{2*d-1} x \code{K1*K2}
#'         specifying the banded part of the covariance in the following format:
#'         \code{B[1+i,1+j,]} is a vector giving the off-diagonal of C
#'         corresponding
#'         to shift i in time and shift j in space, while \code{B[1+i,2*d-j,]} is
#'         a vector giving the off-diagonal of C corresponding to shift i in
#'         time and negative shift j in space, suitably adjunct by zeros.
#'         Examples:
#'         \code{B[1,1,]} gives the diagonal of C;
#'         \code{B[3,2,]} gives covariances between locations on the
#'         domain separated by 2 in time and 1 in space in the same direction
#'         (including e.g. \code{cov(X[1,1], X[3,2])} or
#'         \code{cov(X[2,4], X[4,5])});
#'         \code{B[3,2*d-1,]} gives covariances between locations on the domain
#'         separated by 2 in time and 1 in space but in the opposite directions
#'         (including e.g. \code{cov(X[3,1], X[1,2])}).
#' @export
#'
#' @examples
#' X <- array(runif(20*3*4), c(20,3,4))
#' B <- banded_nonstat_est(X,2)
#' B[1,1,] # diagonal of the empirical covariance
banded_nonstat_est <- function(X,d,A1=NULL,A2=NULL){
  N <- dim(X)[1]
  K1 <- dim(X)[2]
  K2 <- dim(X)[3]
  B <- array(0,c(d, 2*d-1, K1*K2))
  if(length(A1)==0) A1 <- array(0,c(K1,K1))
  if(length(A2)==0) A2 <- array(0,c(K2,K2))
  for(i in 0:(d-1)){
    for(j in 0:(d-1)){
      for(n in 1:N){
        Xshift <- shift_obs(X[n,,],i,j)
        B[i+1,j+1,] <- B[i+1,j+1,] + c(Xshift*X[n,,])
      }
      B[i+1,j+1,] <- B[i+1,j+1,]/N
      a1 <- rep(0,K1)
      a2 <- rep(0,K2)
      ind <- 1:(K1-i)
      a1[ind] <- A1[cbind(ind,ind+i)]
      ind <- 1:(K2-j)
      a2[ind] <- A2[cbind(ind,ind+j)]
      B[i+1,j+1,] <- B[i+1,j+1,] - outer(a1,a2)
    }
  }
  for(i in 1:(d-1)){
    for(j in (-d+1):(-1)){
      for(n in 1:N){
        Xshift <- shift_obs(X[n,,],i,j)
        B[i+1,2*d+j,] <- B[i+1,2*d+j,] + c(Xshift*X[n,,])
      }
      B[i+1,2*d+j,] <- B[i+1,2*d+j,]/N
      a1 <- rep(0,K1)
      a2 <- rep(0,K2)
      ind <- 1:(K1-i)
      a1[ind] <- A1[cbind(ind,ind+i)]
      ind <- 1:(K2+j)
      a2[ind] <- A2[cbind(ind,ind-j)]
      B[i+1,2*d+j,] <- B[i+1,2*d+j,] - outer(a1,a2)
    }
  }
  return(B)
}

#' Shift an Observation to Calculate Off-diagonals of a Banded Covariance
#'
#' @param X observation, a matrix
#' @param i row shift, has to be positive
#' @param j column shift, has to be negative (due to symmetry of the covariance,
#'          it is enough to allow for negative row shift, forcing the column
#'          shift to be positive, or vice versa)
#'
#' @return shifted observation adjunct by zeros to the size of \code{X}
#' @export
#'
#' @examples
#' ( X <- matrix(9:1,3) )
#' shift_obs(X,1,2)
#' shift_obs(X,1,-2)
shift_obs <- function(X,i,j){
  K1 <- dim(X)[1]
  K2 <- dim(X)[2]
  Xshift <- array(0,c(K1,K2))
  if(j < 0){
    Xshift[1:(K1-i),(1-j):K2] <- X[(i+1):K1,1:(K2+j)]
  }else{
    Xshift[1:(K1-i),1:(K2-j)] <- X[(i+1):K1,(j+1):K2]
  }
  return(Xshift)
}
