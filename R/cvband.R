#' Cross-validation for Bandwidth Selection in Separable-plus-banded Model
#'
#' Calculates the objective value of the cross-validation criterion for
#' bandwidth d between mind and maxd.
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
#' N <- 500
#' K1 <- 5
#' K2 <- 7
#' set.seed(517)
#' X <- array(rnorm(N*K1*K2),c(N,K1,K2))
#' A <- matrix(rnorm(K1^2),K1)
#' B <- matrix(rnorm(K2^2),K2)
#' A <- A %*% t(A)
#' B <- B %*% t(B)
#' A <- mat_root(A)
#' B <- mat_root(B)
#' for(n in 1:N){
#'   X[n,,] <- A %*% X[n,,] %*% B + matrix(rnorm(K1*K2),K1)
#' }
#' maxd <- 3 # max range tobe
#' mind <- 1
#' cvscores <- cvband(X,10,maxd,mind) # two rows: fit and norm
#' cvscores[1,]/cvscores[2,] # cross-validation objective: fit/norm
#' cvscores[2,]^2 - 2*cvscores[1,] # alternative objective: norm^2 - 2*fit
#' # choose bandwidth as the smallest local maximum:
#' min(localMaxima(cvscores[1,]/cvscores[2,])) + mind - 1
#' min(which.min(cvscores[2,]^2 - 2*cvscores[1,])) + mind - 1
#'   # "+mind-1" to take care if mind in not equal to one
cvband <- function(X,Folds=10,maxd=20,mind=1){
  N <- dim(X)[1]
  K1 <- dim(X)[2]
  K2 <- dim(X)[3]
  Nnew <- (N %/% Folds)*Folds
  Ind <- matrix(sample(1:Nnew), nrow=Folds)
  Fit <- array(0,c(Folds,maxd-mind+1))
  Norms <- rep(0,maxd-mind+1)
  for(delta in mind:maxd){
    Chat <- spb(X,delta)
    Norms[delta-mind+1] <- frobenius_spb(Chat,relative=F)
  }
  for(fold in 1:Folds){
    cat(paste0("Fold: ",fold,"/",Folds),"\r")
    Xtrain <- X[-Ind[fold,],,]
    Xtest <- X[Ind[fold,],,]
    Ntest <- dim(Xtest)[1]
    for(delta in mind:maxd){
      Chat <- spb(Xtrain,delta)
      akt_err <- rep(0,Ntest)
      for(n in 1:Ntest){
        eigvals <- Re(stats::fft(to_book_format(Chat$B,K1,K2)))
        pom <- Chat$A1 %*% Xtest[n,,] %*% Chat$A2
        pom <- pom + BXfast(eigvals,c(Xtest[n,,]))
        akt_err[n] <- sum(pom*Xtest[n,,])
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

#' Cross-validation for Bandwidth Selection in Separable-plus-banded Model via
#' Prediction
#'
#' Examines the prediction performance of the separable-plus-banded model for
#' different (discrete) bandwidths d=\code{mind},...,\code{maxd} via a
#' cross-validated scheme: for every surface in the current hold-out sample,
#' a part of the surface (a fixed number of last few rows and columns
#' give by \code{perc}) is predicted based on the remainder of the given surface,
#' and relative prediction error is calculated.
#'
#' @param X data, array of size N x K1 x K2
#' @param Folds number of folds
#' @param maxd maximum bandwidth value to check
#' @param mind minimum bandwidth value to check
#' @param perc vector of 2 values in (0,1) giving the ratio of rows and columns
#'             to be predicted
#'
#' @return vector of cross-validated prediction errors for
#'         d=\code{mind},...,\code{maxd}
#' @export
#'
#' @examples
#' N <- 4
#' K1 <- 4
#' K2 <- 4
#' set.seed(517)
#' X <- array(0,c(N,K1,K2))
#' A <- matrix(rnorm(K1^2),K1)
#' B <- matrix(rnorm(K2^2),K2)
#' A <- A %*% t(A)
#' B <- B %*% t(B)
#' A <- mat_root(A)
#' B <- mat_root(B)
#' for(n in 1:N){
#'   X[n,,] <- A %*% X[n,,] %*% B + matrix(rnorm(K1*K2),K1)
#' }
#' ( cvscores <- cvband_pred(X,2,2,0) )
#' min(localMaxima(-cvscores))
cvband_pred <- function(X,Folds=10,maxd=20,mind=1, perc=c(1/4,1/4)){
  N <- dim(X)[1]
  K1 <- dim(X)[2]
  K2 <- dim(X)[3]
  Nnew <- (N %/% Folds)*Folds
  Ind <- matrix(sample(1:Nnew), nrow=Folds)
  ERR <- array(0,c(Folds,maxd-mind+1))
  m1 <- K1 - floor(K1*perc[1])
  m2 <- K2 - floor(K2*perc[2])
  for(fold in 1:Folds){
    cat(paste0("Fold: ",fold,"/",Folds),"\r")
    Xtrain <- X[-Ind[fold,],,]
    Xtest <- X[Ind[fold,],,]
    Ntest <- dim(Xtest)[1]
    for(delta in mind:maxd){
      Chat <- spb(Xtrain,delta)
      Chat <- seppband_psd(Chat)
      akt_err <- rep(0,Ntest)
      for(n in 1:Ntest){
        Xhat <- blup_adi(Xtest[n,,],Chat$A1,Chat$A2,Chat$B,m1,m2)
        akt_err[n] <- norm(Xhat-Xtest[n,,], type="F")/norm(Xtest[n,,], type="F")
      }
      ERR[fold,delta-mind+1] <- mean(akt_err)
    }
  }
  ERR <- colMeans(ERR)
  return(ERR)
}

#' Best Linear Unbiased Prediction in Separable-plus-banded Model
#'
#' @param Xnew new observation to be predicted, a K1 x K2 matrix with last
#'             K1-m1+1 rows and K2-m2+1 columns equal to zero (unobserved)
#' @param A1 temporal kernel of the separable-plus-banded model
#' @param A2 spatial kernel of the separable-plus-banded model
#' @param band symbol of the banded part of the separable-plus-banded model
#' @param m1 \code{Xnew} observed until this row index
#' @param m2 Xnew observed until this column index
#' @param theta regularization for adi()
#'
#' @return prediction of Xnew, i.e. matrix of size K1 x K2 where the last
#'         columns and rows are no longer zero
#' @export
#'
#' @examples
#' N <- 100
#' K1 <- 5
#' K2 <- 7
#' A1 <- matrix(rnorm(K1^2),K1)
#' A2 <- matrix(rnorm(K2^2),K2)
#' A1 <- A1 %*% t(A1)
#' A2 <- A2 %*% t(A2)
#' A1 <- mat_root(A1)
#' A2 <- mat_root(A2)
#' X <- A1 %*% array(rnorm(K1*K2),c(K1,K2)) %*% A2 + matrix(rnorm(K1*K2),K1)
#' Xhat <- blup_adi(X,A1,A2,as.matrix(1),4,5)
#' norm(X-Xhat, type="F")/norm(X, type="F")
blup_adi <- function(Xnew,A1,A2,band,m1,m2,theta=1e-5){
  K1 <- dim(A1)[1]
  K2 <- dim(A2)[2]
  d <- dim(band)[1]
  I <- (m1+1):K1
  J <- (m2+1):K2
  A1obs <- A1[-I,-I,drop=F]
  A2obs <- A2[-J,-J,drop=F]
  if(m1 < d) band <- band[1:m1,,drop=F]
  if(m2 < d) band <- band[,1:m2,drop=F]
  # TODO check whether A1 or A2 are not low-rank, if so, adjust `rho' for adi()
  Xobs <- Xnew[-I,-J]
  z <- adi(A1obs,A2obs,band,c(Xobs),theta)
  z <- z$x
  Z <- array(0,c(K1,K2))
  Z[-I,-J] <- z
  # TODO precalculate `eigvals' once in cvband_pred()
  eigvals <- Re(stats::fft(to_book_format(band,K1,K2)))
  Y <- A1 %*% Z %*% A2 + BXfast(eigvals,c(Z))
  Y[-I,-J] <- Xobs
  return(Y)
}

#' Positivize a Separable-plus-banded Estimator
#'
#' Changes the separable-plus-banded model in list \code{Res} so that all three
#' kernels of the model are positive semi-definite.
#'
#' @param Res list of 3 specifying the separable-plus-banded estimator:
#'            A1 - the temporal kernel;
#'            A2 - the spatial kernel;
#'            B - the symbol of the banded part
#'
#' @return list of 3, the same as \code{Res}
#' @export
#'
#' @examples
#' X <- array(runif(20*3*4),c(20,3,4))
#' Res <- spb(X,1)
#' Res_psd <- seppband_psd(Res)
seppband_psd <- function(Res){
  K1 <- dim(Res$A1)[1]
  K2 <- dim(Res$A2)[2]
  # temporal kernel
  EIG <- eigen(Res$A1)
  EIG$values <- pmax(EIG$values,1e-6)
  Res$A1 <- EIG$vectors %*% diag(EIG$values) %*% t(EIG$vectors)
  # spatial kernel
  EIG <- eigen(Res$A2)
  EIG$values <- pmax(EIG$values,1e-6)
  Res$A2 <- EIG$vectors %*% diag(EIG$values) %*% t(EIG$vectors)
  # banded part
  eigvals <- Re(stats::fft(to_book_format(Res$B,K1,K2)))
  eigvals <- pmax(eigvals,1e-6)
  Res$B <- Re(stats::fft(eigvals,inverse = T)/(4*K1*K2))
  Res$B <- Res$B[1:K1,1:K2]
  return(Res)
}
