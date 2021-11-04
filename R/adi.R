#' Alternating Direction Implicit Method for Separable-plus-banded Covariance
#'
#' Calculates the solution to the inverse problem involving the
#' separable-plus-banded structure as the left-hand side.
#'
#' @param A1 temporal kernel of the separable-plus-banded model
#' @param A2 spatial kernel of the separable-plus-banded model
#' @param band symbol of the banded part of the separable-plus-banded model
#' @param y right-hand side vector
#' @param theta regularization
#' @param rho initial value of the shift parameter, defaults to the square-root
#'            of the smaller one of the condition numbers of \code{A1} and
#'            \code{A2}. It should should be kept as default unless A1 or A2 is
#'            (exactly) low-rank.
#' @param adapt whether the shift should change in between iterations
#' @param maxiter maximum number of iterations
#' @param tol relative tolerance - stopping criterion
#'
#' @return list of 3 elements: the solution as a vector, the number of ADI
#'         iterations and the vector of PCG iterations (one number per one ADI
#'         iteration)
#' @export
#'
#' @examples
#' K <- 30
#' A1 <- brownian_cov(K)
#' A1 <- A2 <- A1/sum(diag(A1))
#' B <- matrix(c(4,2,2,1)/4,2)/K^2 # to have trace 1
#' theta <- 1e-5
#' eigvals <- Re(stats::fft(to_book_format(B,K,K)))
#' B <- Re(stats::fft(eigvals,inverse = TRUE)/(4*K^2))
#' B <- B[1:K,1:K]
#' x <- runif(K^2)
#' X <- matrix(x,ncol=K)
#' y <- c(A1 %*% X %*% A2) + BXfast(eigvals,x) + theta*x
#' x_adi <- adi(A1, A2, B, y, theta, adapt=TRUE,maxiter=100,tol=10^-7)
#' sum((x-x_adi$x)^2)/sum(x^2)
adi <- function(A1, A2, band, y, theta=1e-5, rho=NULL, adapt=TRUE, maxiter=200,
                tol=10^-9){
  K1 <- dim(A1)[1]
  K2 <- dim(A2)[1]
  iter <- 0
  eps <- 1
  x <- stats::runif(K1*K2)
  PCGiter <- rep(0,maxiter)
  # pre-calculate eigendecompositions
  A <- A1; A1 <- A2; A2 <- A # Kronecker product definition's artifact
  EIG <- eigen(A1)
  U <- EIG$vectors
  alpha1 <- EIG$values
  EIG <- eigen(A2)
  V <- EIG$vectors
  alpha2 <- EIG$values
  if(length(rho)==0){
    rho <- sqrt(max(alpha1[1]*alpha1[K2],alpha2[1]*alpha2[K1])) + 1e-15
  }
  while(iter < maxiter && eps > tol){
    x_old <- x
    iter <- iter + 1
    # calculate x^{(k+1/2)}
    eigvals <- Re(stats::fft(to_book_format(band,K1,K2)))
    y_temp <- y - BXfast(eigvals,x_old) + rho*x_old    # RHS for the first half
    Y <- matrix(y_temp, ncol=K2)                       #        of the iteration
    H <- ( (alpha2 %*% t(alpha1)) + rho + theta )^(-1) #
    X <- V %*% ( H * (t(V) %*% Y %*% U) ) %*% t(U)     # analytic solution to
    x_half <- c(X)                                     #    the Stein's equation
    # calculate x^{(k+1)}
    y_temp <- y - c( A2 %*% X %*% A1 ) + rho*x_half
    band_reg <- band
    band_reg[1,1] <- band_reg[1,1] + rho + theta
    eigvals_prec <- calc_prec(to_book_format(band_reg,K1,K2))
    eigvals <- Re(stats::fft(to_book_format(band_reg,K1,K2)))
    PCG <- pcgtt(eigvals,y_temp,stats::runif(K1*K2),eigvals_prec,tol)
    x <- PCG$u
    PCGiter[iter] <- PCG$iter
    eps <- sum((x-x_old)^2)/sum(x_old^2)
    if(!is.finite(eps)){
      eps <- 0
      x <- rep(0,length(x))
    }
    if(adapt){
      rho <- min(rho, eps)
    }
    # print(c(iter,eps))
  }
  return(list(x=x, iter=iter, PCGiter=PCGiter[1:iter]))
}

#' Alternating Direction Implicit Method for Separable-plus-diagonal Covariance
#'
#' Calculates the solution to the inverse problem involving the
#' separable-plus-diagonal structure as the left-hand side.
#'
#' @param A1 temporal kernel of the separable-plus-diagonal model
#' @param A2 spatial kernel of the separable-plus-diagonal model
#' @param B_diag the diagonal of the model
#' @param y right-hand side vector
#' @param theta regularization
#' @param rho initial value of the shift parameter, defaults to the square-root
#'            of the smaller one of the condition numbers of \code{A1} and
#'            \code{A2}. It should should be kept as default unless A1 or A2 is
#'            (exactly) low-rank.
#' @param adapt whether the shift should change in between iterations
#' @param maxiter maximum number of iterations
#' @param tol relative tolerance - stopping criterion
#'
#' @return list of 3 elements: the solution as a vector, the number of ADI
#'         iterations and the vector of PCG iterations (one number per one ADI
#'         iteration)
#' @export
#'
#' @examples
#' K <- 30
#' A1 <- brownian_cov(K)
#' A1 <- A2 <- A1/sum(diag(A1))
#' B_diag <- array(runif(K^2),c(K,K))
#' B_diag <- B_diag/sum(B_diag)
#' theta <- 1e-5
#' x <- runif(K^2)
#' X <- matrix(x,ncol=K)
#' y <- c(A1 %*% X %*% A2) + c(B_diag*X) + theta*x
#' x_adi <- adi_diag(A1, A2, B_diag, y, theta, adapt=TRUE,maxiter=100,tol=10^-7)
#' sum((x-x_adi$x)^2)/sum(x^2)
adi_diag <- function(A1, A2, B_diag, y, theta=0, rho=NULL, adapt=F, maxiter=200,
                     tol=10^-6){
  K1 <- dim(A1)[1]
  K2 <- dim(A2)[1]
  iter <- 0
  eps <- 1
  x <- stats::runif(K1*K2)
  # pre-calculate eigen-decompositions
  A <- A1; A1 <- A2; A2 <- A # Kronecker product definition's artifact
  EIG <- eigen(A1)
  U <- EIG$vectors
  alpha1 <- EIG$values
  EIG <- eigen(A2)
  V <- EIG$vectors
  alpha2 <- EIG$values
  if(length(rho)==0){
    rho <- sqrt(max(alpha1[1]*alpha1[K2],alpha2[1]*alpha2[K1])) + 1e-15
  }
  while(iter < maxiter && eps > tol){
    x_old <- x
    iter <- iter + 1
    # calculate x^{(k+1/2)}
    y_temp <- y - B_diag*x_old + rho*x_old              # RHS for the first
    Y <- matrix(y_temp, ncol=K2)                        #  half of the iteration
    H <- ( (alpha2 %*% t(alpha1)) + rho + theta )^(-1)  #
    X_half <- V %*% ( H * (t(V) %*% Y %*% U) ) %*% t(U) # analytic solution to
    # calculate x^{(k+1)}                               #   the Stein's equation
    y_temp <- y - c( A2 %*% X_half %*% A1 ) + rho*c(X_half)
    Y <- matrix(y_temp, ncol=K2)
    B_akt <- (B_diag + rho + theta)^(-1)
    X <- B_akt*Y
    x <- c(X)
    eps <- sum((x-x_old)^2)/sum(x_old^2)
    if(!is.finite(eps)){
      eps <- 0
      x <- rep(0,length(x))
    }
    if(adapt){
      rho <- min(rho, eps)
    }
    # print(eps)
  }
  return(list(x=x, iter=iter))
}

#' Preconditioned Conjugate Gradient for Two-level Toeplitz Left-hand Side
#'
#' Solves the second ADI equation, i.e. a linear system with two-level Toeplitz
#' left-hand side.
#'
#' @param M symbol of the circulant embedding of a two-level Toeplitz matrix,
#'          generated from the symbol the two-level Toeplitz matrix by
#'          \code{to_book_format()}
#' @param b right-hand side vector
#' @param init initial guess for solution
#' @param eigvals_prec eigenvalues of the circulant preconditioner generated by
#'                     \code{calc_prec()}
#' @param tol relative tolerance - stopping criterion
#'
#' @return the solution as a vector
#' @export
#'
#' @examples
#' eigvals <- matrix(9:1,3)
#' x <- stats::runif(9)
#' M <- to_book_format(eigvals,3,3)
#' eigvals_prec <- calc_prec(M)
#' pcgtt(M, x, runif(9), eigvals_prec)
pcgtt <- function(M, b, init, eigvals_prec, tol=1e-7){
  mmax = 2000; # the maximal number of iterations
  u <- init
  r <- b-BXfast(M,u)
  e <- rep(0,mmax)
  e[1] <- sqrt(sum(r^2))
  if(!is.finite(e[1])){
    e[1] <- 1
    iter <- mmax
    u <- rep(0,length(u))
  }
  iter = 1
  t1 = 1.0
  d = rep(0, length(init))
  while(iter < mmax && e[iter]/e[1] > tol){
    z = solve_prec(eigvals_prec,r)
    t1old = t1
    t1 = sum(z*r)
    beta = t1/t1old
    d = z+beta*d
    s = BXfast(M,d)
    suma = sum(d*s)
    tau = t1/suma
    u = u+tau*d
    r = r-tau*s
    iter = iter+1
    e[iter] = sqrt(sum(r^2))
    # cat('Step', iter, 'relative residual', e[iter]/e[1],'\n')
    if(!is.finite(e[iter])){
      e[iter] <- 0
      u <- rep(0,length(u))
    }
  }
  # if (iter == mmax) cat('Maximum iterations reached')
  # print(iter)
  return(list(u=u, iter=iter))
}

#' Calculates Eigenvalues of the Preconditioner
#'
#' Computes eigenvalues of the blocks of the circulant preconditioner, which is
#' denoted as \code{c^{(1)}} in Chan, R.H.F., & Jin, X.Q. (2007) An
#' introduction to iterative Toeplitz solvers.
#'
#' @param M the circulant symbol matrix generated by \code{to_book_format()}
#'
#' @return matrix of eigenvalues, one block of the preconditioner corresponds to
#'         one column of the eigenvalue matrix
#' @export
#'
#' @examples
#' eigvals <- matrix(9:1,3)
#' M <- to_book_format(eigvals,3,3)
#' calc_prec(M)
calc_prec <- function(M){
  n <- dim(M)[1]/2
  m <- dim(M)[2]
  X = array(0,c(n,m))
  X[1,] = M[1,];
  if(n>1) for (i in 2:n) X[i,] <- ((n-(i-1))*M[i,]+(i-1)*M[n+i,])/n;
  return(Re(stats::mvfft(X))) # FFT applied column by column
}

#' Solve the Preconditioning System
#'
#' Solves the linear system with the circulant preconditioner for the PCG, i.e.
#' it changes back the variables.
#'
#' @param eigvals eigenvalues of the preconditioner
#' @param x the right-hand side vector
#'
#' @return the solution as a matrix
#' @export
#'
#' @examples
#' eigvals <- matrix(9:1,3)
#' x <- stats::runif(9)
#' solve_prec(eigvals,x)
solve_prec <- function(eigvals, x){
  n <- dim(eigvals)[1]
  m <- dim(eigvals)[2]/2
  X = matrix(x,ncol=m)
  X = stats::mvfft(X)      # FFT column by column
  for(i in 1:n){
    A = stats::toeplitz(eigvals[i,1:m])
    X[i,] = t(solve(A, X[i,]))
    # We may solve the Toeplitz systems by the PCG methods or
    # by fast direct methods
  }
  X = stats::mvfft(X, inverse = T)/n # inverse FFT column by column
  return(c(Re(X)))
}
