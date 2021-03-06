% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cvband_nonstat.R
\name{cvband_nonstat}
\alias{cvband_nonstat}
\title{Cross-validation for Bandwidth Selection in Separable-plus-banded Model
without Stationarity}
\usage{
cvband_nonstat(X, Folds = 10, maxd = 2, mind = 0)
}
\arguments{
\item{X}{data set, array of size \code{N} x \code{K1} x \code{K2}}

\item{Folds}{number of folds}

\item{maxd}{maximum bandwidth value to check}

\item{mind}{minimum bandwidth value to check}
}
\value{
matrix with two rows: first row gives the fit, the second gives the
norm, columns correspond to different bandwidths
}
\description{
Similar to \code{cvband}, but the banded part of the model is not assumed to
be stationary here. Because of that this function is computationally more
demanding compared to \code{cvband}. The purpose of this function
is to check whether there might be heteroscedastic noise in the data,
corresponding to a separable-plus-diagonal model, i.e. to \code{d=1}. For
that, one needs to look past the separable-plus-diagonal, i.e. search for
\code{d=0,1,..,maxd}. However, relatively small \code{maxd} should be enough.
}
\examples{
N <- 20
K1 <- 5
K2 <- 7
set.seed(517)
X <- array(0,c(N,K1,K2))
A <- matrix(rnorm(K1^2),K1)
B <- matrix(rnorm(K2^2),K2)
A <- A \%*\% t(A)
B <- B \%*\% t(B)
A <- mat_root(A)
B <- mat_root(B)
for(n in 1:N){
  X[n,,] <- A \%*\% X[n,,] \%*\% B + 2*matrix(rnorm(K1*K2),K1)
}
cvscores <- cvband_nonstat(X,10,4,0) # two rows: fit and norm
cvscores[1,]/cvscores[2,] # cross-validation objective: fit/norm
min(localMaxima(cvscores[1,]/cvscores[2,])) # choose bandwidth as the
                                            # smallest local maximum
}
