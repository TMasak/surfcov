% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SPT.R
\name{diagonal_band}
\alias{diagonal_band}
\title{Estimate the Entry-wise Variance}
\usage{
diagonal_band(X, A1, A2)
}
\arguments{
\item{X}{data array of size \code{N} x \code{K1} x \code{K2}}

\item{A1}{temporal kernel of the separable part}

\item{A2}{spatial kernel of the separable part}
}
\value{
matrix of size K1 x K2 providing the entry-wise variance
}
\description{
Estimates the variance of heteroscedastic white noise after the separable
part of the model has already been estimated.
}
\examples{
X <- array(runif(20*3*4),c(20,3,4))
A1 <- matrix(runif(3^2),3)
A2 <- matrix(runif(4^2),4)
diagonal_band(X,A1,A2)
}
