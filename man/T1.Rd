% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PIP.R
\name{T1}
\alias{T1}
\title{Partial Inner Product w.r.t. the First Argument}
\usage{
T1(X, A, B, sigma, r, N, K1, K2)
}
\arguments{
\item{X}{data set, array of size \code{N} x \code{K1} x \code{K2}}

\item{A}{temporal kernels, array of size (\code{>= r}) x \code{K1} x
\code{K1}}

\item{B}{spatial kernels, array of size (\code{>= r}) x \code{K2} x \code{K2}}

\item{sigma}{separable component scores, vector of size \code{>= r}}

\item{r}{current sought separable rank}

\item{N}{sample size}

\item{K1}{temporal grid size}

\item{K2}{spatial grid size}
}
\value{
\code{r}-th temporal kernel A_r, array of size \code{K1} x \code{K1}
}
\description{
Calculates T_1(X x X, B_r), i.e. the partial inner product of the
covariance of X with the current weight B_r
resulting in a temporal kernel (proxy for A_r), which is then
orthogonalized w.r.t A_1,...,A_{r-1}.
}
\examples{
X <- array(rnorm(20*3*4),c(20,3,4))
A <- array(runif(2*3^2),c(2,3,3))
B <- array(runif(2*4^2),c(2,4,4))
T1(X, A, B, c(1,1), 2, 20, 3, 4)
}
