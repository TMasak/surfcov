
<!-- README.md is generated from README.Rmd. Please edit that file -->

# surfcov

<!-- badges: start -->
<!-- badges: end -->

The purpose of R package `surfcov` is to enable covariance estimation
for random surfaces beyond separability, proposed in the papers
[arXiv:1912.12870](https://arxiv.org/abs/1912.12870) and
[arXiv:2007.12175](https://arxiv.org/abs/2007.12175).

Let
*X*<sub>1</sub>, …, *X*<sub>*N*</sub>
be i.i.d. matrices of size
*K*<sub>1</sub> × *K*<sub>2</sub>
representing discrete measurements (on a grid) of some latent random
surfaces on a 2D domain, and
*C* := *c**o**v*(*X*<sub>1</sub>)
be the covariance operator. The covariance is a tensor of size
*K*<sub>1</sub> × *K*<sub>2</sub> × *K*<sub>1</sub> × *K*<sub>2</sub>
, which becomes problematic to handle for
*K*<sub>1</sub>
and
*K*<sub>2</sub>
as small as 100. The assumption of separability postulates that

*C*\[*i*, *j*, *i*′, *j*′\] = *C*<sub>1</sub>\[*i*, *i*′\]*C*<sub>2</sub>\[*j*, *j*′\],

which reduces statistical and computational burden, but is often
critized as an oversimplification, since it does not allow any
interaction between the two dimensions.

This package allows for efficient estimation and subsequent manipulation
of two alternative models, which are both strict generalizations of
separability:

1.  the separable-plus-banded model:
    *C*\[*i*, *j*, *i*′, *j*′\] = *A*<sub>1</sub>\[*i*, *i*′\]*A*<sub>2</sub>\[*j*, *j*′\] + *B*\[*i*, *j*, *i*′, *j*′\]
    , where
    *B*\[*i*, *j*, *i*′, *j*′\] = 0
    for
    \|*i* − *i*′\| &gt; *d*
    or
    \|*j* − *j*′\| &gt; *d*
    ;
2.  the
    *R*
    -separable model:
    *C*\[*i*, *j*, *i*′, *j*′\] = *A*<sub>1</sub>\[*i*, *i*′\]*B*<sub>1</sub>\[*j*, *j*′\] + … + *A*<sub>*R*</sub>\[*i*, *i*′\]*B*<sub>*R*</sub>\[*j*, *j*′\]
    .

When the data are stored in form of an array `X` of size
*N* × *K*<sub>1</sub> × *K*<sub>2</sub>
, it can be simply checked whether a given model can be potentially
useful by running

1.  `spb(X)` for the separable-plus-banded model; or
2.  `scd(X)` for the
    *R*
    -separable model.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
install.packages("devtools") # only if devtools is not yet installed
library(devtools)
install_github("TMasak/surfcov")
```

## Examples

When the data are stored in form of an array `X` of size
*N* × *K*<sub>1</sub> × *K*<sub>2</sub>
, it can be simply checked whether a given model can be potentially
useful by running

-   `spb(X)` for the separable-plus-banded model; or
-   `scd(X)` for the
    *R*
    -separable model.

In particular,

``` r
# X <- array(runif(100*20*30),c(100,20,30))
spb_est <- spb(X)
spb$d
```

Run cross-validation (CV) in order to pick a value of the bandwidth
*d*
, and then fits the model with the best
*d*
found. If a value larger than 0 is returned, it means that a
separable-plus-banded model can fit the data better compared to a
separable model. If we instead of fit-based CV decide to used
prediction-based CV, it makes sense to check the gains in prediction:

``` r
# X <- array(runif(100*20*30),c(100,20,30))
spb_est <- spb(X,predict=T)
spb$d
spb$cv
```

The same can be done with the
*R*
-separable model, replacing function `spb` by `scd`, see `?spb` and
`?scd` for a more detailed description and examples. Note that by
default, the mean is always estimated empirically, unless other
estimator of the mean is provided.

Validity of a separable-plus-model can also be checked using a bootstrap
test, e.g. by

``` r
test_spb(X,d=1)
```

When one of the models is fitted, a list of functions that can be useful
to the user for subsequent manipulation of the estimated covariances is
as follows:

-   `adi()` solves the inverse problem involving a separable-plus-banded
    covariance
-   `pcg()` solves the inverse problem involving an
    *R*
    -separable covariance
-   `apply_lhs()` applies fast an
    *R*
    -separable covariance

To apply the separable-plus-banded model fast, see `to_book_format` and
`BXfast`. TODO:

-   handle for a fast application of the separable-plus-banded model
-   handle for best linear unbiased predictor (BLUP) using the
    separable-plus-banded or
    *R*
    -separable model.