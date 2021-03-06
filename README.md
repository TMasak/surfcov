
# surfcov

The purpose of R package `surfcov` is to enable covariance estimation
for random surfaces beyond separability, proposed in the papers
[arXiv:1912.12870](https://arxiv.org/abs/1912.12870) and
[arXiv:2007.12175](https://arxiv.org/abs/2007.12175).

Let
![X_1, \ldots, X_N](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;X_1%2C%20%5Cldots%2C%20X_N "X_1, \ldots, X_N")
be i.i.d. matrices of size
![K_1 \times K_2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;K_1%20%5Ctimes%20K_2 "K_1 \times K_2")
representing discrete measurements (on a grid) of some latent random
surfaces on a 2D domain, and
![C := cov(X_1)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;C%20%3A%3D%20cov%28X_1%29 "C := cov(X_1)")
be the covariance operator. The covariance is a tensor of size
![K_1 \times K_2 \times K_1 \times K_2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;K_1%20%5Ctimes%20K_2%20%5Ctimes%20K_1%20%5Ctimes%20K_2 "K_1 \times K_2 \times K_1 \times K_2"),
which becomes problematic to handle for
![K_1](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;K_1 "K_1")
and
![K_2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;K_2 "K_2")
as small as 100. The assumption of separability postulates that

![C\[i,j,i',j'\] = C_1\[i,i'\] C_2\[j,j'\],](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;C%5Bi%2Cj%2Ci%27%2Cj%27%5D%20%3D%20C_1%5Bi%2Ci%27%5D%20C_2%5Bj%2Cj%27%5D%2C "C[i,j,i',j'] = C_1[i,i'] C_2[j,j'],")

which reduces statistical and computational burden, but is often
critized as an oversimplification, since it does not allow any
interaction between the two dimensions.

This package allows for efficient estimation and subsequent manipulation
of two alternative models, which are both strict generalizations of
separability:

1.  the separable-plus-banded model:
    ![\[i,j,i',j'\] = A_1\[i,i'\] A_2\[j,j'\] + B\[i,j,i',j'\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Bi%2Cj%2Ci%27%2Cj%27%5D%20%3D%20A_1%5Bi%2Ci%27%5D%20A_2%5Bj%2Cj%27%5D%20%2B%20B%5Bi%2Cj%2Ci%27%2Cj%27%5D "[i,j,i',j'] = A_1[i,i'] A_2[j,j'] + B[i,j,i',j']"),
    where
    ![B\[i,j,i',j'\] = 0](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;B%5Bi%2Cj%2Ci%27%2Cj%27%5D%20%3D%200 "B[i,j,i',j'] = 0")
    for
    ![\|i-i'\| \> d](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%7Ci-i%27%7C%20%3E%20d "|i-i'| > d")
    or
    ![\|j-j'\| \> d](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%7Cj-j%27%7C%20%3E%20d "|j-j'| > d");
2.  the
    ![R](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;R "R")-separable
    model:
    ![C\[i,j,i',j'\] = A_1\[i,i'\]B_1\[j,j'\] + \ldots + A_R\[i,i'\]B_R\[j,j'\]](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;C%5Bi%2Cj%2Ci%27%2Cj%27%5D%20%3D%20A_1%5Bi%2Ci%27%5DB_1%5Bj%2Cj%27%5D%20%2B%20%5Cldots%20%2B%20A_R%5Bi%2Ci%27%5DB_R%5Bj%2Cj%27%5D "C[i,j,i',j'] = A_1[i,i']B_1[j,j'] + \ldots + A_R[i,i']B_R[j,j']").

When the data are stored in form of an array `X` of size
![N \times K_1 \times K_2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;N%20%5Ctimes%20K_1%20%5Ctimes%20K_2 "N \times K_1 \times K_2"),
it can be simply checked whether a given model can be potentially useful
by running

1.  `spb(X)` for the separable-plus-banded model; or
2.  `scd(X)` for the
    ![R](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;R "R")-separable
    model.

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
![N \times K_1 \times K_2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;N%20%5Ctimes%20K_1%20%5Ctimes%20K_2 "N \times K_1 \times K_2"),
it can be simply checked whether a given model can be potentially useful
by running

-   `spb(X)` for the separable-plus-banded model; or
-   `scd(X)` for the
    ![R](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;R "R")-separable
    model.

In particular,

``` r
# X <- array(runif(100*20*30),c(100,20,30))
spb_est <- spb(X)
spb$d
```

Run cross-validation (CV) in order to pick a value of the bandwidth
![d](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;d "d"),
and then fits the model with the best
![d](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;d "d")
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
![R](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;R "R")-separable
model, replacing function `spb` by `scd`, see `?spb` and `?scd` for a
more detailed description and examples. Note that by default, the mean
is always estimated empirically, unless other estimator of the mean is
provided.

Validity of a separable-plus-model can also be checked using a bootstrap
test, e.g.??by

``` r
test_spb(X,d=1)
```

When one of the models is fitted, a list of functions that can be useful
to the user for subsequent manipulation of the estimated covariances is
as follows:

-   `adi()` solves the inverse problem involving a separable-plus-banded
    covariance
-   `pcg()` solves the inverse problem involving an
    ![R](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;R "R")-separable
    covariance
-   `apply_lhs()` applies fast an
    ![R](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;R "R")-separable
    covariance

To apply the separable-plus-banded model fast, see `to_book_format` and
`BXfast`. TODO:

-   a handle for a fast application of the separable-plus-banded model
-   a handle for best linear unbiased predictor (BLUP) using the
    separable-plus-banded or
    ![R](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;R "R")-separable
    model.
