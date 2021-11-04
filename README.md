
<!-- README.md is generated from README.Rmd. Please edit that file -->

# surfcov

<!-- badges: start -->
<!-- badges: end -->

The purpose of R package `surfcov` is to enable covariance estimation
for random surfaces beyond separability.

Let X<sub>1</sub>,…,X<sub>N</sub> be i.i.d. matrices of size
(K<sub>1</sub> x K<sub>2</sub>) representing discrete measurements (on a
grid) on some latent random surfaces, and C:=cov(X<sub>1</sub>) is the
covariance, i.e. tensor of order 4. The assumption of separability
postulates

<p style="text-align: center;">
C
</p>

*C*\[*i*, *j*\]<sub>2</sub>⟨

## Installation

You can install the released version of surfcov from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("surfcov")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("TMasak/surfcov")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
# library(surfcov)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/master/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.
