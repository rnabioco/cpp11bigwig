
# cpp11bigwig

<!-- badges: start -->

[![R-CMD-check](https://github.com/rnabioco/cpp11bigwig/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/rnabioco/cpp11bigwig/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/rnabioco/cpp11bigwig/graph/badge.svg)](https://app.codecov.io/gh/rnabioco/cpp11bigwig)
<!-- badges: end -->

cpp11bigwig provides read access to bigWig files in R using `libBigWig`
from @dpryan79. Data is read into an R `data.frame`.

## Installation

<div class=".pkgdown-devel">

``` r
# Install development version from GitHub
# install.packages("pak")
pak::pak("rnabioco/cpp11bigwig")
```

</div>

``` r
library(cpp11bigwig)
```
