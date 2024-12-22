
# cpp11bigwig

<!-- badges: start -->

[![R-CMD-check](https://github.com/rnabioco/cpp11bigwig/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/rnabioco/cpp11bigwig/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/rnabioco/cpp11bigwig/graph/badge.svg)](https://app.codecov.io/gh/rnabioco/cpp11bigwig)
<!-- badges: end -->

cpp11bigwig provides read access to bigWig files in R using
[libBigWig](https://github.com/dpryan79/libBigWig). Data is read into an
R `data.frame`.

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

bw = system.file('extdata', 'test.bw', package = 'cpp11bigwig')

read_bigwig(bw)
```

    ##   chrom start end value
    ## 1     1     0   1   0.1
    ## 2     1     1   2   0.2
    ## 3     1     2   3   0.3
    ## 4     1   100 150   1.4
    ## 5     1   150 151   1.5
    ## 6    10   200 300   2.0
