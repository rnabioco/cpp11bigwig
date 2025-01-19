
# cpp11bigwig

<!-- badges: start -->

[![R-CMD-check](https://github.com/rnabioco/cpp11bigwig/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/rnabioco/cpp11bigwig/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/rnabioco/cpp11bigwig/graph/badge.svg)](https://app.codecov.io/gh/rnabioco/cpp11bigwig)
<!-- badges: end -->

cpp11bigwig provides read-only access to bigWig and bigBed files using
libBigWig <https://github.com/dpryan79/libBigWig>.

## Installation

<div class=".pkgdown-release">

``` r
# Install released version from CRAN
install.packages("cpp11bigwig")
```

</div>

<div class=".pkgdown-devel">

``` r
# Install development version from GitHub
# install.packages("pak")
pak::pak("rnabioco/cpp11bigwig")
```

</div>

## Usage

``` r
library(cpp11bigwig)

bw <- system.file("extdata", "test.bw", package = "cpp11bigwig")
read_bigwig(bw)
#> # A tibble: 6 × 4
#>   chrom start   end value
#>   <chr> <int> <int> <dbl>
#> 1 1         0     1 0.100
#> 2 1         1     2 0.200
#> 3 1         2     3 0.300
#> 4 1       100   150 1.40 
#> 5 1       150   151 1.5  
#> 6 10      200   300 2

bb <- system.file("extdata", "test.bb", package = "cpp11bigwig")
read_bigbed(bb)
#> # A tibble: 3 × 12
#>   chrom  start    end name  score strand thickStart thickEnd reserved blockCount
#>   <chr>  <int>  <int> <chr> <chr> <chr>  <chr>      <chr>    <chr>    <chr>     
#> 1 chr1  4.80e6 4.80e6 test… 1     +      4797973    4836816  1        9         
#> 2 chr10 4.85e6 4.85e6 diff… 1     +      4848118    4880877  1        6         
#> 3 chr20 5.07e6 5.07e6 negs… 1     -      5073253    5152630  1        14        
#> # ℹ 2 more variables: blockSizes <chr>, chromStarts <chr>
```
