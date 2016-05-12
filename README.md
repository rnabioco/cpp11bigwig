`bigwrig`: access bigWig files with R
================

`bigwrig` provides read access to bigWig files in R using `libBigWig` from @dpryan79. Data is read into an R `data_frame`.

``` r
library(bigwrig)

url <- ''

read_bigwig(url)

read_bigwig_genome(url)
```
