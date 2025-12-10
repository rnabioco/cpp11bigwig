# Read data from bigWig files.

Read data from bigWig files.

## Usage

``` r
read_bigwig(bwfile, chrom = NULL, start = NULL, end = NULL, as = NULL)
```

## Arguments

- bwfile:

  filename for bigWig file

- chrom:

  read data for specific chromosome

- start:

  start position for data

- end:

  end position for data

- as:

  return data as a specific type. The default is a tibble (`tbl`) or
  GRanges (`gr`)

## Value

`tibble`

## See also

<https://github.com/dpryan79/libBigWig>

<https://github.com/brentp/bw-python>

## Examples

``` r
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

read_bigwig(bw, chrom = "10")
#> # A tibble: 1 × 4
#>   chrom start   end value
#>   <chr> <int> <int> <dbl>
#> 1 10      200   300     2

read_bigwig(bw, chrom = "1", start = 100, end = 130)
#> # A tibble: 1 × 4
#>   chrom start   end value
#>   <chr> <int> <int> <dbl>
#> 1 1       100   130  1.40

read_bigwig(bw, as = "GRanges")
#> GRanges object with 6 ranges and 1 metadata column:
#>       seqnames    ranges strand |     score
#>          <Rle> <IRanges>  <Rle> | <numeric>
#>   [1]        1       0-1      * |       0.1
#>   [2]        1       1-2      * |       0.2
#>   [3]        1       2-3      * |       0.3
#>   [4]        1   100-150      * |       1.4
#>   [5]        1   150-151      * |       1.5
#>   [6]       10   200-300      * |       2.0
#>   -------
#>   seqinfo: 2 sequences from an unspecified genome; no seqlengths
```
