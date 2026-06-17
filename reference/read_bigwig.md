# Read data from bigWig files.

Read data from bigWig files.

## Usage

``` r
read_bigwig(
  bwfile,
  chrom = NULL,
  start = NULL,
  end = NULL,
  as = NULL,
  fill = 0
)
```

## Arguments

- bwfile:

  path or URL for a bigWig file. Remote files (`http://`, `https://`,
  `ftp://`) are supported when the package was installed with libcurl
  available.

- chrom:

  read data for specific chromosome

- start:

  start position for data

- end:

  end position for data

- as:

  return data as a specific type. One of `"tbl"` (the default tibble),
  `"GRanges"`, or `"Rle"`. `"Rle"` returns a per-base run-length-encoded
  vector spanning the requested range (see Details).

- fill:

  value used for bases with no data when `as = "Rle"`. Defaults to `0`
  (the convention for coverage); use `NA` to mark uncovered bases as
  missing. Ignored for other `as` values.

## Value

A `tibble`, `GRanges`, or `Rle`/`RleList` depending on `as`.

## Details

When `as = "Rle"`, the result is an
[S4Vectors::Rle](https://rdrr.io/pkg/S4Vectors/man/Rle-class.html) whose
expanded length equals the queried range, i.e. `end - start` when both
are supplied, otherwise the extent of the returned data for each
chromosome. Bases with no data in the file are set to `fill`. bigWig
coordinates are 0-based and half-open, so element `i` corresponds to
genomic position `start + i - 1`. A single-chromosome query returns a
bare `Rle`; a multi-chromosome query returns a named
[IRanges::RleList](https://rdrr.io/pkg/IRanges/man/AtomicList-class.html).

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

read_bigwig(bw, chrom = "1", start = 100, end = 130, as = "Rle")
#> numeric-Rle of length 30 with 1 run
#>   Lengths:  30
#>   Values : 1.4
```
