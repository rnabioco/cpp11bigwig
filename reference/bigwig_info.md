# Report header metadata and summary statistics for a bigWig file.

Reads the bigWig header without loading any intervals. The summary
statistics (`min`, `max`, `mean`, `std`) are the file-level values
stored in the header and computed over all covered bases.

## Usage

``` r
bigwig_info(bwfile)
```

## Arguments

- bwfile:

  path or URL for a bigWig file. Remote files (`http://`, `https://`,
  `ftp://`) are supported when the package was installed with libcurl
  available.

## Value

A named list with elements `version`, `n_levels`, `n_chroms`,
`n_bases_covered`, `min`, `max`, `mean`, and `std`.

## See also

[`read_bigwig()`](https://rnabioco.github.io/cpp11bigwig/reference/read_bigwig.md),
[`bigbed_info()`](https://rnabioco.github.io/cpp11bigwig/reference/bigbed_info.md)

## Examples

``` r
bw <- system.file("extdata", "test.bw", package = "cpp11bigwig")

bigwig_info(bw)
#> $version
#> [1] 4
#> 
#> $n_levels
#> [1] 1
#> 
#> $n_chroms
#> [1] 2
#> 
#> $n_bases_covered
#> [1] 154
#> 
#> $min
#> [1] 0.1
#> 
#> $max
#> [1] 2
#> 
#> $mean
#> [1] 1.766883
#> 
#> $std
#> [1] 0.356945
#> 
```
