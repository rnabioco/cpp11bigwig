# Report header metadata for a bigBed file.

Reads the bigBed header without loading any intervals. This is useful
for identifying the BED variant a file holds before reading it: a
genuine BED12 has `defined_field_count == 12`, whereas a `bed9+3` file
(9 standard BED columns plus 3 custom fields) has
`defined_field_count == 9` and `field_count == 12`.

## Usage

``` r
bigbed_info(bbfile)
```

## Arguments

- bbfile:

  path or URL for a bigBed file. Remote files (`http://`, `https://`,
  `ftp://`) are supported when the package was installed with libcurl
  available.

## Value

A named list with elements `version`, `n_chroms`, `field_count`,
`defined_field_count`, `n_bases_covered`, and `autosql` (the embedded
autoSql schema string, or `""` when the file has none).

## See also

[`read_bigbed()`](https://rnabioco.github.io/cpp11bigwig/reference/read_bigbed.md),
[`bigwig_info()`](https://rnabioco.github.io/cpp11bigwig/reference/bigwig_info.md)

## Examples

``` r
bb <- system.file("extdata", "test.bb", package = "cpp11bigwig")

info <- bigbed_info(bb)
info$defined_field_count
#> [1] 12
```
