# Read data from bigBed files.

Columns are automatically typed based on the autoSql schema embedded in
the bigBed file. Integer types (`uint`, `int`) become R integers,
floating point types (`float`, `double`) become R doubles, and all other
types (including array types like `int[blockCount]`) remain as character
strings.

## Usage

``` r
read_bigbed(bbfile, chrom = NULL, start = NULL, end = NULL)
```

## Arguments

- bbfile:

  filename for bigBed file

- chrom:

  read data for specific chromosome

- start:

  start position for data

- end:

  end position for data

## Value

`tibble`

## See also

<https://github.com/dpryan79/libBigWig>

<https://github.com/brentp/bw-python>

## Examples

``` r
bb <- system.file("extdata", "test.bb", package = "cpp11bigwig")

read_bigbed(bb)
#> # A tibble: 3 × 12
#>   chrom  start    end name  score strand thickStart thickEnd reserved blockCount
#>   <chr>  <int>  <int> <chr> <int> <chr>       <int>    <int>    <int>      <int>
#> 1 chr1  4.80e6 4.84e6 test…     1 +         4797973  4836816        1          9
#> 2 chr10 4.85e6 4.88e6 diff…     1 +         4848118  4880877        1          6
#> 3 chr20 5.07e6 5.15e6 negs…     1 -         5073253  5152630        1         14
#> # ℹ 2 more variables: blockSizes <chr>, chromStarts <chr>

read_bigbed(bb, chrom = "chr10")
#> # A tibble: 1 × 12
#>   chrom  start    end name  score strand thickStart thickEnd reserved blockCount
#>   <chr>  <int>  <int> <chr> <int> <chr>       <int>    <int>    <int>      <int>
#> 1 chr10 4.85e6 4.88e6 diff…     1 +         4848118  4880877        1          6
#> # ℹ 2 more variables: blockSizes <chr>, chromStarts <chr>
```
