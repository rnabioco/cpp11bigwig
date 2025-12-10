#' Read data from bigWig files.
#'
#' @param bwfile filename for bigWig file
#' @param chrom read data for specific chromosome
#' @param start start position for data
#' @param end end position for data
#' @param as return data as a specific type.
#'  The default is a tibble (`tbl`) or GRanges (`gr`)
#'
#' @return \code{tibble}
#'
#' @seealso \url{https://github.com/dpryan79/libBigWig}
#' @seealso \url{https://github.com/brentp/bw-python}
#'
#' @examples
#' bw <- system.file("extdata", "test.bw", package = "cpp11bigwig")
#'
#' read_bigwig(bw)
#'
#' read_bigwig(bw, chrom = "10")
#'
#' read_bigwig(bw, chrom = "1", start = 100, end = 130)
#'
#' read_bigwig(bw, as = "GRanges")
#'
#' @export
read_bigwig <- function(
  bwfile,
  chrom = NULL,
  start = NULL,
  end = NULL,
  as = NULL
) {
  if (!file.exists(bwfile)) {
    stop("File does not exist: ", bwfile)
  }

  if ((!is.null(start) && start < 0) || (!is.null(end) && end < 0)) {
    stop("`start` and `end` must both be >= 0")
  }

  if (!is.null(as) && !as %in% c("GRanges", "tbl")) {
    stop("`as` must be one of 'GRanges' or 'tbl' (the default)")
  }

  res <- read_bigwig_cpp(bwfile, chrom, start, end)

  if (!is.null(as) && as == "GRanges") {
    return(as_granges(res))
  } else {
    return(as_tibble(res))
  }
}

#' convert to GRanges
#' @noRd
as_granges <- function(x) {
  GRanges(
    seqnames = x$chrom,
    ranges = IRanges(start = x$start, end = x$end),
    score = x$value
  )
}

#' Read data from bigBed files.
#'
#' Columns are automatically typed based on the autoSql schema embedded
#' in the bigBed file. Integer types (`uint`, `int`) become R integers,
#' floating point types (`float`, `double`) become R doubles, and all
#' other types (including array types like `int[blockCount]`) remain
#' as character strings.
#'
#' @param bbfile filename for bigBed file
#' @param chrom read data for specific chromosome
#' @param start start position for data
#' @param end end position for data
#'
#' @return \code{tibble}
#'
#' @seealso \url{https://github.com/dpryan79/libBigWig}
#' @seealso \url{https://github.com/brentp/bw-python}
#'
#' @examples
#' bb <- system.file("extdata", "test.bb", package = "cpp11bigwig")
#'
#' read_bigbed(bb)
#'
#' read_bigbed(bb, chrom = "chr10")
#'
#' @export
read_bigbed <- function(
  bbfile,
  chrom = NULL,
  start = NULL,
  end = NULL
) {
  if (!file.exists(bbfile)) {
    stop("File does not exist: ", bbfile)
  }

  if ((!is.null(start) && start < 0) || (!is.null(end) && end < 0)) {
    stop("`start` and `end` must both be >= 0")
  }

  res <- read_bigbed_cpp(bbfile, chrom, start, end)
  as_tibble(res)
}
