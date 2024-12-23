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
#' bw <- system.file('extdata', 'test.bw', package = 'cpp11bigwig')
#'
#' read_bigwig(bw)
#'
#' read_bigwig(bw, chrom = "10")
#'
#' read_bigwig(bw, chrom = "1", start = 100, end = 130)
#'
#' read_bigwig(bw, as = 'gr')
#'
#' @export
read_bigwig <- function(bwfile, chrom = NULL, start = NULL, end = NULL, as = NULL) {

  if (!file.exists(bwfile)) {
    stop("File does not exist: ", bwfile)
  }

  if ((!is.null(start) && start < 0) || (!is.null(end) && end < 0)) {
    stop("`start` and `end` must both be >= 0")
  }

  if (!is.null(as) && !as %in% c("gr", "tbl")) {
    stop("`as` must be one of 'gr' or 'tbl' (the default)")
  }

  res <- read_bigwig_cpp(bwfile, chrom, start, end)

  if (!is.null(as) && as == 'gr') {
    return(as_gr(res))
  } else {
    return(as_tibble(res))
  }
}

#' convert to GRanges
#' @noRd
as_gr <- function(x) {
  GRanges(
    seqnames = x$chrom,
    ranges = IRanges(start = x$start, end = x$end),
    score = x$value
  )
}