#' Read data from bigWig files.
#'
#' @param bwfile filename for bigWig file
#' @param chrom read data for specific chromosome
#' @param start start position for data
#' @param end end position for data
#'
#' @return \code{data.frame}
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
#' @export
read_bigwig <- function(bwfile, chrom = NULL, start = NULL, end = NULL) {

  if (!file.exists(bwfile)) {
    stop("File does not exist: ", bwfile)
  }

  if ((!is.null(start) && start < 0) || (!is.null(end) && end < 0)) {
    stop("`start` and `end` must both be >= 0")
  }

  read_bigwig_cpp(bwfile, chrom, start, end)
}
