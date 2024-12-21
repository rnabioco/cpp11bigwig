#' Read data from bigWig files.
#'
#' @param bwfile filename or URL for bigWig file
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
#' bw <- system.file('extdata', 'test.bw', package = 'bigwrig')
#' read_bigwig(bw)
#'
#' @export
read_bigwig <- function(bwfile, genome, chrom = '', start = -1, end = -1) {
  read_bigwig_impl(bwfile, chrom, start, end)
}
