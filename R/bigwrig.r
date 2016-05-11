#' bigwrig: read data from bigWig files
#'
#' bigwrig provides methods to read data from bigWig files. bigWrig uses
#' Rcpp to wrap libBigWig from @@dpryan79.
#'
#' @author Jay Hesselberth <jay.hesselberth@gmail.com>
#'
#' @docType package
#' @name bigwrig
#'
#' @import dplyr
#' @useDynLib bigwrig
#' @exportPattern "^[[:alpha:]]+"
NULL

#' Read data from bigWig files.
#'
#' @param bwfile filename or URL for bigWig file
#' @param chrom read data for specific chromsome
#'
#' @return \code{data_frame}
#'
#' @seealso \url{https://github.com/dpryan79/libBigWig}
#' @seealso \url{https://github.com/brentp/bw-python}
#'
#' @examples
#'
#' @export
read_bigwig <- function() {}

#' @param bwfile filename or URL for bigWig file
#'
#' @return \code{data_frame}
#'
#' @export
read_bigwig_genome <- function() {}









