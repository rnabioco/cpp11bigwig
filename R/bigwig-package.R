#' cpp11bigwig: read data from bigWig files
#'
#' bigwig provides methods to read data from bigWig files. bigwig uses
#' cpp11 to wrap libBigWig from @@dpryan79.
#'
#' <https://github.com/dpryan79/libBigWig>
#'
#' @author Jay Hesselberth <jay.hesselberth@gmail.com>
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @useDynLib cpp11bigwig, .registration = TRUE
#' @importFrom tibble as_tibble
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
## usethis namespace: end
NULL
