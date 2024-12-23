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

#' Read data from bigBed files.
#'
#' @param bbfile filename for bigBed file
#' @param chrom read data for specific chromosome
#' @param start start position for data
#' @param end end position for data
#' @param convert convert bigBed values to individual columns
#'
#' @return \code{tibble}
#'
#' @seealso \url{https://github.com/dpryan79/libBigWig}
#' @seealso \url{https://github.com/brentp/bw-python}
#'
#' @examples
#' bb <- system.file('extdata', 'test.bb', package = 'cpp11bigwig')
#'
#' read_bigbed(bb)
#'
#' read_bigbed(bb, chrom = "chr10")
#'
#' @export
read_bigbed <- function(bbfile, chrom = NULL, start = NULL, end = NULL, convert = TRUE) {

  if (!file.exists(bbfile)) {
    stop("File does not exist: ", bbfile)
  }

  if ((!is.null(start) && start < 0) || (!is.null(end) && end < 0)) {
    stop("`start` and `end` must both be >= 0")
  }

  res <- read_bigbed_cpp(bbfile, chrom, start, end)

  if (!convert) {
    return(as_tibble(res))
  }

  vals <- do.call(rbind, strsplit(res[["value"]], "\t"))
  # merge chrom, start, end with new values
  res_new <- cbind(res[, 1:3], vals)

  fnames <- bigbed_sql_fields(bbfile)
  # drop chrom, start, end
  fnames <- fnames[-(1:3)]

  colnames(res_new) <- c("chrom", "start", "end", fnames)
  return(as_tibble(res_new))
}

#' @examples
#' bb <- system.file('extdata', 'test.bb', package = 'cpp11bigwig')
#' bigbed_sql_fields(bb)
#'
#' @noRd
bigbed_sql_fields <- function(bbfile) {

  res <- bigbed_sql_cpp(bbfile)

  # parse the autoSql
  lines <- unlist(strsplit(res, "\n"))
  fields <- lines[grep(";", lines)]

  unlist(
    lapply(
      fields, function(line) {
        field <- sub("^\\s*\\S+\\s+(\\S+);.*", "\\1", line)
        field
      }
    )
  )
}
