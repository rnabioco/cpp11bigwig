#' Read data from bigWig files.
#'
#' @param bwfile path or URL for a bigWig file. Remote files
#'  (`http://`, `https://`, `ftp://`) are supported when the package was
#'  installed with libcurl available.
#' @param chrom read data for specific chromosome
#' @param start start position for data
#' @param end end position for data
#' @param as return data as a specific type. One of `"tbl"` (the default
#'  tibble), `"GRanges"`, or `"Rle"`. `"Rle"` returns a per-base
#'  run-length-encoded vector spanning the requested range (see Details).
#' @param fill value used for bases with no data when `as = "Rle"`.
#'  Defaults to `0` (the convention for coverage); use `NA` to mark
#'  uncovered bases as missing. Ignored for other `as` values.
#'
#' @details
#' When `as = "Rle"`, the result is an [S4Vectors::Rle] whose expanded
#' length equals the queried range, i.e. `end - start` when both are
#' supplied, otherwise the extent of the returned data for each
#' chromosome. Bases with no data in the file are set to `fill`. bigWig
#' coordinates are 0-based and half-open, so element `i` corresponds to
#' genomic position `start + i - 1`. A single-chromosome query returns a
#' bare `Rle`; a multi-chromosome query returns a named
#' [IRanges::RleList].
#'
#' @return A `tibble`, `GRanges`, or `Rle`/`RleList` depending on `as`.
#'
#' @seealso \url{https://github.com/dpryan79/libBigWig}
#' @seealso \url{https://github.com/brentp/bw-python}
#'
#' @importFrom S4Vectors Rle
#' @importFrom IRanges RleList
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
#' read_bigwig(bw, chrom = "1", start = 100, end = 130, as = "Rle")
#'
#' @export
read_bigwig <- function(
  bwfile,
  chrom = NULL,
  start = NULL,
  end = NULL,
  as = NULL,
  fill = 0
) {
  check_bigwig_file(bwfile)

  if ((!is.null(start) && start < 0) || (!is.null(end) && end < 0)) {
    stop("`start` and `end` must both be >= 0")
  }

  if (!is.null(as) && !as %in% c("GRanges", "tbl", "Rle")) {
    stop("`as` must be one of 'GRanges', 'Rle', or 'tbl' (the default)")
  }

  res <- read_bigwig_cpp(bwfile, chrom, start, end)

  if (!is.null(as) && as == "GRanges") {
    return(as_granges(res))
  } else if (!is.null(as) && as == "Rle") {
    return(as_rle(res, start, end, fill))
  } else {
    return(as_tibble(res))
  }
}

#' is `x` a remote (http/https/ftp) URL?
#' @noRd
is_remote <- function(x) {
  grepl("^(https?|ftp)://", x, ignore.case = TRUE)
}

#' validate a bigWig/bigBed file path or URL
#' @noRd
check_bigwig_file <- function(file) {
  if (is_remote(file)) {
    if (!bigwig_has_curl_cpp()) {
      stop(
        "Remote files require libcurl, which was not available when ",
        "cpp11bigwig was installed. Reinstall with libcurl development ",
        "headers to enable remote access."
      )
    }
  } else if (!file.exists(file)) {
    stop("File does not exist: ", file)
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

#' build a per-base Rle covering [lo, hi) from sorted, clipped runs
#' @noRd
runs_to_rle <- function(s, e, v, lo, hi, fill) {
  # clip runs to the window and drop those entirely outside it
  keep <- e > lo & s < hi
  s <- pmax(s[keep], lo)
  e <- pmin(e[keep], hi)
  v <- v[keep]

  o <- order(s)
  s <- s[o]
  e <- e[o]
  v <- v[o]

  n <- length(s)
  if (n == 0L) {
    return(Rle(fill, max(hi - lo, 0L)))
  }

  # interleave a `fill` gap before each run with the run itself, then a
  # trailing `fill` gap out to `hi`; gaps of length 0 are dropped below
  prev_end <- c(lo, e[-n])
  values <- c(as.vector(rbind(rep(fill, n), v)), fill)
  lengths <- c(as.vector(rbind(s - prev_end, e - s)), hi - e[n])

  pos <- lengths > 0L
  Rle(values[pos], lengths[pos])
}

#' convert to an Rle (single chrom) or RleList (multiple chroms)
#' @noRd
as_rle <- function(x, start, end, fill) {
  chroms <- unique(x$chrom)

  # no data: with an explicit window, return a fill Rle of that length
  if (length(chroms) == 0L) {
    if (!is.null(start) && !is.null(end)) {
      return(Rle(fill, max(end - start, 0L)))
    }
    return(Rle(numeric(0)))
  }

  rles <- lapply(chroms, function(ch) {
    rows <- x$chrom == ch
    lo <- if (!is.null(start)) start else min(x$start[rows])
    hi <- if (!is.null(end)) end else max(x$end[rows])
    runs_to_rle(x$start[rows], x$end[rows], x$value[rows], lo, hi, fill)
  })

  if (length(rles) == 1L) {
    return(rles[[1]])
  }

  names(rles) <- chroms
  RleList(rles, compress = FALSE)
}

#' Read data from bigBed files.
#'
#' Columns are automatically typed based on the autoSql schema embedded
#' in the bigBed file. Integer types (`uint`, `int`) become R integers,
#' floating point types (`float`, `double`) become R doubles, and all
#' other types (including array types like `int[blockCount]`) remain
#' as character strings.
#'
#' @param bbfile path or URL for a bigBed file. Remote files
#'  (`http://`, `https://`, `ftp://`) are supported when the package was
#'  installed with libcurl available.
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
  check_bigwig_file(bbfile)

  if ((!is.null(start) && start < 0) || (!is.null(end) && end < 0)) {
    stop("`start` and `end` must both be >= 0")
  }

  res <- read_bigbed_cpp(bbfile, chrom, start, end)
  as_tibble(res)
}
