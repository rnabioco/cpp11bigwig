#' Read data from bigWig files.
#'
#' @param bwfile path or URL for a bigWig file. Remote files
#'  (`http://`, `https://`, `ftp://`) are supported when the package was
#'  installed with libcurl available.
#' @param chrom chromosome(s) to read. Either a character vector of
#'  chromosome names, or a [GenomicRanges::GRanges] of query regions (in
#'  which case `start`/`end` are ignored; see Details).
#' @param start start position(s) for data. May be a vector, recycled
#'  against `chrom`/`end` to describe several ranges.
#' @param end end position(s) for data. May be a vector, recycled
#'  against `chrom`/`start` to describe several ranges.
#' @param as return data as a specific type. One of `"tbl"` (the default
#'  tibble), `"GRanges"`, or `"Rle"`. `"Rle"` returns a per-base
#'  run-length-encoded vector spanning the requested range (see Details).
#' @param fill value used for bases with no data when `as = "Rle"`.
#'  Defaults to `0` (the convention for coverage); use `NA` to mark
#'  uncovered bases as missing. Ignored for other `as` values.
#'
#' @details
#' Multiple ranges can be queried in one call by passing equal-length
#' (or length-1, recycled) `chrom`, `start`, and `end` vectors, where
#' range `i` is `(chrom[i], start[i], end[i])`. Alternatively, pass a
#' [GenomicRanges::GRanges] as `chrom`; its regions are used directly.
#' Because `GRanges` is 1-based and inclusive while bigWig is 0-based and
#' half-open, a region is converted as `start(gr) - 1` to `end(gr)`.
#'
#' When `as = "Rle"`, the result is an [S4Vectors::Rle] whose expanded
#' length equals the queried range, i.e. `end - start` when both are
#' supplied, otherwise the extent of the returned data for each
#' chromosome. Bases with no data in the file are set to `fill`. bigWig
#' coordinates are 0-based and half-open, so element `i` corresponds to
#' genomic position `start + i - 1`. A single-range query returns a bare
#' `Rle`; a multi-range (or multi-chromosome) query returns a named
#' [IRanges::RleList] with one element per range.
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
#' # query several ranges in one call with equal-length vectors
#' read_bigwig(bw, chrom = c("1", "10"), start = c(0, 0), end = c(50, 50))
#'
#' # multiple windows on the same chromosome (chrom recycles)
#' read_bigwig(bw, chrom = "1", start = c(0, 100), end = c(50, 130))
#'
#' # a multi-range "Rle" query returns a named RleList, one element per range
#' read_bigwig(bw, chrom = "1", start = c(0, 100), end = c(50, 130), as = "Rle")
#'
#' # pass a GRanges of regions; 1-based coords are converted automatically
#' gr <- GenomicRanges::GRanges(
#'   c("1", "10"),
#'   IRanges::IRanges(start = c(1, 1), end = c(50, 50))
#' )
#' read_bigwig(bw, chrom = gr)
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

  if (!is.null(as) && !as %in% c("GRanges", "tbl", "Rle")) {
    stop("`as` must be one of 'GRanges', 'Rle', or 'tbl' (the default)")
  }

  ranges <- normalize_ranges(chrom, start, end)
  multi <- !is.null(ranges)
  if (!multi) {
    ranges <- single_range(chrom, start, end)
  }
  check_ranges_nonneg(ranges)

  # one C++ call opens the file once and returns a per-range list of frames
  reslist <- read_ranges(read_bigwig_cpp, bwfile, ranges)

  if (!is.null(as) && as == "Rle") {
    # single range / whole file: split by chromosome (bare Rle, or an
    # RleList named by chromosome when several chroms are present)
    if (!multi) {
      return(as_rle(reslist[[1]], start, end, fill))
    }
    # multiple ranges: one Rle per requested range
    rles <- lapply(seq_len(nrow(ranges)), function(i) {
      x <- reslist[[i]]
      lo <- if (!is.na(ranges$start[i])) {
        ranges$start[i]
      } else if (length(x$start)) {
        min(x$start)
      } else {
        0L
      }
      hi <- if (!is.na(ranges$end[i])) {
        ranges$end[i]
      } else if (length(x$end)) {
        max(x$end)
      } else {
        lo
      }
      runs_to_rle(x$start, x$end, x$value, lo, hi, fill)
    })
    # a single requested range collapses to a bare Rle
    if (length(rles) == 1L) {
      return(rles[[1]])
    }
    names(rles) <- sprintf(
      "%s:%s-%s",
      ranges$chrom,
      format(ranges$start, trim = TRUE),
      format(ranges$end, trim = TRUE)
    )
    return(RleList(rles, compress = FALSE))
  }

  combined <- do.call(rbind, reslist)
  if (!is.null(as) && as == "GRanges") {
    return(as_granges(combined))
  }
  as_tibble(combined)
}

#' normalize range inputs into a data.frame of (chrom, start, end) rows,
#' or NULL for the single-range / whole-file legacy path.
#'
#' `chrom` may be a GRanges (start/end taken from it, ignoring the
#' `start`/`end` args) or a character vector recycled against `start`/`end`.
#' Unspecified start/end become NA (passed to C++ as NULL per range).
#' @noRd
normalize_ranges <- function(chrom, start, end) {
  if (inherits(chrom, "GRanges")) {
    df <- as.data.frame(chrom)
    return(data.frame(
      chrom = as.character(df$seqnames),
      start = df$start - 1L, # 1-based inclusive -> 0-based half-open
      end = df$end,
      stringsAsFactors = FALSE
    ))
  }

  n <- max(length(chrom), length(start), length(end))
  if (n <= 1L) {
    return(NULL)
  }

  recycle <- function(x) {
    if (is.null(x)) {
      return(rep(NA, n))
    }
    if (length(x) == 1L) {
      return(rep(x, n))
    }
    if (length(x) != n) {
      stop(
        "`chrom`, `start`, and `end` must have the same length (or length 1)"
      )
    }
    x
  }

  data.frame(
    chrom = as.character(recycle(chrom)),
    start = as.numeric(recycle(start)),
    end = as.numeric(recycle(end)),
    stringsAsFactors = FALSE
  )
}

#' error if any range has a negative start or end
#' @noRd
check_ranges_nonneg <- function(ranges) {
  bad_start <- !is.na(ranges$start) & ranges$start < 0
  bad_end <- !is.na(ranges$end) & ranges$end < 0
  if (any(bad_start) || any(bad_end)) {
    stop("`start` and `end` must both be >= 0")
  }
}

#' build the single-row range frame for the whole-file / single-query path,
#' mapping unspecified (NULL) chrom/start/end to NA (every chromosome / the
#' chromosome bounds in C++).
#' @noRd
single_range <- function(chrom, start, end) {
  data.frame(
    chrom = if (is.null(chrom)) NA_character_ else as.character(chrom),
    start = if (is.null(start)) NA_real_ else as.numeric(start),
    end = if (is.null(end)) NA_real_ else as.numeric(end),
    stringsAsFactors = FALSE
  )
}

#' call a vectorized C++ reader once for all ranges; NA chrom selects every
#' chromosome and NA start/end default to the chromosome bounds. Returns a
#' list with one data frame per range.
#' @noRd
read_ranges <- function(fun, file, ranges) {
  fun(
    file,
    as.character(ranges$chrom),
    as.integer(ranges$start),
    as.integer(ranges$end)
  )
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
#' @param chrom chromosome(s) to read. Either a character vector of
#'  chromosome names, or a [GenomicRanges::GRanges] of query regions (in
#'  which case `start`/`end` are ignored). As with [read_bigwig()],
#'  `GRanges` 1-based coordinates are converted to bigBed's 0-based
#'  half-open coordinates.
#' @param start start position(s) for data. May be a vector describing
#'  several ranges, recycled against `chrom`/`end`.
#' @param end end position(s) for data. May be a vector describing
#'  several ranges, recycled against `chrom`/`start`.
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
#' # query several chromosomes in one call
#' read_bigbed(bb, chrom = c("chr1", "chr10"))
#'
#' # restrict each query to a window
#' read_bigbed(bb, chrom = c("chr1", "chr10"), start = c(0, 0), end = c(5e6, 5e6))
#'
#' # pass a GRanges of regions; 1-based coords are converted automatically
#' gr <- GenomicRanges::GRanges(
#'   c("chr1", "chr10"),
#'   IRanges::IRanges(start = 1, width = 1e7)
#' )
#' read_bigbed(bb, chrom = gr)
#'
#' @export
read_bigbed <- function(
  bbfile,
  chrom = NULL,
  start = NULL,
  end = NULL
) {
  check_bigwig_file(bbfile)

  ranges <- normalize_ranges(chrom, start, end)
  if (is.null(ranges)) {
    ranges <- single_range(chrom, start, end)
  }
  check_ranges_nonneg(ranges)

  # one C++ call opens the file once and returns a per-range list of frames
  reslist <- read_ranges(read_bigbed_cpp, bbfile, ranges)
  as_tibble(do.call(rbind, reslist))
}
