test_that("results have expected shape", {
  bw <- test_path("data/test.bw")

  res <- read_bigwig(bw)
  expect_equal(ncol(res), 4)
  expect_equal(nrow(res), 6)
  expect_true("tbl_df" %in% class(res))

  # interval sizes
  expect_true(all(res$end - res$start == c(1, 1, 1, 50, 1, 100)))
  # interval values
  expect_equal(sum(res$value), 5.5)

  # params work
  expect_equal(nrow(read_bigwig(bw, chrom = "1")), 5)
  expect_equal(nrow(read_bigwig(bw, start = 100)), 3)
  expect_equal(nrow(read_bigwig(bw, end = 3)), 3)

  # GRanges
  res <- read_bigwig(bw, as = "GRanges")
  expect_true("GRanges" %in% class(res))

  # bigbed
  bb <- test_path("data/test.bb")
  res <- read_bigbed(bb)

  expect_equal(ncol(res), 12)
  expect_equal(nrow(res), 3)
})

test_that("as = 'Rle' returns a per-base run-length vector", {
  bw <- test_path("data/test.bw")

  # values are stored as 32-bit floats, so compare with tolerance
  tol <- 1e-6

  # windowed query: length equals end - start, gaps filled with 0
  r <- read_bigwig(bw, chrom = "1", start = 0, end = 5, as = "Rle")
  expect_s4_class(r, "Rle")
  expect_equal(length(r), 5L)
  expect_equal(as.numeric(r), c(0.1, 0.2, 0.3, 0, 0), tolerance = tol)

  # uncovered bases inside the data range are also filled
  r <- read_bigwig(bw, chrom = "1", start = 98, end = 103, as = "Rle")
  expect_equal(as.numeric(r), c(0, 0, 1.4, 1.4, 1.4), tolerance = tol)

  # fill = NA marks uncovered bases as missing
  r <- read_bigwig(bw, chrom = "1", start = 0, end = 5, as = "Rle", fill = NA)
  expect_equal(as.numeric(r), c(0.1, 0.2, 0.3, NA, NA), tolerance = tol)

  # without an explicit window the Rle spans the data extent of the chrom
  r <- read_bigwig(bw, chrom = "1", as = "Rle")
  expect_equal(length(r), 151L)

  # multiple chromosomes return a named RleList
  rl <- read_bigwig(bw, as = "Rle")
  expect_s4_class(rl, "RleList")
  expect_equal(names(rl), c("1", "10"))
})

test_that("read_bigwig accepts multiple ranges", {
  bw <- test_path("data/test.bw")

  # vectorized chrom/start/end: result is the rbind of the single queries
  multi <- read_bigwig(bw, chrom = c("1", "10"), start = c(0, 0), end = c(50, 50))
  s1 <- read_bigwig(bw, chrom = "1", start = 0, end = 50)
  s2 <- read_bigwig(bw, chrom = "10", start = 0, end = 50)
  expect_equal(nrow(multi), nrow(s1) + nrow(s2))
  expect_equal(
    as.data.frame(multi),
    as.data.frame(rbind(s1, s2))
  )

  # two ranges on the same chromosome are both returned
  same <- read_bigwig(bw, chrom = "1", start = c(0, 100), end = c(50, 150))
  expect_equal(unique(same$chrom), "1")
  expect_gt(nrow(same), nrow(s1))

  # scalar chrom recycles against vector start/end
  expect_equal(
    nrow(read_bigwig(bw, chrom = "1", start = c(0, 100), end = c(50, 150))),
    nrow(same)
  )

  # mismatched lengths are an error
  expect_error(
    read_bigwig(bw, chrom = c("1", "10"), start = c(0, 0, 0), end = 50),
    "same length"
  )

  # negative coordinates still error in multi-range mode
  expect_error(
    read_bigwig(bw, chrom = c("1", "10"), start = c(-1, 0), end = c(5, 5)),
    ">= 0"
  )
})

test_that("read_bigwig accepts a GRanges of regions", {
  skip_if_not_installed("GenomicRanges")
  bw <- test_path("data/test.bw")

  # GRanges is 1-based inclusive; bigWig is 0-based half-open. A GRanges
  # region start(gr):end(gr) must match the vectorized (start(gr) - 1, end(gr)).
  gr <- GenomicRanges::GRanges(
    c("1", "10"),
    IRanges::IRanges(start = c(1, 1), end = c(50, 50))
  )
  from_gr <- read_bigwig(bw, chrom = gr)
  from_vec <- read_bigwig(bw, chrom = c("1", "10"), start = c(0, 0), end = c(50, 50))
  expect_equal(as.data.frame(from_gr), as.data.frame(from_vec))

  # GRanges combines with as = "GRanges" too
  g <- read_bigwig(bw, chrom = gr, as = "GRanges")
  expect_s4_class(g, "GRanges")
})

test_that("multi-range as = 'Rle' returns a per-range RleList", {
  bw <- test_path("data/test.bw")

  rl <- read_bigwig(
    bw,
    chrom = c("1", "10"),
    start = c(0, 0),
    end = c(50, 50),
    as = "Rle"
  )
  expect_s4_class(rl, "RleList")
  expect_equal(length(rl), 2L)
  expect_equal(names(rl), c("1:0-50", "10:0-50"))
  # each element spans its requested window
  expect_equal(lengths(rl), c("1:0-50" = 50L, "10:0-50" = 50L))

  # a single requested range still collapses to a bare Rle
  r1 <- read_bigwig(bw, chrom = "1", start = 0, end = 50, as = "Rle")
  expect_s4_class(r1, "Rle")
})

test_that("read_bigbed accepts multiple ranges", {
  bb <- test_path("data/test.bb")

  multi <- read_bigbed(
    bb,
    chrom = c("chr1", "chr10"),
    start = c(0, 0),
    end = c(1e7, 1e7)
  )
  s1 <- read_bigbed(bb, chrom = "chr1", start = 0, end = 1e7)
  s2 <- read_bigbed(bb, chrom = "chr10", start = 0, end = 1e7)
  expect_s3_class(multi, "tbl_df")
  expect_equal(nrow(multi), nrow(s1) + nrow(s2))

  # GRanges input works for bigBed as well
  skip_if_not_installed("GenomicRanges")
  gr <- GenomicRanges::GRanges(
    c("chr1", "chr10"),
    IRanges::IRanges(start = c(1, 1), end = c(1e7, 1e7))
  )
  from_gr <- read_bigbed(bb, chrom = gr)
  expect_equal(nrow(from_gr), nrow(multi))
})

test_that("missing file causes error", {
  expect_snapshot_error(read_bigwig("missing.bw"))
})

test_that("negative coords causes error", {
  expect_snapshot_error(read_bigwig(test_path("data/test.bw"), start = -1))
  expect_snapshot_error(read_bigwig(test_path("data/test.bw"), end = -1))
})

test_that("is_remote detects URLs", {
  expect_true(is_remote("https://example.com/x.bw"))
  expect_true(is_remote("http://example.com/x.bw"))
  expect_true(is_remote("ftp://example.com/x.bw"))
  expect_true(is_remote("HTTPS://EXAMPLE.COM/x.bw"))
  expect_false(is_remote("/local/path/x.bw"))
  expect_false(is_remote("x.bw"))
})

test_that("remote bigWig reads match local", {
  skip_on_cran()
  skip_if_not(bigwig_has_curl_cpp(), "built without libcurl")

  url <- paste0(
    "https://raw.githubusercontent.com/rnabioco/cpp11bigwig/",
    "main/inst/extdata/test.bw"
  )
  local <- read_bigwig(system.file("extdata", "test.bw", package = "cpp11bigwig"))
  # a network/remote-host failure should skip, not fail the suite
  remote <- tryCatch(
    read_bigwig(url),
    error = function(e) skip(paste("remote fetch unavailable:", conditionMessage(e)))
  )

  expect_equal(remote, local)
})

test_that("remote bigWig larger than the read buffer reads a windowed range (#18)", {
  # Regression test: the HTTP Range header was not being set, so servers
  # returned the entire file. Files larger than the internal read buffer
  # (1 << 17 bytes) overran it and crashed R or failed to open. This reads a
  # small window from an 81 MB file, which only succeeds with working range
  # requests (it must not download the whole file).
  skip_on_cran()
  skip_if_not(bigwig_has_curl_cpp(), "built without libcurl")

  url <- "https://genome.ucsc.edu/goldenPath/help/examples/bigWigExample.bw"
  remote <- tryCatch(
    read_bigwig(url, chrom = "chr21", start = 33031597, end = 33041570),
    error = function(e) skip(paste("remote fetch unavailable:", conditionMessage(e)))
  )

  expect_s3_class(remote, "tbl_df")
  expect_true(nrow(remote) > 0)
  # a windowed query must return only the window, not the whole 81 MB file
  expect_true(all(remote$chrom == "chr21"))
  expect_true(min(remote$start) >= 33031597)
  expect_true(max(remote$end) <= 33041570)
  # known values from this stable UCSC example file
  expect_equal(remote$start[1], 33031597L)
  expect_equal(remote$value[1], 40)
})

test_that("read_bigbed correctly parses interval coordinates", {
  # Test with the sample bigBed file
  bb_file <- test_path("data/test.bb")
  bb_data <- read_bigbed(bb_file)

  # Basic structure tests
  expect_s3_class(bb_data, "tbl_df")
  expect_true(nrow(bb_data) > 0)
  expect_true(all(c("chrom", "start", "end") %in% names(bb_data)))

  # Test that intervals have proper genomic spans
  # All intervals should have end > start
  expect_true(all(bb_data$end > bb_data$start))

  # Test specific known intervals from your test file
  # Based on your session output, the first feature should span from 4797973 to a much larger end
  first_row <- bb_data[1, ]
  expect_equal(first_row$chrom, "chr1")
  expect_equal(first_row$start, 4797973)
  expect_equal(first_row$end, 4836816) # Should be much larger than start + 1

  # Test that we can parse block structure correctly
  # The testgene should have multiple blocks
  if ("blockSizes" %in% names(bb_data)) {
    block_sizes <- as.numeric(strsplit(first_row$blockSizes, ",")[[1]])
    block_starts <- as.numeric(strsplit(first_row$chromStarts, ",")[[1]])

    expect_gt(length(block_sizes), 1) # Should have multiple blocks
    expect_equal(length(block_sizes), length(block_starts))

    # Calculate total genomic span from blocks
    last_block_end <- first_row$start +
      tail(block_starts, 1) +
      tail(block_sizes, 1)
    expect_equal(first_row$end, last_block_end)
  }

  # Test different chromosomes are properly handled
  unique_chroms <- unique(bb_data$chrom)
  expect_gt(length(unique_chroms), 1)
  expect_true("chr1" %in% unique_chroms)

  # Test with chromosome filtering
  chr1_data <- read_bigbed(bb_file, chrom = "chr1")
  expect_true(all(chr1_data$chrom == "chr1"))
  expect_lt(nrow(chr1_data), nrow(bb_data)) # Should be subset

  # Test with coordinate filtering
  # Get a small region around the first feature
  region_data <- read_bigbed(
    bb_file,
    chrom = "chr1",
    start = 4797000,
    end = 4798000
  )
  expect_true(all(region_data$chrom == "chr1"))
  expect_true(all(region_data$start >= 4797000 | region_data$end <= 4798000))
})

test_that("read_bigbed handles a bigBed with no embedded autoSql schema", {
  # Regression: bbGetSQL() returns NULL when a bigBed has no autoSql schema
  # (sqlOffset == 0). Constructing a std::string from NULL previously crashed
  # R. test_noschema.bb is a bed3 file with its header sqlOffset zeroed.
  bb <- test_path("data/test_noschema.bb")

  res <- read_bigbed(bb)
  expect_s3_class(res, "tbl_df")
  # with no schema there are no extra typed fields: just chrom/start/end
  expect_equal(names(res), c("chrom", "start", "end"))
  expect_gt(nrow(res), 0)
  expect_true(all(res$end > res$start))

  # windowed / per-chrom queries still work without a schema
  chr1 <- read_bigbed(bb, chrom = "chr1")
  expect_true(all(chr1$chrom == "chr1"))
})

test_that("read_bigbed coerces columns based on autoSql types", {
  bb_file <- test_path("data/test.bb")
  bb_data <- read_bigbed(bb_file)

  # Coordinate columns should be integers
  expect_type(bb_data$start, "integer")
  expect_type(bb_data$end, "integer")

  # String columns should be character
  expect_type(bb_data$chrom, "character")
  expect_type(bb_data$name, "character")

  # uint/int columns should be integer (based on BED12 autoSql schema)
  expect_type(bb_data$score, "integer")
  expect_type(bb_data$thickStart, "integer")
  expect_type(bb_data$thickEnd, "integer")
  expect_type(bb_data$reserved, "integer")
  expect_type(bb_data$blockCount, "integer")

  # char[1] should be character
  expect_type(bb_data$strand, "character")

  # Array types (int[blockCount]) should remain character (comma-separated)
  expect_type(bb_data$blockSizes, "character")
  expect_type(bb_data$chromStarts, "character")

  # Verify integer values are correct (not NA or corrupted)
  expect_true(all(bb_data$score >= 0))
  expect_true(all(bb_data$blockCount > 0))
  expect_true(all(bb_data$thickStart >= bb_data$start))
  expect_true(all(bb_data$thickEnd <= bb_data$end))
})
