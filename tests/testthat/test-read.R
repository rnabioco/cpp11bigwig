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

test_that("remote bigwig_info matches local", {
  skip_on_cran()
  skip_if_not(bigwig_has_curl_cpp(), "built without libcurl")

  url <- paste0(
    "https://raw.githubusercontent.com/rnabioco/cpp11bigwig/",
    "main/inst/extdata/test.bw"
  )
  local <- bigwig_info(system.file("extdata", "test.bw", package = "cpp11bigwig"))
  # a network/remote-host failure should skip, not fail the suite
  remote <- tryCatch(
    bigwig_info(url),
    error = function(e) skip(paste("remote fetch unavailable:", conditionMessage(e)))
  )

  expect_equal(remote, local)
})

test_that("remote bigbed_info matches local", {
  skip_on_cran()
  skip_if_not(bigwig_has_curl_cpp(), "built without libcurl")

  url <- paste0(
    "https://raw.githubusercontent.com/rnabioco/cpp11bigwig/",
    "main/inst/extdata/test.bb"
  )
  local <- bigbed_info(system.file("extdata", "test.bb", package = "cpp11bigwig"))
  # a network/remote-host failure should skip, not fail the suite
  remote <- tryCatch(
    bigbed_info(url),
    error = function(e) skip(paste("remote fetch unavailable:", conditionMessage(e)))
  )

  expect_equal(remote, local)
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

  res <- suppressMessages(read_bigbed(bb))
  expect_s3_class(res, "tbl_df")
  # with no schema there are no extra typed fields: just chrom/start/end
  expect_equal(names(res), c("chrom", "start", "end"))
  expect_gt(nrow(res), 0)
  expect_true(all(res$end > res$start))

  # windowed / per-chrom queries still work without a schema
  chr1 <- suppressMessages(read_bigbed(bb, chrom = "chr1"))
  expect_true(all(chr1$chrom == "chr1"))
})

test_that("read_bigbed recovers BED columns for a schema-less bed12", {
  # A bed12 file written by `bedToBigBed` without `-as` has no autoSql schema
  # but still stores 12 fields per record. The reader should fall back to the
  # header's fieldCount/definedFieldCount and name the columns with the standard
  # BED field names instead of dropping everything past chrom/start/end.
  # test_noschema_bed12.bb is test.bb with its header sqlOffset zeroed.
  noschema <- test_path("data/test_noschema_bed12.bb")
  schema <- test_path("data/test.bb")

  res <- suppressMessages(read_bigbed(noschema))
  expect_s3_class(res, "tbl_df")
  expect_equal(ncol(res), 12)
  expect_equal(
    names(res),
    c(
      "chrom", "start", "end", "name", "score", "strand",
      "thickStart", "thickEnd", "itemRgb", "blockCount",
      "blockSizes", "blockStarts"
    )
  )

  # schema-less output should match the schema-backed file column-for-column,
  # apart from the (schema-derived) column names
  ref <- read_bigbed(schema)
  expect_equal(nrow(res), nrow(ref))
  expect_equal(unname(as.list(res)), unname(as.list(ref)))
})

test_that("read_bigbed messages when a bigBed has no autoSql, stays silent otherwise", {
  # column names for a schema-less file are inferred, so the user is told
  expect_message(
    read_bigbed(test_path("data/test_noschema_bed12.bb")),
    "autoSql"
  )
  expect_message(
    read_bigbed(test_path("data/test_noschema.bb")),
    "autoSql"
  )
  # a schema-backed file declares its own names: no message
  expect_no_message(read_bigbed(test_path("data/test.bb")))
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

test_that("bigbed_info reports header metadata and autoSql", {
  info <- bigbed_info(test_path("data/test.bb"))

  expect_named(
    info,
    c(
      "version", "n_chroms", "field_count", "defined_field_count",
      "n_bases_covered", "autosql"
    )
  )
  # test.bb is a genuine BED12 with an embedded schema
  expect_equal(info$field_count, 12L)
  expect_equal(info$defined_field_count, 12L)
  expect_gt(info$n_chroms, 0L)
  expect_match(info$autosql, "blockSizes")
})

test_that("bigbed_info reads field counts from a schema-less bed12 header", {
  # no embedded autoSql (sqlOffset zeroed) but the header still records 12 fields
  info <- bigbed_info(test_path("data/test_noschema_bed12.bb"))

  expect_equal(info$field_count, 12L)
  expect_equal(info$defined_field_count, 12L)
  expect_identical(info$autosql, "")
})

test_that("bigbed_info distinguishes a bed3 from a bed12", {
  bed3 <- bigbed_info(test_path("data/test_noschema.bb"))
  expect_equal(bed3$defined_field_count, 3L)
})

test_that("bigwig_info reports header metadata and summary stats", {
  info <- bigwig_info(test_path("data/test.bw"))

  expect_named(
    info,
    c(
      "version", "n_levels", "n_chroms", "n_bases_covered",
      "min", "max", "mean", "std"
    )
  )
  expect_gt(info$n_chroms, 0L)
  expect_gt(info$n_bases_covered, 0)
  expect_true(info$max >= info$min)
})

test_that("info functions reject the wrong file type and missing files", {
  bb <- test_path("data/test.bb")
  bw <- test_path("data/test.bw")

  expect_error(bigbed_info(bw), "Not a bigBed")
  expect_error(bigwig_info(bb), "Not a bigWig")
  expect_error(bigbed_info("does-not-exist.bb"), "does not exist")
})
