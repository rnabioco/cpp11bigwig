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

test_that("missing file causes error", {
  expect_snapshot_error(read_bigwig("missing.bw"))
})

test_that("negative coords causes error", {
  expect_snapshot_error(read_bigwig(test_path("data/test.bw"), start = -1))
  expect_snapshot_error(read_bigwig(test_path("data/test.bw"), end = -1))
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
