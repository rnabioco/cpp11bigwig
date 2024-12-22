
test_that("results have expected shape", {
  # -- raw output -------
  bw <- test_path("data/test.bw")

  res <- read_bigwig(bw)
  expect_equal(ncol(res), 4)
  expect_equal(nrow(res), 154)

  # chrom
  expect_equal(nrow(read_bigwig(bw, chrom = "1")), 54)
  expect_equal(nrow(read_bigwig(bw, chrom = "10")), 100)

  # `chrom` must be a string
  expect_snapshot_error(read_bigwig(bw, chrom = 10))

  # 0 rows, chrom doesn't exist
  expect_equal(nrow(read_bigwig(bw, chrom = "100")), 0)

  # start/end
  expect_equal(nrow(read_bigwig(bw, chrom = "1", start = 1, end = 100)), 2)
  expect_equal(nrow(read_bigwig(bw, chrom = "1", start = 100, end = 100)), 0)
  expect_equal(nrow(read_bigwig(bw, chrom = "10", start = 100, end = 1000)), 100)
  expect_equal(nrow(read_bigwig(bw, chrom = "10", start = 200, end = 250)), 50)

  # values
  res <- read_bigwig(bw, chrom = "10", start = 100, end = 1000)
  expect_equal(length(unique(res$value)), 1)
  expect_equal(unique(res$value), 2)

  # -- bedgraph output -----
})

test_that("missing file causes error", {
  expect_snapshot_error(read_bigwig("missing.bw"))
})

test_that("negative coords causes error", {
  expect_snapshot_error(read_bigwig(test_path("data/test.bw"), start = -1))
  expect_snapshot_error(read_bigwig(test_path("data/test.bw"), end = -1))
})
