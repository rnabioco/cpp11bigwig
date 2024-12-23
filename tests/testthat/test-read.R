
test_that("results have expected shape", {
  bw <- test_path("data/test.bw")

  res <- read_bigwig(bw)
  expect_equal(ncol(res), 4)
  expect_equal(nrow(res), 6)
  expect_true("tbl_df" %in% class(res))

  # interval sizes
  expect_true(all(res$end - res$start == c(1,1,1,50,1,100)))
  # interval values
  expect_equal(sum(res$value), 5.5)

  # params work
  expect_equal(nrow(read_bigwig(bw, chrom = "1")), 5)
  expect_equal(nrow(read_bigwig(bw, start = 100)), 3)
  expect_equal(nrow(read_bigwig(bw, end = 3)), 3)

  # GRanges
  res <- read_bigwig(bw, as = 'gr')
  expect_true("GRanges" %in% class(res))

})

test_that("missing file causes error", {
  expect_snapshot_error(read_bigwig("missing.bw"))
})

test_that("negative coords causes error", {
  expect_snapshot_error(read_bigwig(test_path("data/test.bw"), start = -1))
  expect_snapshot_error(read_bigwig(test_path("data/test.bw"), end = -1))
})
