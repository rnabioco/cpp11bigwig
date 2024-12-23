
test_that("results have expected shape", {
  bw <- test_path("data/test.bw")

  res <- read_bigwig(bw)
  expect_equal(ncol(res), 4)
  expect_equal(nrow(res), 6)

  # interval sizes
  expect_true(all(res$end - res$start == c(1,1,1,50,1,100)))
  # interval values
  expect_equal(sum(res$value), 5.5)

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
