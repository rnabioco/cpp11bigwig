
test_that("output have expected shape", {
  out <- read_bigwig(test_path("data/test.bw"))
  expect_equal(ncol(out), 4)
  expect_equal(nrow(out), 154)
})

test_that("missing file causes error", {
  expect_error(read_bigwig("missing.bw"))
})

test_that("negative coords causes error", {
  expect_error(read_bigwig(test_path("data/test.bw"), start = -1))
  expect_error(read_bigwig(test_path("data/test.bw"), end = -1))
})
