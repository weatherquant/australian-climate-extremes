testthat::test_that("canonical constants are internally consistent", {
  source("R/utilities/constants.R")
  testthat::expect_true(all(c("raw", "homogenised") %in% DATA_VERSIONS))
  testthat::expect_true(all(c("tmax", "tmin", "rainfall") %in% CANONICAL_VARIABLES))
})
