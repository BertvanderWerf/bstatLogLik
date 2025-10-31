test_that("weibull_pdf returns correct values", {
  # Compare with standard dweibull (2-parameter case, c=0)
  x <- seq(0, 10, by = 0.5)
  k <- 2
  lambda <- 5

  result <- weibull_pdf(x, k, lambda, c = 0)
  expected <- dweibull(x, shape = k, scale = lambda)

  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("weibull_cdf returns correct values", {
  # Compare with standard pweibull
  x <- seq(0, 10, by = 0.5)
  k <- 2
  lambda <- 5

  result <- weibull_cdf(x, k, lambda, c = 0)
  expected <- pweibull(x, shape = k, scale = lambda)

  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("weibull_pdf integrates to 1", {
  k <- 2
  lambda <- 5

  # Numerical integration over support
  integral <- integrate(weibull_pdf, lower = 0, upper = Inf, k = k, lambda = lambda, c = 0)

  expect_equal(integral$value, 1, tolerance = 1e-4)
})

test_that("weibull_mean matches theoretical value", {
  k <- 2
  lambda <- 5
  c <- 0

  # Theoretical mean
  expected <- lambda * gamma(1 + 1/k) - c
  result <- weibull_mean(k, lambda, c)

  expect_equal(result, expected)
})

test_that("weibull_median is correct", {
  k <- 2
  lambda <- 5
  c <- 0

  # Median should satisfy CDF(median) = 0.5
  median_val <- weibull_median(k, lambda, c)
  cdf_at_median <- weibull_cdf(median_val, k, lambda, c)

  expect_equal(cdf_at_median, 0.5, tolerance = 1e-10)
})

test_that("weibull_mode is 0 for k <= 1", {
  expect_equal(weibull_mode(k = 0.5, lambda = 5, c = 0), 0)
  expect_equal(weibull_mode(k = 1, lambda = 5, c = 0), 0)
})

test_that("weibull_mode is positive for k > 1", {
  mode_val <- weibull_mode(k = 2, lambda = 5, c = 0)
  expect_true(mode_val > 0)
})

test_that("weibull_hazard is non-negative", {
  x <- seq(0, 10, by = 0.5)
  h <- weibull_hazard(x, k = 2, lambda = 5, c = 0)

  expect_true(all(h >= 0))
})

test_that("weibull_hazard is 0 for negative x", {
  h <- weibull_hazard(x = -5, k = 2, lambda = 5, c = 0)
  expect_equal(h, 0)
})

test_that("weibull_init returns reasonable estimates", {
  set.seed(123)
  x <- rweibull(1000, shape = 2, scale = 5)

  init <- weibull_init(x)

  expect_named(init, c("k", "lambda", "c"))
  expect_true(init["k"] > 0)
  expect_true(init["lambda"] > 0)
  expect_equal(init["c"], c(c=0))

  # Should be reasonably close to true values
  expect_equal(init["k"], c(k=2), tolerance = 0.3)
  expect_equal(init["lambda"], c(lambda=5), tolerance = 1)
})

test_that("weibull_init handles difftime objects", {
  x <- as.difftime(c(1, 2, 3, 4, 5), units = "hours")
  init <- weibull_init(x)

  expect_named(init, c("k", "lambda", "c"))
  expect_true(all(is.finite(init)))
})

test_that("weibull_init handles weights", {
  x <- c(1, 2, 3, 4, 5)
  w <- c(10, 1, 1, 1, 10)

  init_unweighted <- weibull_init(x)
  init_weighted <- weibull_init(x, weights = w)

  # Weighted version should differ due to different effective mean/variance
  expect_false(identical(init_unweighted, init_weighted))
})

test_that("drule contains correct structure", {
  expect_type(drule, "list")
  #expect_named(drule, c("weibull_pdf", "weibull_cdf"))

  expect_length(drule$weibull_pdf, 4)
  expect_named(drule$weibull_pdf, c("x", "k", "lambda", "c"))

  expect_length(drule$weibull_cdf, 4)
  expect_named(drule$weibull_cdf, c("x", "k", "lambda", "c"))
})

test_that("location parameter c shifts distribution correctly", {
  x <- seq(0, 10, by = 1)
  k <- 2
  lambda <- 5
  c_shift <- 2

  # PDF with c=0 at x should equal PDF with c=c_shift at x-c_shift
  pdf_no_shift <- weibull_pdf(x, k, lambda, c = 0)
  pdf_shifted <- weibull_pdf(x + c_shift, k, lambda, c = 0)
  pdf_with_c <- weibull_pdf(x, k, lambda, c = c_shift)

  expect_equal(pdf_shifted, pdf_with_c, tolerance = 1e-10)
})

