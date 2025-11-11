test_that("gamma_pdf returns correct values", {
  x <- seq(0.1, 10, by = 0.5)
  alpha <- 2
  beta <- 1

  result <- gamma_pdf(x, alpha, beta)
  expected <- dgamma(x, shape = alpha, rate = beta)

  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("gamma_pdf is zero for negative x", {
  result <- gamma_pdf(x = -5, alpha = 2, beta = 1)
  expect_equal(result, 0)
})

test_that("gamma_cdf matches pgamma", {
  x <- seq(0.1, 10, by = 0.5)
  alpha <- 2
  beta <- 1

  result <- gamma_cdf(x, alpha, beta)
  expected <- pgamma(x, shape = alpha, rate = beta)

  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("gamma_mean is correct", {
  alpha <- 3
  beta <- 2

  expected <- alpha / beta
  result <- gamma_mean(alpha, beta)

  expect_equal(result, expected)
})

test_that("gamma_mode is 0 for alpha <= 1", {
  expect_equal(gamma_mode(alpha = 0.5, beta = 1), 0)
  expect_equal(gamma_mode(alpha = 1, beta = 1), 0)
})

test_that("gamma_init returns reasonable estimates", {
  set.seed(123)
  x <- rgamma(1000, shape = 2, rate = 1)

  init <- gamma_init(x)

  expect_named(init, c("alpha", "beta"))
  expect_equal(init["alpha"], c(alpha=2), tolerance = 0.2)
  expect_equal(init["beta"], c(beta=1), tolerance = 0.2)
})

test_that("drule contains correct structure", {
  expect_named(drule, c("gamma_pdf", "gamma_cdf", "weibull_pdf", "weibull_cdf"))
  expect_length(drule$gamma_pdf, 3)
  expect_named(drule$gamma_pdf, c("x", "alpha", "beta"))
})
