test_that("lognormal_pdf returns correct values", {
  x <- c(0.5, 1, 2, 5, 10)
  mu <- 0
  sigma <- 1

  result <- lognormal_pdf(x, mu, sigma)
  expected <- dlnorm(x, meanlog = mu, sdlog = sigma)

  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("lognormal_pdf is zero for non-positive x", {
  expect_equal(lognormal_pdf(x = -5, mu = 0, sigma = 1), 0)
  expect_equal(lognormal_pdf(x = 0, mu = 0, sigma = 1), 0)
})

test_that("lognormal_cdf matches plnorm", {
  x <- c(0.5, 1, 2, 5, 10)
  mu <- 0
  sigma <- 1

  result <- lognormal_cdf(x, mu, sigma)
  expected <- plnorm(x, meanlog = mu, sdlog = sigma)

  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("lognormal_median is correct", {
  mu <- 0
  median_val <- lognormal_median(mu)

  # Median should satisfy CDF(median) = 0.5
  cdf_at_median <- lognormal_cdf(median_val, mu, sigma = 1)
  expect_equal(cdf_at_median, 0.5, tolerance = 1e-10)
})

test_that("lognormal_mean is correct", {
  mu <- 0
  sigma <- 1

  expected <- exp(mu + 0.5 * sigma^2)
  result <- lognormal_mean(mu, sigma)

  expect_equal(result, expected)
})

test_that("lognormal_init returns reasonable estimates", {
  set.seed(123)
  x <- rlnorm(1000, meanlog = 0, sdlog = 1)

  init <- lognormal_init(x)

  expect_named(init, c("mu", "sigma"))
  expect_equal(init["mu"], c(mu=0), tolerance = 0.2)
  expect_equal(init["sigma"], c(sigma=1), tolerance = 0.2)
})

test_that("lognormal_init rejects non-positive values", {
  x <- c(-1, 0, 1, 2, 3)
  expect_error(lognormal_init(x), "No valid positive observations")
})
