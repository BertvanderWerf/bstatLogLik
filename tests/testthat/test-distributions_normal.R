test_that("normal_pdf returns correct values", {
  x <- -3:3
  mu <- 0
  sigma <- 1

  result <- normal_pdf(x, mu, sigma)
  expected <- dnorm(x, mean = mu, sd = sigma)

  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("normal_cdf matches pnorm", {
  x <- -3:3
  mu <- 0
  sigma <- 1

  result <- normal_cdf(x, mu, sigma)
  expected <- pnorm(x, mean = mu, sd = sigma)

  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("normal_mean returns mu", {
  result <- normal_mean(mu = 5)
  expect_equal(result, 5)
})

test_that("normal_median equals mean", {
  result <- normal_median(mu = 5)
  expect_equal(result, 5)
})

test_that("normal_mode equals mean", {
  result <- normal_mode(mu = 5)
  expect_equal(result, 5)
})

test_that("normal_init returns reasonable estimates", {
  set.seed(123)
  x <- rnorm(1000, mean = 5, sd = 2)

  init <- normal_init(x)

  expect_named(init, c("mu", "sigma"))
  expect_equal(init["mu"], 5, tolerance = 0.2)
  expect_equal(init["sigma"], 2, tolerance = 0.2)
})
