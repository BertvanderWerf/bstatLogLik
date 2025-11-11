test_that("evaluate_function computes cdfWeibull correctly", {
  cdfWeibull <- function(x, k, lambda, c) { 1 - exp(-((x + c)/lambda)^k) }
  testdata <- data.frame(Age = 0:10, Sex=factor(rep(c('Male', 'Female'), c(5,6))))
  parameter_info <- create_parameter_info(
    c(k = 2, lambda = 5, c = 0),
    link = list(lambda='log'),
    formulas = list(k=~0+Sex),
    model_data = testdata
  )
  #library(bstatUtils)
  testdata$cdf <- evaluate_function(fun = cdfWeibull, vars = c(x = "Age"), parameter_info = parameter_info, data = testdata)
  expect_true(is.numeric(testdata$cdf))
  expect_equal(length(testdata$cdf), nrow(testdata))
  expect_true(all(testdata$cdf >= 0 & testdata$cdf <= 1))
})

test_that("error raised for missing parameters", {
  parameter_info <- create_parameter_info(c(a = 1), model_data = data.frame(x = 1:3))
  expect_error(
    evaluate_function(function(x) x + b, vars = c(x = "x"), parameter_info = parameter_info, data = data.frame(x = 1:3))
  )
})

test_that("error raised for missing variable", {
  parameter_info <- create_parameter_info(c(b = 2), model_data = data.frame(x = 1:3))
  expect_error(
    evaluate_function(function(x) x + b, vars = c(y = "x"), parameter_info = parameter_info, data = data.frame(x = 1:3))
  )
})
