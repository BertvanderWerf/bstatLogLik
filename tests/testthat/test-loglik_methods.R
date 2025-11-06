library(testthat)

test_that("coef.loglik extracts coefficients", {
  # Mock loglik object
  obj <- structure(list(coefficients = c(a = 1, b = 2)), class = "loglik")
  expect_equal(coef(obj), c(a = 1, b = 2))
})

test_that("nobs.loglik returns correct value", {
  obj <- structure(list(nobs = 100), class = "loglik")
  expect_equal(nobs(obj), 100)
})

test_that("df.residual.loglik computes correctly", {
  obj <- structure(list(nobs = 100, npar_estimated = 3), class = "loglik")
  expect_equal(df.residual(obj), 97)
})

test_that("logLik.loglik returns proper logLik object", {
  obj <- structure(list(lnlik = -50, npar_estimated = 3, nobs = 100), class = "loglik")
  ll <- logLik(obj)

  expect_s3_class(ll, "logLik")
  expect_equal(as.numeric(ll), -50)
  expect_equal(attr(ll, "df"), 3)
  expect_equal(attr(ll, "nobs"), 100)
})

test_that("AICc, AIC and BIC compute correctly", {
  obj <- structure(list(lnlik = -50, npar_estimated = 3, nobs = 100), class = "loglik")

  expect_equal(AIC(obj, corrected=TRUE), -2 * (-50) + 2 * 3 + 2*3*4/(100-4))
  expect_equal(AIC(obj, corrected=FALSE), -2 * (-50) + 2 * 3)
  expect_equal(BIC(obj), -2 * (-50) + log(100) * 3)
})

test_that("converged.loglik returns logical", {
  obj1 <- structure(list(converged = TRUE), class = "loglik")
  obj2 <- structure(list(converged = FALSE), class = "loglik")

  expect_true(converged(obj1))
  expect_false(converged(obj2))
})
