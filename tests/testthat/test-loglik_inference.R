test_that("lrt correctly computes LR statistic", {
  # Mock nested models
  model1 <- structure(list(
    lnlik = -100,
    npar_estimated = 2,
    nobs = 100,
    converged = TRUE
  ), class = "loglik")

  model2 <- structure(list(
    lnlik = -95,
    npar_estimated = 4,
    nobs = 100,
    converged = TRUE
  ), class = "loglik")

  result <- lrt(model1, model2)

  # LR = -2 * (-100 - (-95)) = 10
  expect_equal(result$LR_stat[2], 10)
  expect_equal(result$Df_diff[2], 2)

  # P-value should be pchisq(10, df=2, lower.tail=FALSE)
  expect_equal(result$Pr_Chi[2], pchisq(10, df = 2, lower.tail = FALSE))
})

test_that("lrt orders models correctly", {
  # Models given in wrong order
  model_complex <- structure(list(
    lnlik = -95, npar_estimated = 4, nobs = 100, converged = TRUE
  ), class = "loglik")

  model_simple <- structure(list(
    lnlik = -100, npar_estimated = 2, nobs = 100, converged = TRUE
  ), class = "loglik")

  # Should auto-order by df
  result <- lrt(model_complex, model_simple)

  expect_equal(result$Df[1], 2)
  expect_equal(result$Df[2], 4)
})

test_that("lrt detects non-nested models", {
  model1 <- structure(list(
    lnlik = -95, npar_estimated = 3, nobs = 100, converged = TRUE
  ), class = "loglik")

  model2 <- structure(list(
    lnlik = -100, npar_estimated = 3, nobs = 100, converged = TRUE
  ), class = "loglik")

  expect_warning(lrt(model1, model2), "same df")
})

test_that("lrt detects different sample sizes", {
  model1 <- structure(list(
    lnlik = -100, npar_estimated = 2, nobs = 100, converged = TRUE
  ), class = "loglik")

  model2 <- structure(list(
    lnlik = -95, npar_estimated = 4, nobs = 150, converged = TRUE
  ), class = "loglik")

  expect_error(lrt(model1, model2), "same number of observations")
})

test_that("anova.loglik calls lrt", {
  model1 <- structure(list(
    lnlik = -100, npar_estimated = 2, nobs = 100, converged = TRUE
  ), class = "loglik")

  model2 <- structure(list(
    lnlik = -95, npar_estimated = 4, nobs = 100, converged = TRUE
  ), class = "loglik")

  result <- anova(model1, model2)

  expect_s3_class(result, "lrt")
  expect_equal(nrow(result), 2)
})

test_that("anova.loglik rejects non-LRT tests", {
  model1 <- structure(list(
    lnlik = -100, npar_estimated = 2, nobs = 100, converged = TRUE
  ), class = "loglik")

  model2 <- structure(list(
    lnlik = -95, npar_estimated = 4, nobs = 100, converged = TRUE
  ), class = "loglik")

  expect_error(anova(model1, model2, test = "F"), "Only likelihood ratio test")
})
