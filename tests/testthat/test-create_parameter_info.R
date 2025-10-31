test_that("create_parameter_info handles minimal input", {
  info <- suppressWarnings({
    create_parameter_info(parameter_values = c(alpha = 2))
  })
  expect_s3_class(info, "parameter_info")
  expect_equal(info$parameter_values$alpha, 2)
  expect_equal(info$lower_bounds$alpha, -Inf)
  expect_equal(info$upper_bounds$alpha, Inf)
  expect_equal(info$estimate_flag$alpha, 1)
  expect_equal(info$link_functions$alpha, "identity")
  expect_true(class(info$formulas$alpha)=='formula')
  expect_null(info$design_matrices$alpha)
})

test_that("create_parameter_info supports updating slots from existing info", {
  df <- data.frame(A = 1:3)
  info1 <- create_parameter_info(parameter_values = c(b = 1), formulas = list(b = ~A), model_data = df)
  info2 <- create_parameter_info(parameter_values = info1, formulas = list(b = ~1))
  expect_equal(info2$formulas$b, ~1)
  expect_equal(info2$parameter_values$b, c(Intercept=1))
  expect_identical(info2$model_data, df)
})

test_that("create_parameter_info detects missing or duplicate names", {
  expect_error(create_parameter_info(c(1, 2)), "must be a named vector")
  expect_error(create_parameter_info(c(a=1, a=2)), "must be a named vector")
})

test_that("create_parameter_info handles partial and named updates", {
  df <- data.frame(A = 1:2)
  info <- create_parameter_info(
    parameter_values = c(u=1, v=2),
    formulas = list(u = ~A),
    model_data = df
  )
  expect_true(class(info$formulas$u)=='formula')
  expect_true(class(info$formulas$v)=='formula')
})

test_that("output contains all slots and model_data is retained", {
  df <- data.frame(B = 1)
  info <- create_parameter_info(parameter_values = c(x = 3), model_data = df)
  expect_named(info, c("parameter_values", "lower_bounds", "upper_bounds", "estimate_flag",
                       "link_functions", "formulas", "design_matrices", "model_data", "sep"))
  expect_identical(info$model_data, df)
})

test_that("create_parameter_info allows changing separator", {
  info <- suppressWarnings({
    create_parameter_info(parameter_values = c(y = 1), sep = ":")
  })
  expect_equal(info$sep, ":")
})

test_that("create_parameter_info sets estimate_flag", {
  df <- expand.grid(Capital = c("A", "B"), Lower = c("a", "b"))
  param_info <- create_parameter_info(
    parameter_values = c(k = 2, lambda = 10, c = 0),
    estimate_flag = c(k = 1, lambda = 1, c = 0),
    formulas = list(k = ~0+Capital, lambda = ~1, c = ~Capital*Lower),
    link = c(lambda = 'log'),
    model_data = df
  )
  expect_equal(param_info$estimate_flag$k[1], 1)
  expect_equal(param_info$estimate_flag$lambda[1], 1)
  expect_equal(param_info$estimate_flag$c[1], 0)
})

test_that("create_parameter_info returns expected structure", {
  df <- expand.grid(Capital = c("A", "B"), Lower = c("a", "b"))
  result <- create_parameter_info(
    parameter_values = c(k = 2, lambda = 10, c = 0),
    formulas = list(k = ~0+Capital, lambda = ~1, c = ~Capital*Lower),
    model_data = df
  )
  expect_s3_class(result, "parameter_info")
  expect_named(result$parameter_values, c("k", "lambda", "c"))
  expect_true(!is.null(result$design_matrices$k))
  expect_true(ncol(result$design_matrices$c) >= 1)
})

test_that("summary.parameter_info returns info and parameters", {
  df <- expand.grid(Capital = c("A", "B"), Lower = c("a", "b"))
  param_info <- create_parameter_info(
    parameter_values = c(k = 2, lambda = 10, c = 0),
    formulas = list(k = ~0+Capital, lambda = ~1, c = ~Capital*Lower),
    model_data = df
  )
  result <- summary.parameter_info(param_info)
  expect_true("info" %in% names(result))
  expect_true("parameters" %in% names(result))
  expect_equal(nrow(result$info), 3)
})

