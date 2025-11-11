# Unit Tests for set_parameter_values Function
# File: tests/testthat/test_set_parameter_values.R

describe("set_parameter_values", {

  # ==========================================================================
  # Test 1: Basic Functionality
  # ==========================================================================

  it("updates simple parameter values", {
    expect_warning(
    param_info <- create_parameter_info(
      parameter_values = c(k = 2, lambda = 10),
      estimate_flag = c(k = 1, lambda = 1)
    ),
    "model_data is NULL; no design matrices are created. Pass model_data to expand parameters via formulas.")

    param_info_updated <- set_parameter_values(
      param_info,
      newvals = c(k = 3, lambda = 12)
    )

    expect_equal(param_info_updated$parameter_values$k, 3)
    expect_equal(param_info_updated$parameter_values$lambda, 12)
  })

  it("updates single parameter", {
    param_info <- create_parameter_info(
      parameter_values = c(k = 2, lambda = 10)
    )

    param_info_updated <- set_parameter_values(
      param_info,
      newvals = c(k = 5)
    )

    expect_equal(param_info_updated$parameter_values$k, 5)
    expect_equal(param_info_updated$parameter_values$lambda, 10)  # Unchanged
  })

  it("returns unchanged param_info when newvals is NULL", {
    param_info <- create_parameter_info(
      parameter_values = c(k = 2, lambda = 10)
    )

    param_info_updated <- set_parameter_values(param_info, newvals = NULL)

    expect_equal(param_info_updated$parameter_values$k, 2)
    expect_equal(param_info_updated$parameter_values$lambda, 10)
  })

  # ==========================================================================
  # Test 2: With Formula Expansion
  # ==========================================================================

  it("updates expanded parameters from formulas", {
    data <- data.frame(group = factor(c("A", "B", "A", "B")))

    param_info <- create_parameter_info(
      parameter_values = c(k = 2, lambda = 10),
      formulas = list(k = ~0 + group, lambda = ~1),
      model_data = data
    )

    # Should have kA, kB, lambda
    param_info_updated <- set_parameter_values(
      param_info,
      newvals = c(kA = 2.5, kB = 3.5)
    )

    expect_equal(param_info_updated$parameter_values$k[["A"]], 2.5)
    expect_equal(param_info_updated$parameter_values$k[["B"]], 3.5)
    expect_equal(param_info_updated$parameter_values$lambda, 10)  # Unchanged
  })

  it("updates all expanded parameters", {
    data <- data.frame(group = factor(c("A", "B", "A", "B")))

    param_info <- create_parameter_info(
      parameter_values = c(k = 2, lambda = 10),
      formulas = list(k = ~0 + group, lambda = ~0 + group),
      model_data = data
    )

    param_info_updated <- set_parameter_values(
      param_info,
      newvals = c(kA = 1, kB = 2, lambdaA = 8, lambdaB = 12)
    )

    expect_equal(param_info_updated$parameter_values$k[["A"]], 1)
    expect_equal(param_info_updated$parameter_values$k[["B"]], 2)
    expect_equal(param_info_updated$parameter_values$lambda[["A"]], 8)
    expect_equal(param_info_updated$parameter_values$lambda[["B"]], 12)
  })

  # ==========================================================================
  # Test 3: Input Validation
  # ==========================================================================

  it("rejects non-parameter_info input", {
    expect_error(
      set_parameter_values(
        list(parameter_values = c(k = 2)),
        newvals = c(k = 3)
      ),
      "must be of class 'parameter_info'"
    )
  })

  it("rejects non-numeric newvals", {
    param_info <- create_parameter_info(
      parameter_values = c(k = 2)
    )

    expect_error(
      set_parameter_values(param_info, newvals = c(k = "three")),
      "must be numeric"
    )
  })

  it("rejects unnamed newvals", {
    param_info <- create_parameter_info(
      parameter_values = c(k = 2, lambda = 10)
    )

    expect_error(
      set_parameter_values(param_info, newvals = c(3, 12)),
      "must be a named vector"
    )
  })

  it("rejects newvals with empty names", {
    param_info <- create_parameter_info(
      parameter_values = c(k = 2, lambda = 10)
    )

    newvals <- c(k = 3, 12)
    names(newvals) <- c("k", "")

    expect_error(
      set_parameter_values(param_info, newvals = newvals),
      "must be a named vector"
    )
  })

  it("rejects non-existent parameter names", {
    param_info <- create_parameter_info(
      parameter_values = c(k = 2, lambda = 10)
    )

    expect_error(
      set_parameter_values(param_info, newvals = c(theta = 5)),
      "does not exist"
    )
  })

  it("rejects non-existent sub-parameter names", {
    data <- data.frame(group = factor(c("A", "B")))

    param_info <- create_parameter_info(
      parameter_values = c(k = 2),
      formulas = list(k = ~0 + group),
      model_data = data
    )

    # kC doesn't exist (only kA and kB)
    expect_error(
      set_parameter_values(param_info, newvals = c(kC = 3)),
      "does not exist"
    )
  })

  # ==========================================================================
  # Test 4: Partial Updates
  # ==========================================================================

  it("allows partial updates with formulas", {
    data <- data.frame(group = factor(c("A", "B", "C")))

    param_info <- create_parameter_info(
      parameter_values = c(k = 2),
      formulas = list(k = ~0 + group),
      model_data = data
    )

    # Update only kA and kC, leave kB unchanged
    param_info_updated <- set_parameter_values(
      param_info,
      newvals = c(kA = 1.5, kC = 3.5)
    )

    expect_equal(param_info_updated$parameter_values$k[["A"]], 1.5)
    expect_equal(param_info_updated$parameter_values$k[["B"]], 2)  # Unchanged
    expect_equal(param_info_updated$parameter_values$k[["C"]], 3.5)
  })

  # ==========================================================================
  # Test 5: Multiple Parameters
  # ==========================================================================

  it("updates multiple parameters simultaneously", {
    param_info <- create_parameter_info(
      parameter_values = c(k = 2, lambda = 10, c = 0, sigma = 1)
    )

    param_info_updated <- set_parameter_values(
      param_info,
      newvals = c(k = 3, lambda = 15, c = 0.5, sigma = 1.2)
    )

    expect_equal(param_info_updated$parameter_values$k, 3)
    expect_equal(param_info_updated$parameter_values$lambda, 15)
    expect_equal(param_info_updated$parameter_values$c, 0.5)
    expect_equal(param_info_updated$parameter_values$sigma, 1.2)
  })

  # ==========================================================================
  # Test 6: Custom Separator
  # ==========================================================================

  it("handles custom separator", {
    data <- data.frame(group = factor(c("A", "B")))

    param_info <- create_parameter_info(
      parameter_values = c(k = 2),
      formulas = list(k = ~0 + group),
      model_data = data,
      sep = "_"  # Custom separator
    )

    # With "_" separator: k_A, k_B
    param_info_updated <- set_parameter_values(
      param_info,
      newvals = c(k_A = 2.5, k_B = 3.5),
      sep = "_"
    )

    expect_equal(param_info_updated$parameter_values$k[["A"]], 2.5)
    expect_equal(param_info_updated$parameter_values$k[["B"]], 3.5)
  })

  # ==========================================================================
  # Test 7: Numeric Values
  # ==========================================================================

  it("accepts negative values", {
    param_info <- create_parameter_info(
      parameter_values = c(k = 2, lambda = 10)
    )

    param_info_updated <- set_parameter_values(
      param_info,
      newvals = c(k = -5, lambda = -2)
    )

    expect_equal(param_info_updated$parameter_values$k, -5)
    expect_equal(param_info_updated$parameter_values$lambda, -2)
  })

  it("accepts zero values", {
    param_info <- create_parameter_info(
      parameter_values = c(k = 2, lambda = 10)
    )

    param_info_updated <- set_parameter_values(
      param_info,
      newvals = c(k = 0, lambda = 0)
    )

    expect_equal(param_info_updated$parameter_values$k, 0)
    expect_equal(param_info_updated$parameter_values$lambda, 0)
  })

  it("accepts very small values", {
    param_info <- create_parameter_info(
      parameter_values = c(k = 2)
    )

    param_info_updated <- set_parameter_values(
      param_info,
      newvals = c(k = 1e-10)
    )

    expect_equal(param_info_updated$parameter_values$k, 1e-10)
  })

  it("accepts very large values", {
    param_info <- create_parameter_info(
      parameter_values = c(k = 2)
    )

    param_info_updated <- set_parameter_values(
      param_info,
      newvals = c(k = 1e10)
    )

    expect_equal(param_info_updated$parameter_values$k, 1e10)
  })

  # ==========================================================================
  # Test 8: Integration with Design of Experiments
  # ==========================================================================

  it("works with parameter design output", {
    param_info <- create_parameter_info(
      parameter_values = c(k = 2, lambda = 10),
      lower_bounds = c(k = 1, lambda = 5),
      upper_bounds = c(k = 4, lambda = 15),
      estimate_flag = c(k = 1, lambda = 1)
    )

    # Create design
    design <- create_parameter_design(param_info)

    # Simulate selecting best point
    best_point <- design[5, ]  # Arbitrary point

    # Extract parameter values
    new_values <- c(
      k = best_point$k,
      lambda = best_point$lambda
    )

    # Update param_info
    param_info_updated <- set_parameter_values(param_info, new_values)

    expect_equal(param_info_updated$parameter_values$k, new_values[["k"]])
    expect_equal(param_info_updated$parameter_values$lambda, new_values[["lambda"]])
  })

  # ==========================================================================
  # Test 9: Return Value
  # ==========================================================================

  it("returns parameter_info object", {
    param_info <- create_parameter_info(
      parameter_values = c(k = 2)
    )

    result <- set_parameter_values(param_info, newvals = c(k = 3))

    expect_s3_class(result, "parameter_info")
  })

  it("does not modify original param_info", {
    param_info <- create_parameter_info(
      parameter_values = c(k = 2, lambda = 10)
    )

    original_k <- param_info$parameter_values$k

    # Update returns new object
    param_info_updated <- set_parameter_values(
      param_info,
      newvals = c(k = 5)
    )

    # Original unchanged (if R copy-on-modify works)
    expect_equal(param_info$parameter_values$k, original_k)
    expect_equal(param_info_updated$parameter_values$k, 5)
  })
})
