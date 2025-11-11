describe("create_parameter_design", {

  # ==========================================================================
  # Test 1: Basic Functionality
  # ==========================================================================

  it("creates design for simple 2-parameter case", {
    expect_warning(
      param_info <- create_parameter_info(
        parameter_values = c(k = 2, lambda = 10),
        lower_bounds = c(k = 1, lambda = 5),
        upper_bounds = c(k = 4, lambda = 15),
        estimate_flag = c(k = 1, lambda = 1)
      ),
      "model_data is NULL; no design matrices are created. Pass model_data to expand parameters via formulas.")

    design <- create_parameter_design(param_info)

    expect_s3_class(design, "parameter_design")
    expect_s3_class(design, "data.frame")
    expect_true("k" %in% names(design))
    expect_true("lambda" %in% names(design))
    expect_true("point_type" %in% names(design))
  })

  it("creates correct number of design points", {
    expect_warning(
    param_info <- create_parameter_info(
      parameter_values = c(k = 2, lambda = 10),
      lower_bounds = c(k = 1, lambda = 5),
      upper_bounds = c(k = 4, lambda = 15),
      estimate_flag = c(k = 1, lambda = 1)
    ),
    "model_data is NULL; no design matrices are created. Pass model_data to expand parameters via formulas.")

    design <- create_parameter_design(param_info, center_points = 3)

    # 2 params: 4 factorial + 4 axial + 3 center = 11
    expect_equal(nrow(design), 11)
  })

  # ==========================================================================
  # Test 2: Parameter Scaling
  # ==========================================================================

  it("transforms parameters to correct scale", {
    param_info <- create_parameter_info(
      parameter_values = c(k = 2.5, lambda = 10),
      lower_bounds = c(k = 1, lambda = 5),
      upper_bounds = c(k = 4, lambda = 15),
      estimate_flag = c(k = 1, lambda = 1),
      model_data=data.frame(x=1)
    )

    design <- create_parameter_design(param_info)

    # Center point should be at midpoint of bounds
    center_points <- design[design$point_type == "center", ]

    # Center = (lower + upper) / 2
    expect_equal(center_points$k[1], 2.5, tolerance = 1e-6)
    expect_equal(center_points$lambda[1], 10, tolerance = 1e-6)
  })

  it("factorial points are at bounds for 2-level factorial", {
    param_info <- create_parameter_info(
      parameter_values = c(k = 2.5),
      lower_bounds = c(k = 1),
      upper_bounds = c(k = 4),
      estimate_flag = c(k = 1),
      model_data=data.frame(x=1)
    )

    design <- create_parameter_design(param_info, alpha = "faces")

    factorial_points <- design[design$point_type == "factorial", "k"]

    # Should be at lower and upper bounds for face-centered
    expect_true(all(factorial_points %in% c(1, 4)))
  })

  # ==========================================================================
  # Test 3: With Generators (Fractional Design)
  # ==========================================================================

  it("accepts generators for fractional designs", {
    param_info <- create_parameter_info(
      parameter_values = c(k = 2, lambda = 10, c = 0, d = 1),
      lower_bounds = c(k = 1, lambda = 5, c = -1, d = 0),
      upper_bounds = c(k = 4, lambda = 15, c = 1, d = 2),
      estimate_flag = c(k = 1, lambda = 1, c = 1, d = 1),
      model_data=data.frame(x=1)
    )

    design <- create_parameter_design(
      param_info,
      generators = "D=ABC"  # Using coded factor names
    )

    # 2^(4-1) = 8 factorial points
    expect_equal(sum(design$point_type == "factorial"), 8)
  })

  # ==========================================================================
  # Test 4: Fixed Parameters
  # ==========================================================================

  it("handles fixed parameters correctly", {
    param_info <- create_parameter_info(
      parameter_values = c(k = 2, lambda = 10, c = 0),
      lower_bounds = c(k = 1, lambda = 5, c = -1),
      upper_bounds = c(k = 4, lambda = 15, c = 1),
      estimate_flag = c(k = 1, lambda = 1, c = 0),  # c is fixed
      model_data=data.frame(x=1)
    )

    design <- create_parameter_design(param_info)

    # Should only design over k and lambda
    expect_equal(attr(design, "num_parameters"), 2)

    # c should be constant at 0
    expect_true(all(design$c == 0))
  })

  it("includes fixed parameters in output", {
    expect_warning(
      param_info <- create_parameter_info(
      parameter_values = c(k = 2, lambda = 10, c = 5),
      estimate_flag = c(k = 1, lambda = 1, c = 0)
    ),
    "model_data is NULL; no design matrices are created. Pass model_data to expand parameters via formulas.")

    expect_warning(
    design <- create_parameter_design(param_info),
    "Some parameters have infinite bounds. Creating ranges based on current values and infinite_range_factor = 2")

    expect_true("c" %in% names(design))
    expect_true(all(design$c == 5))
  })

  # ==========================================================================
  # Test 5: Infinite Bounds Handling
  # ==========================================================================
  it("handles infinite lower bound", {
    expect_warning(
      param_info <- create_parameter_info(
      parameter_values = c(k = 2),
      lower_bounds = c(k = -Inf),
      upper_bounds = c(k = 10),
      estimate_flag = c(k = 1)
    ),
    "model_data is NULL; no design matrices are created. Pass model_data to expand parameters via formulas.")

    expect_warning(
      design <- create_parameter_design(param_info),
      "Some parameters have infinite bounds. Creating ranges based on current values and infinite_range_factor = 2"
    )

    bounds <- attr(design, "bounds_used")
    expect_false(is.infinite(bounds$lower))
  })

  it("handles infinite upper bound", {
    param_info <- create_parameter_info(
      parameter_values = c(k = 2),
      lower_bounds = c(k = 0),
      upper_bounds = c(k = Inf),
      estimate_flag = c(k = 1),
      model_data = data.frame(x=1)
    )

    expect_warning(
      design <- create_parameter_design(param_info),
      "infinite bounds"
    )

    bounds <- attr(design, "bounds_used")
    expect_false(is.infinite(bounds$upper))
  })

  it("handles both bounds infinite", {
    param_info <- create_parameter_info(
      parameter_values = c(k = 2),
      lower_bounds = c(k = -Inf),
      upper_bounds = c(k = Inf),
      estimate_flag = c(k = 1),
      model_data=data.frame(x=1)
    )

    expect_warning(
      design <- create_parameter_design(param_info, infinite_range_factor = 3),
      "Some parameters have infinite bounds. Creating ranges based on current values and infinite_range_factor = 3"
    )

    bounds <- attr(design, "bounds_used")

    # Should create range [2/3, 2*3] = [0.667, 6]
    expect_equal(bounds$lower, 2/3, tolerance = 1e-6)
    expect_equal(bounds$upper, 2*3, tolerance = 1e-6)
  })

  # ==========================================================================
  # Test 6: Parameter Expansion (Formulas)
  # ==========================================================================

  it("works with expanded parameters from formulas", {
    data <- data.frame(group = factor(c("A", "B", "A", "B")))

    param_info <- create_parameter_info(
      parameter_values = c(k = 2, lambda = 10),
      formulas = list(k = ~0 + group, lambda = ~1),
      lower_bounds = c(k = 1, lambda = 5),
      upper_bounds = c(k = 4, lambda = 15),
      estimate_flag = c(k = 1, lambda = 1),
      model_data = data
    )

    design <- create_parameter_design(param_info)

    # Should have kA, kB, lambda columns
    expect_true("k.groupA" %in% names(design))
    expect_true("k.groupB" %in% names(design))
    expect_true("lambda.Intercept" %in% names(design))

    # 3 parameters in design
    expect_equal(attr(design, "num_parameters"), 3)
  })

  # ==========================================================================
  # Test 7: Input Validation
  # ==========================================================================

  it("rejects non-parameter_info input", {
    expect_error(
      create_parameter_design(data.frame(k = 1)),
      "must be an object of class 'parameter_info'"
    )
  })

  it("rejects when no parameters to estimate", {
    param_info <- create_parameter_info(
      parameter_values = c(k = 2, lambda = 10),
      estimate_flag = c(k = 0, lambda = 0),  # All fixed
      model_data=data.frame(x=1)
    )

    expect_error(
      create_parameter_design(param_info),
      "No parameters marked for estimation"
    )
  })

  it("rejects when lower >= upper", {
    param_info <- create_parameter_info(
      parameter_values = c(k = 2),
      lower_bounds = c(k = 5),
      upper_bounds = c(k = 1),  # Upper < lower
      estimate_flag = c(k = 1),
      model_data = data.frame(x=1)

    )

    expect_error(
      create_parameter_design(param_info),
      "Lower bounds must be less than upper"
    )
  })

  # ==========================================================================
  # Test 8: Design Attributes
  # ==========================================================================

  it("sets correct design attributes", {
    expect_warning(
      param_info <- create_parameter_info(
        parameter_values = c(k = 2, lambda = 10),
        lower_bounds = c(k = 1, lambda = 5),
        upper_bounds = c(k = 4, lambda = 15),
        estimate_flag = c(k = 1, lambda = 1)
      ),
      "model_data is NULL; no design matrices are created. Pass model_data to expand parameters via formulas.")

    design <- create_parameter_design(param_info)

    expect_equal(attr(design, "design_type"), "parameter_space_ccd")
    expect_equal(attr(design, "num_parameters"), 2)
    expect_equal(attr(design, "parameter_names"), c("k", "lambda"))
    expect_true(!is.null(attr(design, "bounds_used")))
  })

  it("stores bounds information", {
    expect_warning(
    param_info <- create_parameter_info(
      parameter_values = c(k = 2, lambda = 10),
      lower_bounds = c(k = 1, lambda = 5),
      upper_bounds = c(k = 4, lambda = 15),
      estimate_flag = c(k = 1, lambda = 1)
    ),
    "model_data is NULL; no design matrices are created. Pass model_data to expand parameters via formulas.")


    design <- create_parameter_design(param_info)
    bounds <- attr(design, "bounds_used")

    expect_s3_class(bounds, "data.frame")
    expect_equal(nrow(bounds), 2)
    expect_true("lower" %in% names(bounds))
    expect_true("upper" %in% names(bounds))
    expect_true("center" %in% names(bounds))
  })

  # ==========================================================================
  # Test 9: Alpha Options
  # ==========================================================================

  it("creates rotatable design", {
    expect_warning(
    param_info <- create_parameter_info(
      parameter_values = c(k = 2, lambda = 10),
      lower_bounds = c(k = 1, lambda = 5),
      upper_bounds = c(k = 4, lambda = 15),
      estimate_flag = c(k = 1, lambda = 1)
    ),
    "model_data is NULL; no design matrices are created. Pass model_data to expand parameters via formulas.")


    design <- create_parameter_design(param_info, alpha = "rotatable")

    expect_true(attr(design, "rotatable"))
  })

  it("creates face-centered design", {
    expect_warning(
    param_info <- create_parameter_info(
      parameter_values = c(k = 2, lambda = 10),
      lower_bounds = c(k = 1, lambda = 5),
      upper_bounds = c(k = 4, lambda = 15),
      estimate_flag = c(k = 1, lambda = 1)
    ),
    "model_data is NULL; no design matrices are created. Pass model_data to expand parameters via formulas.")


    design <- create_parameter_design(param_info, alpha = "faces")

    expect_equal(attr(design, "alpha"), 1)
  })

  it("accepts custom alpha", {
    expect_warning(
    param_info <- create_parameter_info(
      parameter_values = c(k = 2),
      lower_bounds = c(k = 1),
      upper_bounds = c(k = 4),
      estimate_flag = c(k = 1)
    ),
    "model_data is NULL; no design matrices are created. Pass model_data to expand parameters via formulas.")


    design <- create_parameter_design(param_info, alpha = 1.5)

    expect_equal(attr(design, "alpha"), 1.5)
  })

  # ==========================================================================
  # Test 10: Warnings
  # ==========================================================================

  it("warns about wide parameter ranges", {
    expect_warning(
    param_info <- create_parameter_info(
      parameter_values = c(k = 500),
      lower_bounds = c(k = 0),
      upper_bounds = c(k = 2000),  # Range = 2000
      estimate_flag = c(k = 1)
    ),
    "model_data is NULL; no design matrices are created. Pass model_data to expand parameters via formulas.")


    expect_warning(
      create_parameter_design(param_info),
      "Very wide parameter ranges"
    )
  })

  it("warns about narrow parameter ranges", {
    expect_warning(
    param_info <- create_parameter_info(
      parameter_values = c(k = 2),
      lower_bounds = c(k = 1.999),
      upper_bounds = c(k = 2.001),  # Range = 0.011
      estimate_flag = c(k = 1)
    ),
    "model_data is NULL; no design matrices are created. Pass model_data to expand parameters via formulas.")


    expect_warning(
      create_parameter_design(param_info),
      "Very narrow parameter ranges for: k. Consider fixing these parameters."
    )
  })

  # ==========================================================================
  # Test 11: Column Order
  # ==========================================================================

  it("preserves parameter order in output", {
    expect_warning(
    param_info <- create_parameter_info(
      parameter_values = c(k = 2, lambda = 10, c = 0),
      estimate_flag = c(k = 1, lambda = 1, c = 0)
    ),
    "model_data is NULL; no design matrices are created. Pass model_data to expand parameters via formulas.")

    expect_warning(
    design <- create_parameter_design(param_info),
    "Some parameters have infinite bounds. Creating ranges based on current values and infinite_range_factor = 2")

    param_cols <- names(design)[names(design) != "point_type"]
    expect_equal(param_cols, c("k", "lambda", "c"))
  })

  # ==========================================================================
  # Test 12: Data Integrity
  # ==========================================================================

  it("has no missing values", {
    expect_warning(
    param_info <- create_parameter_info(
      parameter_values = c(k = 2, lambda = 10),
      lower_bounds = c(k = 1, lambda = 5),
      upper_bounds = c(k = 4, lambda = 15),
      estimate_flag = c(k = 1, lambda = 1)
    ),
    "model_data is NULL; no design matrices are created. Pass model_data to expand parameters via formulas.")


    design <- create_parameter_design(param_info)

    expect_false(anyNA(design))
  })

  it("values are within specified bounds", {
    expect_warning(
    param_info <- create_parameter_info(
      parameter_values = c(k = 2, lambda = 10),
      lower_bounds = c(k = 1, lambda = 5),
      upper_bounds = c(k = 4, lambda = 15),
      estimate_flag = c(k = 1, lambda = 1)
    ),
    "model_data is NULL; no design matrices are created. Pass model_data to expand parameters via formulas.")


    design <- create_parameter_design(param_info, alpha = "faces")

    # For face-centered, all points should be within bounds
    # (for rotatable, star points may exceed bounds)
    expect_true(all(design$k >= 1 & design$k <= 4))
    expect_true(all(design$lambda >= 5 & design$lambda <= 15))
  })
})
