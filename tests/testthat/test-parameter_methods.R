
# Note: These tests assume the bstatUtils package is available with
# update_by_name() and assign_list_circular() functions

describe(\"create_parameter_info\", {

  # ============================================================================
  # Test 1: Basic parameter creation with defaults
  # ============================================================================
  #
  it(\"creates basic parameter_info object with minimal inputs\", {
    pars <- c(k = 2, lambda = 10)
    info <- create_parameter_info(parameter_values = pars)

    expect_s3_class(info, \"parameter_info\")
    expect_named(info, c(
      \"parameter_values\", \"lower_bounds\", \"upper_bounds\",
      \"estimate_flag\", \"link_functions\", \"formulas\",
      \"design_matrices\", \"model_data\", \"sep\"
    ))
    expect_equal(names(info$parameter_values), c(\"k\", \"lambda\"))
    expect_equal(length(info$parameter_values), 2)
  })

  # ============================================================================
  # Test 2: Parameter values are properly stored and named
  # ============================================================================
  #
  it(\"preserves parameter names and values\", {
    pars <- c(k = 2.5, lambda = 10.3, c = 0.1)
    info <- create_parameter_info(parameter_values = pars)

    expect_equal(info$parameter_values$k[[1]], 2.5)
    expect_equal(info$parameter_values$lambda[[1]], 10.3)
    expect_equal(info$parameter_values$c[[1]], 0.1)
  })

  # ============================================================================
  # Test 3: Default bounds are -Inf and Inf
  # ============================================================================
  #
  it(\"applies default bounds (±Inf) when not specified\", {
    pars <- c(k = 2, lambda = 10)
    info <- create_parameter_info(parameter_values = pars)

    expect_true(all(sapply(info$lower_bounds, function(x) all(is.infinite(x) & x < 0))))
    expect_true(all(sapply(info$upper_bounds, function(x) all(is.infinite(x) & x > 0))))
  })

  # ============================================================================
  # Test 4: Custom bounds are applied correctly
  # ============================================================================
  #
  it(\"applies custom bounds correctly\", {
    pars <- c(k = 2, lambda = 10, c = 0)
    info <- create_parameter_info(
      parameter_values = pars,
      lower_bounds = c(k = 0, lambda = 0, c = -1),
      upper_bounds = c(k = Inf, lambda = 100, c = 1)
    )

    expect_equal(info$lower_bounds$k, 0)
    expect_equal(info$lower_bounds$lambda, 0)
    expect_equal(info$lower_bounds$c, -1)
    expect_equal(info$upper_bounds$k, Inf)
    expect_equal(info$upper_bounds$lambda, 100)
    expect_equal(info$upper_bounds$c, 1)
  })

  # ============================================================================
  # Test 5: Estimate flags default to 1 (estimate all)
  # ============================================================================
  #
  it(\"defaults to estimate all parameters (flag = 1)\", {
    pars <- c(k = 2, lambda = 10)
    info <- create_parameter_info(parameter_values = pars)

    expect_true(all(sapply(info$estimate_flag, function(x) all(x == 1))))
  })

  # ============================================================================
  # Test 6: Estimate flags can be set to fix parameters
  # ============================================================================
  #
  it(\"correctly sets estimate flags (1 for estimate, 0 for fix)\", {
    pars <- c(k = 2, lambda = 10, c = 0)
    info <- create_parameter_info(
      parameter_values = pars,
      estimate_flag = c(k = 1, lambda = 1, c = 0)
    )

    expect_equal(info$estimate_flag$k, 1)
    expect_equal(info$estimate_flag$lambda, 1)
    expect_equal(info$estimate_flag$c, 0)
  })

  # ============================================================================
  # Test 7: Default formulas are intercept-only (~1)
  # ============================================================================
  #
  it(\"defaults to intercept-only formulas (~1)\", {
    pars <- c(k = 2, lambda = 10)
    info <- create_parameter_info(parameter_values = pars)

    expect_true(all(sapply(info$formulas, function(f) {
      identical(f, ~1)
    })))
  })

  # ============================================================================
  # Test 8: Custom formulas are applied correctly
  # ============================================================================
  #
  it(\"applies custom formulas to each parameter\", {
    pars <- c(k = 2, lambda = 10, c = 0)
    df <- expand.grid(Group = c(\"A\", \"B\"), Rep = c(1, 2))

    info <- create_parameter_info(
      parameter_values = pars,
      formulas = list(k = ~Group, lambda = ~1, c = ~Group:Rep),
      model_data = df
    )

    expect_equal(info$formulas$k, ~Group)
    expect_equal(info$formulas$lambda, ~1)
    expect_equal(info$formulas$c, ~Group:Rep)
  })

  # ============================================================================
  # Test 9: Design matrices are created from formulas
  # ============================================================================
  #
  it(\"generates design matrices from formulas and model_data\", {
    pars <- c(k = 2, lambda = 10)
    df <- expand.grid(Group = c(\"A\", \"B\"))

    info <- create_parameter_info(
      parameter_values = pars,
      formulas = list(k = ~Group, lambda = ~1),
      model_data = df
    )

    # k formula ~Group should produce 2 columns
    expect_equal(ncol(info$design_matrices$k), 2)
    # lambda formula ~1 should produce 1 column
    expect_equal(ncol(info$design_matrices$lambda), 1)
  })

  # ============================================================================
  # Test 10: Parameter values are expanded to match design matrix columns
  # ============================================================================
  #
  it(\"expands parameter values to match design matrix dimensions\", {
    pars <- c(k = 2, lambda = 10)
    df <- expand.grid(Group = c(\"A\", \"B\"))

    info <- create_parameter_info(
      parameter_values = pars,
      formulas = list(k = ~0 + Group, lambda = ~1),
      model_data = df
    )

    # k has 2 design matrix columns, so parameter should be length 2
    expect_equal(length(info$parameter_values$k), 2)
    # lambda has 1 column, so parameter should be length 1
    expect_equal(length(info$parameter_values$lambda), 1)
  })

  # ============================================================================
  # Test 11: Default link functions are 'identity'
  # ============================================================================
  #
  it(\"defaults to identity link functions\", {
    pars <- c(k = 2, lambda = 10)
    info <- create_parameter_info(parameter_values = pars)

    expect_true(all(sapply(info$link_functions, function(x) x == \"identity\")))
  })

  # ============================================================================
  # Test 12: Custom link functions are applied
  # ============================================================================
  #
  it(\"applies custom link functions\", {
    pars <- c(k = 2, lambda = 10, c = 0)
    info <- create_parameter_info(
      parameter_values = pars,
      link_functions = c(\"log\", \"log\", \"identity\")
    )

    expect_equal(info$link_functions$k, \"log\")
    expect_equal(info$link_functions$lambda, \"log\")
    expect_equal(info$link_functions$c, \"identity\")
  })

  # ============================================================================
  # Test 13: Error if parameter_values is not named
  # ============================================================================
  #
  it(\"raises error if parameter_values are not named\", {
    expect_error(
      create_parameter_info(parameter_values = c(2, 10)),
      \"parameter_values must be a named vector\"
    )
  })

  # ============================================================================
  # Test 14: Error if parameter names are duplicated
  # ============================================================================
  #
  it(\"raises error if parameter names are duplicated\", {
    expect_error(
      create_parameter_info(parameter_values = c(k = 2, k = 10)),
      \"parameter_values must be a named vector\"
    )
  })

  # ============================================================================
  # Test 15: Error if argument lengths don't match number of parameters
  # ============================================================================
  #
  it(\"raises error if estimate_flag length doesn't match parameters\", {
    pars <- c(k = 2, lambda = 10, c = 0)
    expect_error(
      create_parameter_info(
        parameter_values = pars,
        estimate_flag = c(1, 0)  # Only 2 values for 3 parameters
      ),
      \"not equal to the number of parameters\"
    )
  })

  # ============================================================================
  # Test 16: Compositional update - preserve unmodified components
  # ============================================================================
  #
  it(\"preserves unmodified components when updating existing parameter_info\", {
    pars <- c(k = 2, lambda = 10)
    info1 <- create_parameter_info(
      parameter_values = pars,
      estimate_flag = c(k = 1, lambda = 0),
      lower_bounds = c(k = 0, lambda = 1)
    )

    # Update only the upper bounds
    info2 <- create_parameter_info(
      parameter_values = info1,
      upper_bounds = c(k = 100, lambda = 200)
    )

    # Check that other components are preserved
    expect_equal(info2$estimate_flag$k, 1)
    expect_equal(info2$estimate_flag$lambda, 0)
    expect_equal(info2$lower_bounds$k, 0)
    expect_equal(info2$lower_bounds$lambda, 1)
    # Check that upper bounds are updated
    expect_equal(info2$upper_bounds$k, 100)
    expect_equal(info2$upper_bounds$lambda, 200)
  })

  # ============================================================================
  # Test 17: Compositional update - modify formula
  # ============================================================================
  #
  it(\"allows compositional formula updates\", {
    pars <- c(k = 2, lambda = 10)
    df <- expand.grid(Group = c(\"A\", \"B\"))

    info1 <- create_parameter_info(
      parameter_values = pars,
      formulas = list(k = ~1, lambda = ~1),
      model_data = df
    )

    # Update only the k formula
    info2 <- create_parameter_info(
      parameter_values = info1,
      formulas = list(k = ~Group)
    )

    expect_equal(info2$formulas$k, ~Group)
    expect_equal(info2$formulas$lambda, ~1)
  })

  # ============================================================================
  # Test 18: Warning when model_data is NULL
  # ============================================================================
  #
  it(\"warns when model_data is NULL\", {
    pars <- c(k = 2, lambda = 10)
    expect_warning(
      create_parameter_info(parameter_values = pars, model_data = NULL),
      \"model_data is NULL\"
    )
  })

  # ============================================================================
  # Test 19: Design matrices are NULL when model_data is NULL
  # ============================================================================
  #
  it(\"creates NULL design matrices when model_data is NULL\", {
    pars <- c(k = 2, lambda = 10)
    info <- create_parameter_info(parameter_values = pars, model_data = NULL)

    expect_true(all(sapply(info$design_matrices, is.null)))
  })

  # ============================================================================
  # Test 20: Separator is correctly stored
  # ============================================================================
  #
  it(\"stores custom separator in parameter_info object\", {
    pars <- c(k = 2, lambda = 10)
    info <- create_parameter_info(parameter_values = pars, sep = \"_\")

    expect_equal(info$sep, \"_\")
  })
})

# ============================================================================
# Summary method tests
# ============================================================================
#
describe(\"summary.parameter_info\", {

  # ============================================================================
  # Test 21: Summary returns list with info and parameters
  # ============================================================================
  #
  it(\"returns a list with 'info' and 'parameters' components\", {
    pars <- c(k = 2, lambda = 10)
    info <- create_parameter_info(parameter_values = pars)
    result <- summary(info)

    expect_type(result, \"list\")
    expect_named(result, c(\"info\", \"parameters\"))
  })

  # ============================================================================
  # Test 22: Summary info includes formulas with link functions
  # ============================================================================
  #
  it(\"includes formulas with link functions in info\", {
    pars <- c(k = 2, lambda = 10)
    info <- create_parameter_info(
      parameter_values = pars,
      link_functions = c(\"log\", \"identity\")
    )
    result <- summary(info)

    expect_true(any(grepl(\"log\", result$info$formula)))
    expect_true(any(grepl(\"identity\", result$info$formula)))
  })

  # ============================================================================
  # Test 23: Summary parameters include all components
  # ============================================================================
  #
  it(\"includes all parameter components in parameters table\", {
    pars <- c(k = 2, lambda = 10)
    info <- create_parameter_info(
      parameter_values = pars,
      estimate_flag = c(k = 1, lambda = 0),
      lower_bounds = c(k = 0, lambda = 1),
      upper_bounds = c(k = 100, lambda = 200)
    )
    result <- summary(info)

    params <- result$parameters
    expect_true(\"parameter_values\" %in% colnames(params))
    expect_true(\"estimate_flag\" %in% colnames(params))
    expect_true(\"lower_bounds\" %in% colnames(params))
    expect_true(\"upper_bounds\" %in% colnames(params))
  })

  # ============================================================================
  # Test 24: Summary with formulas shows expanded parameters
  # ============================================================================
  #
  it(\"shows expanded parameter names for design matrix columns\", {
    pars <- c(k = 2, lambda = 10)
    df <- expand.grid(Group = c(\"A\", \"B\"))

    info <- create_parameter_info(
      parameter_values = pars,
      formulas = list(k = ~0 + Group, lambda = ~1),
      model_data = df
    )
    result <- summary(info)

    # k should have 2 rows (GroupA and GroupB)
    k_rows <- result$parameters[result$parameters$parameter == \"k\", ]
    expect_equal(nrow(k_rows), 2)
  })

})

# ============================================================================
# Edge case tests
# ============================================================================
#
describe(\"Edge cases\", {

  # ============================================================================
  # Test 25: Single parameter
  # ============================================================================
  #
  it(\"handles single parameter correctly\", {
    info <- create_parameter_info(parameter_values = c(k = 2))

    expect_equal(length(info$parameter_values), 1)
    expect_named(info$parameter_values, \"k\")
  })

  # ============================================================================
  # Test 26: Many parameters
  # ============================================================================
  #
  it(\"handles many parameters correctly\", {
    pars <- setNames(1:10, paste0(\"p\", 1:10))
    info <- create_parameter_info(parameter_values = pars)

    expect_equal(length(info$parameter_values), 10)
    expect_equal(names(info$parameter_values), paste0(\"p\", 1:10))
  })

  # ============================================================================
  # Test 27: Complex formula with interactions
  # ============================================================================
  #
  it(\"handles complex formulas with interactions\", {
    pars <- c(k = 2)
    df <- expand.grid(A = c(\"a1\", \"a2\"), B = c(\"b1\", \"b2\"))

    info <- create_parameter_info(
      parameter_values = pars,
      formulas = list(k = ~A * B),
      model_data = df
    )

    # A * B interaction should produce 4 columns (2x2 design)
    expect_equal(ncol(info$design_matrices$k), 4)
  })

  # ============================================================================
  # Test 28: Parameter values recycled when design matrix larger
  # ============================================================================
  #
  it(\"recycles scalar parameter values to match design matrix size\", {
    pars <- c(k = 2)  # Single value
    df <- expand.grid(Group = c(\"A\", \"B\", \"C\"))

    info <- create_parameter_info(
      parameter_values = pars,
      formulas = list(k = ~0 + Group),
      model_data = df
    )

    # k parameter should be recycled to length 3
    expect_equal(length(info$parameter_values$k), 3)
    # All values should be 2
    expect_true(all(info$parameter_values$k == 2))
  })

  # ============================================================================
  # Test 29: Zero bounds are valid
  # ============================================================================
  #
  it(\"allows zero as a bound value\", {
    pars <- c(k = 2, lambda = 10)
    info <- create_parameter_info(
      parameter_values = pars,
      lower_bounds = c(k = 0, lambda = 0),
      upper_bounds = c(k = 0, lambda = Inf)
    )

    expect_equal(info$lower_bounds$k, 0)
    expect_equal(info$upper_bounds$k, 0)
  })

  # ============================================================================
  # Test 30: Negative parameter values
  # ============================================================================
  #
  it(\"accepts negative parameter values\", {
    pars <- c(k = -2.5, lambda = -10)
    info <- create_parameter_info(parameter_values = pars)

    expect_equal(info$parameter_values$k, -2.5)
    expect_equal(info$parameter_values$lambda, -10)
  })
})
