
describe("augment_with_star_center_points", {

  # ==========================================================================
  # Test 1: Basic Functionality
  # ==========================================================================

  it("augments 2^3 design with star and center points", {
    factorial <- fractional_factorial_design(num_factors = 3)
    ccd <- augment_with_star_center_points(factorial)

    # Should have factorial + axial + center points
    expect_equal(nrow(ccd), 8 + 6 + 5)  # 8 factorial + 6 axial + 5 center
    expect_equal(ncol(ccd), 4)  # 3 factors + point_type
    expect_true("point_type" %in% names(ccd))
  })

  it("creates correct number of axial points", {
    factorial <- fractional_factorial_design(num_factors = 4)
    ccd <- augment_with_star_center_points(factorial, center_points = 5)

    # Axial points = 2 * num_factors
    expect_equal(sum(ccd$point_type == "axial"), 8)  # 2 * 4
  })

  it("creates specified number of center points", {
    factorial <- fractional_factorial_design(num_factors = 2)
    ccd <- augment_with_star_center_points(factorial, center_points = 10)

    expect_equal(sum(ccd$point_type == "center"), 10)
  })

  # ==========================================================================
  # Test 2: Rotatable Designs
  # ==========================================================================

  it("computes rotatable alpha correctly for 2^2 design", {
    factorial <- fractional_factorial_design(num_factors = 2)
    ccd <- augment_with_star_center_points(factorial, alpha = "rotatable")

    # For 4 factorial points: alpha = 4^0.25 = sqrt(2) ≈ 1.414
    expected_alpha <- 4^0.25
    expect_equal(attr(ccd, "alpha"), expected_alpha, tolerance = 1e-6)
    expect_true(attr(ccd, "rotatable"))
  })

  it("computes rotatable alpha correctly for 2^3 design", {
    factorial <- fractional_factorial_design(num_factors = 3)
    ccd <- augment_with_star_center_points(factorial, alpha = "rotatable")

    # For 8 factorial points: alpha = 8^0.25 ≈ 1.682
    expected_alpha <- 8^0.25
    expect_equal(attr(ccd, "alpha"), expected_alpha, tolerance = 1e-6)
  })

  it("computes rotatable alpha for fractional design", {
    factorial <- fractional_factorial_design(
      num_factors = 4,
      generators = "D=ABC"
    )
    ccd <- augment_with_star_center_points(factorial, alpha = "rotatable")

    # 8 factorial points: alpha = 8^0.25 ≈ 1.682
    expected_alpha <- 8^0.25
    expect_equal(attr(ccd, "alpha"), expected_alpha, tolerance = 1e-6)
  })

  # ==========================================================================
  # Test 3: Face-Centered Designs
  # ==========================================================================

  it("creates face-centered design with alpha = 1", {
    factorial <- fractional_factorial_design(num_factors = 3)
    ccf <- augment_with_star_center_points(factorial, alpha = "faces")

    expect_equal(attr(ccf, "alpha"), 1)
    expect_false(attr(ccf, "rotatable"))
  })

  it("handles 'face' and 'face-centered' alpha options", {
    factorial <- fractional_factorial_design(num_factors = 2)

    ccf1 <- augment_with_star_center_points(factorial, alpha = "face")
    ccf2 <- augment_with_star_center_points(factorial, alpha = "face-centered")

    expect_equal(attr(ccf1, "alpha"), 1)
    expect_equal(attr(ccf2, "alpha"), 1)
  })

  # ==========================================================================
  # Test 4: Custom Alpha Values
  # ==========================================================================

  it("accepts custom numeric alpha", {
    factorial <- fractional_factorial_design(num_factors = 3)
    ccd <- augment_with_star_center_points(factorial, alpha = 1.5)

    expect_equal(attr(ccd, "alpha"), 1.5)
  })

  it("rejects non-positive alpha", {
    factorial <- fractional_factorial_design(num_factors = 2)

    expect_error(
      augment_with_star_center_points(factorial, alpha = 0),
      "must be positive"
    )

    expect_error(
      augment_with_star_center_points(factorial, alpha = -1),
      "must be positive"
    )
  })

  # ==========================================================================
  # Test 5: Axial Points Structure
  # ==========================================================================

  it("creates axial points at correct positions", {
    factorial <- fractional_factorial_design(num_factors = 2)
    ccd <- augment_with_star_center_points(factorial, alpha = 2)

    axial <- ccd[ccd$point_type == "axial", c("A", "B")]

    # Should have 4 axial points for 2 factors
    expect_equal(nrow(axial), 4)

    # Check specific axial points
    # (+/-alpha, 0) and (0, +/-alpha)
    expected_axial <- matrix(
      c(2, 0, -2, 0, 0, 2, 0, -2),
      nrow = 4, ncol = 2, byrow = TRUE
    )

    # Sort both for comparison
    axial_sorted <- axial[order(axial$A, axial$B), ]
    expected_sorted <- expected_axial[order(expected_axial[,1], expected_axial[,2]), ]
    colnames(expected_sorted) <- c("A","B")
    rownames(expected_sorted) <- c(6,8,7,5)

    expect_equal(as.matrix(axial_sorted), expected_sorted, tolerance = 1e-10)
  })

  it("axial points have zeros in all but one position", {
    factorial <- fractional_factorial_design(num_factors = 4)
    ccd <- augment_with_star_center_points(factorial, alpha = "rotatable")

    axial <- ccd[ccd$point_type == "axial", 1:4]

    # Each axial point should have exactly 3 zeros (for 4 factors)
    zeros_per_row <- apply(axial, 1, function(row) sum(abs(row) < 1e-10))
    expect_true(all(zeros_per_row == 3))
  })

  # ==========================================================================
  # Test 6: Center Points Structure
  # ==========================================================================

  it("center points are all zeros", {
    factorial <- fractional_factorial_design(num_factors = 3)
    ccd <- augment_with_star_center_points(factorial, center_points = 5)

    center <- ccd[ccd$point_type == "center", c("A", "B", "C")]

    # All should be zero
    expect_true(all(abs(center) < 1e-10))
    expect_equal(nrow(center), 5)
  })

  # ==========================================================================
  # Test 7: Blocking
  # ==========================================================================

  it("adds block column when requested", {
    factorial <- fractional_factorial_design(num_factors = 2)
    ccd <- augment_with_star_center_points(factorial, add_blocks = TRUE)

    expect_true("block" %in% names(ccd))
    expect_equal(length(unique(ccd$block)), 3)  # 3 blocks
  })

  it("block assignments are correct", {
    factorial <- fractional_factorial_design(num_factors = 2)
    ccd <- augment_with_star_center_points(factorial, add_blocks = TRUE)

    # Factorial in block 1
    expect_true(all(ccd$block[ccd$point_type == "factorial"] == "1"))

    # Axial in block 2
    expect_true(all(ccd$block[ccd$point_type == "axial"] == "2"))

    # Center in block 3
    expect_true(all(ccd$block[ccd$point_type == "center"] == "3"))
  })

  # ==========================================================================
  # Test 8: Input Validation
  # ==========================================================================

  it("rejects non-data.frame input", {
    expect_error(
      augment_with_star_center_points(c(1, 2, 3)),
      "must be a data frame"
    )
  })

  it("rejects empty data frame", {
    expect_error(
      augment_with_star_center_points(data.frame()),
      "is empty"
    )
  })

  it("rejects non-numeric columns", {
    bad_design <- data.frame(
      A = c(1, -1),
      B = c("high", "low")
    )

    expect_error(
      augment_with_star_center_points(bad_design),
      "only numeric columns"
    )
  })

  it("rejects invalid alpha strings", {
    factorial <- fractional_factorial_design(num_factors = 2)

    expect_error(
      augment_with_star_center_points(factorial, alpha = "invalid"),
      "must be 'rotatable', 'orthogonal', 'faces'"
    )
  })

  it("rejects non-positive center points", {
    factorial <- fractional_factorial_design(num_factors = 2)

    expect_error(
      augment_with_star_center_points(factorial, center_points = 0),
      "must be a positive integer"
    )

    expect_error(
      augment_with_star_center_points(factorial, center_points = -5),
      "must be a positive integer"
    )
  })

  it("rejects non-integer center points", {
    factorial <- fractional_factorial_design(num_factors = 2)

    expect_error(
      augment_with_star_center_points(factorial, center_points = 3.5),
      "must be a positive integer"
    )
  })

  # ==========================================================================
  # Test 9: Design Attributes
  # ==========================================================================

  it("sets correct design attributes", {
    factorial <- fractional_factorial_design(num_factors = 3)
    ccd <- augment_with_star_center_points(factorial, center_points = 5)

    expect_equal(attr(ccd, "design_type"), "central_composite")
    expect_equal(attr(ccd, "num_factors"), 3)
    expect_equal(attr(ccd, "num_factorial"), 8)
    expect_equal(attr(ccd, "num_axial"), 6)
    expect_equal(attr(ccd, "num_center"), 5)
  })

  # ==========================================================================
  # Test 10: Point Type Column
  # ==========================================================================

  it("creates point_type as factor", {
    factorial <- fractional_factorial_design(num_factors = 2)
    ccd <- augment_with_star_center_points(factorial)

    expect_s3_class(ccd$point_type, "factor")
    expect_equal(levels(ccd$point_type), c("factorial", "axial", "center"))
  })

  it("point_type counts are correct", {
    factorial <- fractional_factorial_design(num_factors = 3)
    ccd <- augment_with_star_center_points(factorial, center_points = 7)

    type_counts <- table(ccd$point_type)

    expect_equal(type_counts["factorial"], c(factorial=8))
    expect_equal(type_counts["axial"], c(axial=6))
    expect_equal(type_counts["center"], c(center=7))
  })

  # ==========================================================================
  # Test 11: Integration with Fractional Designs
  # ==========================================================================

  it("works with 2^(4-1) fractional design", {
    factorial <- fractional_factorial_design(
      num_factors = 4,
      generators = "D=ABC"
    )
    ccd <- augment_with_star_center_points(factorial)

    # 8 factorial + 8 axial + center
    expect_true(nrow(ccd) >= 16)
    expect_equal(sum(ccd$point_type == "factorial"), 8)
    expect_equal(sum(ccd$point_type == "axial"), 8)
  })

  it("works with 2^(5-2) fractional design", {
    expect_warning(
      factorial <- fractional_factorial_design(
        num_factors = 5,
        generators = c("D=ABC", "E=BCD")
      ),
      "Confounding detected: Main effect A is confounded with main effect E")
    ccd <- augment_with_star_center_points(factorial)

    # Should have 10 axial points (2 * 5 factors)
    expect_equal(sum(ccd$point_type == "axial"), 10)
  })

  # ==========================================================================
  # Test 12: Design Size Calculations
  # ==========================================================================

  it("total design size is correct", {
    factorial <- fractional_factorial_design(num_factors = 2)
    ccd <- augment_with_star_center_points(factorial, center_points = 5)

    # 4 factorial + 4 axial + 5 center = 13
    expect_equal(nrow(ccd), 13)
  })

  it("matches standard CCD sizes", {
    # Standard 2^3 CCD with 6 center points
    factorial <- fractional_factorial_design(num_factors = 3)
    ccd <- augment_with_star_center_points(factorial, center_points = 6)

    # 8 + 6 + 6 = 20 (standard CCD for k=3)
    expect_equal(nrow(ccd), 20)
  })

  # ==========================================================================
  # Test 13: Data Integrity
  # ==========================================================================

  it("preserves original factorial points", {
    factorial <- fractional_factorial_design(num_factors = 2)
    ccd <- augment_with_star_center_points(factorial)

    # Extract factorial points from CCD
    factorial_in_ccd <- ccd[ccd$point_type == "factorial", c("A", "B")]

    # Should match original
    expect_equal(
      factorial_in_ccd,
      factorial[, c("A", "B")],
      ignore_attr = FALSE
    )
  })

  it("has no missing values", {
    factorial <- fractional_factorial_design(num_factors = 4)
    ccd <- augment_with_star_center_points(factorial)

    expect_false(anyNA(ccd[, 1:4]))  # Factor columns
  })

  # ==========================================================================
  # Test 14: Default Center Points
  # ==========================================================================

  it("uses sensible default for center points", {
    factorial <- fractional_factorial_design(num_factors = 2)
    ccd <- augment_with_star_center_points(factorial)

    # Should default to 4-6 center points
    n_center <- sum(ccd$point_type == "center")
    expect_true(n_center >= 4 && n_center <= 6)
  })

  it("default center points scales with num_factors", {
    factorial2 <- fractional_factorial_design(num_factors = 2)
    factorial5 <- fractional_factorial_design(num_factors = 5)

    ccd2 <- augment_with_star_center_points(factorial2)
    ccd5 <- augment_with_star_center_points(factorial5)

    n_center2 <- sum(ccd2$point_type == "center")
    n_center5 <- sum(ccd5$point_type == "center")

    # More factors should have more (or equal) center points
    expect_true(n_center5 >= n_center2)
  })

  # ==========================================================================
  # Test 15: Messages
  # ==========================================================================

  it("displays informative messages", {
    factorial <- fractional_factorial_design(num_factors = 3)

    expect_message(
      augment_with_star_center_points(factorial),
      "Central Composite Design Summary"
    )
  })
})
