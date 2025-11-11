
describe("fractional_factorial_design", {

  # ==========================================================================
  # Test 1: Full Factorial Designs
  # ==========================================================================

  it("creates 2^3 full factorial design (8 runs)", {
    design <- fractional_factorial_design(num_factors = 3)

    expect_s3_class(design, "data.frame")
    expect_equal(nrow(design), 8)  # 2^3
    expect_equal(ncol(design), 3)  # 3 factors
    expect_named(design, c("A", "B", "C"))
    expect_true(all(unlist(design) %in% c(-1, 1)))
  })

  it("creates 2^4 full factorial design (16 runs)", {
    design <- fractional_factorial_design(num_factors = 4)

    expect_equal(nrow(design), 16)  # 2^4
    expect_equal(ncol(design), 4)
    expect_named(design, c("A", "B", "C", "D"))
  })

  it("creates 2^2 minimal full factorial (4 runs)", {
    design <- fractional_factorial_design(num_factors = 2)

    expect_equal(nrow(design), 4)  # 2^2
    expect_equal(ncol(design), 2)
    expect_named(design, c("A", "B"))
  })

  # ==========================================================================
  # Test 2: Fractional Factorial Designs
  # ==========================================================================

  it("creates 2^(4-1) fractional design with generator D=ABC", {
    design <- fractional_factorial_design(
      num_factors = 4,
      generators = "D=ABC"
    )

    expect_equal(nrow(design), 8)   # 2^(4-1) = 8 runs
    expect_equal(ncol(design), 4)   # 4 factors
    expect_named(design, c("A", "B", "C", "D"))

    # Check that D = A*B*C
    expect_equal(
      design$D,
      design$A * design$B * design$C
    )
  })

  it("creates 2^(5-1) fractional design", {
    design <- fractional_factorial_design(
      num_factors = 5,
      generators = "E=ABCD"
    )

    expect_equal(nrow(design), 16)  # 2^(5-1) = 16 runs
    expect_equal(ncol(design), 5)

    # E should equal product of A,B,C,D
    expect_equal(
      design$E,
      design$A * design$B * design$C * design$D
    )
  })

  it("creates 2^(6-2) fractional design with two generators", {
    design <- fractional_factorial_design(
      num_factors = 6,
      generators = c("E=ABC", "F=BCD")
    )

    expect_equal(nrow(design), 16)  # 2^(6-2) = 16 runs
    expect_equal(ncol(design), 6)

    # Check generators
    expect_equal(design$E, design$A * design$B * design$C)
    expect_equal(design$F, design$B * design$C * design$D)
  })

  # ==========================================================================
  # Test 3: Negative Generators
  # ==========================================================================

  it("handles negative generators correctly", {
    design <- fractional_factorial_design(
      num_factors = 4,
      generators = "D=-ABC"
    )

    expect_equal(nrow(design), 8)
    # D should equal -(A*B*C)
    expect_equal(
      design$D,
      -1 * design$A * design$B * design$C
    )
  })

  it("handles mix of positive and negative generators", {
    expect_warning(design <- fractional_factorial_design(
      num_factors = 5,
      generators = c("D=ABC", "E=-ABD")),
      "Confounding detected: Main effect -C is confounded with main effect E")

    expect_equal(nrow(design), 8)  # 2^(5-2) = 8 runs

    # Check generators
    expect_equal(design$D, design$A * design$B * design$C)
    expect_equal(design$E, -1 * design$A * design$B * design$D)
  })

  # ==========================================================================
  # Test 4: Input Validation
  # ==========================================================================

  it("rejects levels other than 2", {
    expect_error(
      fractional_factorial_design(num_factors = 3, levels = 3),
      "only 2-level designs are supported"
    )
  })

  it("rejects non-positive num_factors", {
    expect_error(
      fractional_factorial_design(num_factors = 0),
      "must be a positive integer"
    )

    expect_error(
      fractional_factorial_design(num_factors = -1),
      "must be a positive integer"
    )
  })

  it("rejects non-integer num_factors", {
    expect_error(
      fractional_factorial_design(num_factors = 3.5),
      "must be a positive integer"
    )
  })

  it("rejects when num_generators >= num_factors", {
    expect_error(
      fractional_factorial_design(
        num_factors = 3,
        generators = c("B=A", "C=A", "D=AB")
      ),
      "must be less than"
    )
  })

  # ==========================================================================
  # Test 5: Generator Validation
  # ==========================================================================

  it("rejects invalid generator with no RHS factors", {
    expect_error(
      fractional_factorial_design(
        num_factors = 4,
        generators = "D=E"  # E doesn't exist yet
      ),
      "not properly defined"
    )
  })

  it("rejects generator with only one RHS factor", {
    expect_error(
      fractional_factorial_design(
        num_factors = 3,
        generators = "C=A"  # Only one factor on RHS
      ),
      "not properly defined"
    )
  })

  # it("rejects generator with multiple LHS factors", {
  #   expect_error(
  #     fractional_factorial_design(
  #       num_factors = 5,
  #       generators = "DE=ABC"  # Two new factors on LHS
  #     ),
  #     "not properly defined"
  #   )
  # })

  # ==========================================================================
  # Test 6: Design Properties
  # ==========================================================================

  it("produces orthogonal columns in full factorial", {
    design <- fractional_factorial_design(num_factors = 3)

    # Check orthogonality: column dot products should be 0
    expect_equal(sum(design$A * design$B), 0)
    expect_equal(sum(design$A * design$C), 0)
    expect_equal(sum(design$B * design$C), 0)
  })

  it("produces balanced design (equal -1 and +1 counts)", {
    design <- fractional_factorial_design(num_factors = 4)

    for (col in names(design)) {
      expect_equal(sum(design[[col]] == -1), nrow(design) / 2)
      expect_equal(sum(design[[col]] == +1), nrow(design) / 2)
    }
  })

  it("produces orthogonal columns in fractional factorial", {
    design <- fractional_factorial_design(
      num_factors = 5,
      generators = "E=ABCD"
    )

    # Base factors should be orthogonal
    expect_equal(sum(design$A * design$B), 0)
    expect_equal(sum(design$A * design$C), 0)
  })

  # ==========================================================================
  # Test 7: Factor Naming
  # ==========================================================================

  it("uses uppercase letters for factor names", {
    design <- fractional_factorial_design(num_factors = 5)

    expect_named(design, c("A", "B", "C", "D", "E"))
  })

  it("extends to lowercase letters when needed", {
    design <- fractional_factorial_design(num_factors = 27)

    expect_true("a" %in% names(design))
    expect_equal(names(design)[27], "a")
  })

  # ==========================================================================
  # Test 8: Generator Parsing
  # ==========================================================================

  it("handles generators with spaces", {
    design <- fractional_factorial_design(
      num_factors = 4,
      generators = "D = A B C"  # Spaces included
    )

    expect_equal(nrow(design), 8)
    expect_equal(design$D, design$A * design$B * design$C)
  })

  it("handles generators without equals sign", {
    design <- fractional_factorial_design(
      num_factors = 4,
      generators = "DABC"  # No equals sign
    )

    expect_equal(nrow(design), 8)
  })

  it("handles lowercase generator input", {
    design <- fractional_factorial_design(
      num_factors = 4,
      generators = "d=abc"  # Lowercase
    )

    expect_equal(nrow(design), 8)
    expect_named(design, c("A", "B", "C", "D"))  # Still uppercase output
  })

  # ==========================================================================
  # Test 9: Confounding Detection
  # ==========================================================================

  it("warns about main effect confounding", {
    expect_warning(
      fractional_factorial_design(
        num_factors = 5,
        generators = c("D=ABC", "E=ABC")  # D and E confounded
      ),
      "Main effect.*confounded with main effect"
    )
  })

  # it("warns about interaction confounding", {
  #   expect_warning(
  #     fractional_factorial_design(
  #       num_factors = 6,
  #       generators = c("E=ABC", "F=ABD")  # Generates confounding
  #     ),
  #     "confounded"
  #   )
  # })

  # ==========================================================================
  # Test 10: Defining Relation Display
  # ==========================================================================

  it("displays defining relation for single generator", {
    expect_message(
      fractional_factorial_design(
        num_factors = 4,
        generators = "D=ABC"
      ),
      "Defining relation: I = ABCD"
    )
  })

  it("displays complete defining relation for multiple generators", {
    expect_warning(expect_message(
      design <- fractional_factorial_design(
        num_factors = 5,
        generators = c("D=ABC", "E=ABD")
      ),
      "Defining relation: I = ABCD = ABDE = CE"
    ),
    "Confounding detected: Main effect C is confounded with main effect E")
  })

  # ==========================================================================
  # Test 11: Common Standard Designs
  # ==========================================================================

  it("creates standard 2^(3-1) Resolution III design", {
    design <- fractional_factorial_design(
      num_factors = 3,
      generators = "C=AB"
    )

    expect_equal(nrow(design), 4)  # 2^(3-1) = 4 runs
    expect_equal(ncol(design), 3)
  })

  it("creates standard 2^(7-4) Resolution III design", {
    expect_message(
      expect_warning(
        expect_warning(
          expect_warning(
            design <- fractional_factorial_design(
            num_factors = 7,
            generators = c("D=AB", "E=AC", "F=BC", "G=ABC")
          ),
          "Confounding detected: Main effect C is confounded with interaction D:G"),
        "Confounding detected: Main effect B is confounded with interaction E:G"),
      "Confounding detected: Main effect A is confounded with interaction F:G"),
    "Defining relation: I = ABD = ACE = BCF = ABCG = BCDE = ACDF = CDG = ABEF = BEG = AFG")

    expect_equal(nrow(design), 8)  # 2^(7-4) = 8 runs
    expect_equal(ncol(design), 7)
  })

  # ==========================================================================
  # Test 12: Edge Cases
  # ==========================================================================

  it("handles single factor", {
    design <- fractional_factorial_design(num_factors = 1)

    expect_equal(nrow(design), 2)  # 2^1
    expect_equal(ncol(design), 1)
    expect_named(design, "A")
  })

  it("handles large design (2^8 = 256 runs)", {
    design <- fractional_factorial_design(num_factors = 8)

    expect_equal(nrow(design), 256)
    expect_equal(ncol(design), 8)
  })

  it("handles fractional design with many generators", {
    expect_message(
      expect_warning(
        expect_warning(
          expect_warning(
            expect_warning(design <- fractional_factorial_design(
              num_factors = 10,
              generators = c(
                "E=ABC",
                "F=BCD",
                "G=ACD",
                "H=ABD",
                "I=AB",
                "J=AC"
              )
            ),
          "Confounding detected: Main effect C is confounded with interaction E:I"
          ),
        "Confounding detected: Main effect B is confounded with interaction E:J"
        ),
        "Confounding detected: Main effect D is confounded with interaction G:J"),
      "Confounding detected: Main effect D is confounded with interaction H:I"),
    "Defining relation: I = ABCE = BCDF = ACDG = ABDH = ABI = ACJ = ADEF = BDEG = CDEH = CEI = BEJ = ABFG = ACFH = ACDFI = ABDFJ = BCGH = BCDGI = DGJ = DHI = BCDHJ = BCIJ")

    expect_equal(nrow(design), 16)  # 2^(10-6) = 16 runs
    expect_equal(ncol(design), 10)
  })

  # ==========================================================================
  # Test 13: Data Integrity
  # ==========================================================================

  it("contains only -1 and +1 values", {
    design <- fractional_factorial_design(
      num_factors = 5,
      generators = "E=ABCD"
    )

    expect_true(all(unlist(design) %in% c(-1, 1)))
  })

  it("has no missing values", {
    design <- fractional_factorial_design(
      num_factors = 6,
      generators = c("E=ABC", "F=BCD")
    )

    expect_false(anyNA(design))
  })

  it("produces unique rows (no duplicates)", {
    design <- fractional_factorial_design(
      num_factors = 4,
      generators = "D=ABC"
    )

    expect_equal(nrow(design), nrow(unique(design)))
  })
})
