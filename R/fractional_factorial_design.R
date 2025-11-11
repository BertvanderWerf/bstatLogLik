#' Create Fractional Factorial Design Matrix
#'
#' Generates a two-level fractional factorial design using specified generators
#' (defining relations) to reduce the number of experimental runs while maintaining
#' orthogonality and balance.
#'
#' @details
#' This function creates 2^(n-p) fractional factorial designs where:
#' \itemize{
#'   \item \code{n} = number of factors
#'   \item \code{p} = number of generators (defining relations)
#'   \item Design requires 2^(n-p) runs instead of 2^n runs
#' }
#'
#' The function supports:
#' \itemize{
#'   \item Full factorial designs (when \code{generators = NULL})
#'   \item Fractional factorial designs via generator specification
#'   \item Automatic confounding detection and warnings
#'   \item Standard orthogonal encoding (-1, +1)
#' }
#'
#' **Generators** (defining relations) specify which higher-order interactions
#' are confounded with main effects or lower-order interactions. For example:
#' \itemize{
#'   \item "D=ABC" means factor D is confounded with the ABC interaction
#'   \item This creates a 2^(4-1) = 8-run design from a 2^4 = 16-run design
#' }
#'
#' **Confounding** occurs when effects cannot be separated. The function
#' automatically detects and warns about confounded effects based on the
#' generators provided.
#'
#' **Resolution** describes aliasing structure:
#' \itemize{
#'   \item Resolution III: Main effects confounded with 2-factor interactions
#'   \item Resolution IV: Main effects free, some 2-factor interactions confounded
#'   \item Resolution V: Main effects and 2-factor interactions free
#' }
#'
#' @param num_factors Integer. Number of factors in the design. Must be > 0.
#' @param generators Character vector or NULL. Generator definitions for fractional
#'   designs. Each element specifies a defining relation like "D=ABC" or "E=-ABD".
#'   Use "-" prefix for negative generators (e.g., "D=-ABC"). If NULL, creates
#'   full factorial design.
#' @param levels Integer. Number of levels per factor. Currently only 2 is supported.
#'   For 3+ level designs, see NIST Engineering Statistics Handbook.
#'
#' @return Data frame with 2^(num_factors - num_generators) rows and num_factors
#'   columns. Factor levels are coded as -1 (low) and +1 (high). Column names
#'   are A, B, C, ... (uppercase then lowercase letters).
#'
#' @note
#' **Generator Notation:**
#' \itemize{
#'   \item "D=ABC" creates I = ABCD defining relation
#'   \item "D=-ABC" creates I = -ABCD defining relation
#'   \item Spaces and "=" are automatically removed for parsing
#' }
#'
#' **Confounding Structure:**
#' When generators create confounding, the function warns which effects cannot
#' be separated. For example, with generator "D=ABC":
#' \itemize{
#'   \item Main effect A confounded with BCD interaction
#'   \item Main effect B confounded with ACD interaction
#'   \item etc.
#' }
#'
#' @references
#' Montgomery, D.C. (2017) Design and Analysis of Experiments, 9th ed. Wiley.
#'
#' Box, G.E.P., Hunter, J.S., and Hunter, W.G. (2005) Statistics for
#' Experimenters, 2nd ed. Wiley.
#'
#' NIST/SEMATECH e-Handbook of Statistical Methods:
#' https://www.itl.nist.gov/div898/handbook/pri/section3/pri3347.htm
#'
#' @examples
#' # Full factorial design: 2^3 = 8 runs
#' design_full <- fractional_factorial_design(num_factors = 3)
#' print(design_full)
#'
#' # Fractional design: 2^(4-1) = 8 runs instead of 16
#' design_frac <- fractional_factorial_design(
#'   num_factors = 4,
#'   generators = "D=ABC"
#' )
#' print(design_frac)
#'
#' # Multiple generators: 2^(6-2) = 16 runs instead of 64
#' design_multi <- fractional_factorial_design(
#'   num_factors = 6,
#'   generators = c("E=ABC", "F=BCD")
#' )
#'
#' # Negative generator
#' design_neg <- fractional_factorial_design(
#'   num_factors = 5,
#'   generators = c("D=ABC", "E=-ABD")
#' )
#'
#' @export
fractional_factorial_design <- function(num_factors, generators = NULL, levels = 2) {

  # ============================================================================
  # INPUT VALIDATION
  # ============================================================================

  # Validate levels (currently only 2-level designs supported)
  if (levels != 2) {
    stop(
      'Currently only 2-level designs are supported. ',
      'For 3+ level designs, see NIST Engineering Statistics Handbook.',
      call. = FALSE
    )
  }

  # Validate num_factors
  if (!is.numeric(num_factors) || num_factors <= 0 || num_factors != round(num_factors)) {
    stop('num_factors must be a positive integer', call. = FALSE)
  }

  # ============================================================================
  # FULL FACTORIAL DESIGN (no generators)
  # ============================================================================

  if (is.null(generators)) {
    # Create full 2^n factorial design
    # Each factor has levels {-1, +1}
    factor_levels_list <- lapply(
      vector("list", num_factors),
      function(x) c(-1, 1)
    )

    # Generate all combinations via expand.grid
    design_matrix <- do.call(expand.grid, factor_levels_list)

    # Assign factor names: A, B, C, ..., Z, a, b, c, ...
    factor_names <- c(LETTERS, letters)[seq_len(num_factors)]
    names(design_matrix) <- factor_names

    return(design_matrix)
  }

  # ============================================================================
  # FRACTIONAL FACTORIAL DESIGN (with generators)
  # ============================================================================

  # Validate generators
  bstatErr::check_string_vector(generators)

  # Clean and parse generators
  # Remove spaces and "=" signs, convert to uppercase
  generators <- toupper(generators)
  generators_clean <- gsub('=', '', generators, fixed = TRUE)
  generators_clean <- gsub(' ', '', generators_clean, fixed = TRUE)

  # Split each generator into individual letters
  generators_split <- lapply(generators_clean, function(gen) {
    strsplit(gen, split = '')[[1]]
  })

  num_generators <- length(generators_split)

  # Check that num_generators < num_factors
  if (num_generators >= num_factors) {
    stop(
      'The number of generators (', num_generators, ') must be less than ',
      'num_factors (', num_factors, ')',
      call. = FALSE
    )
  }

  # ============================================================================
  # BUILD BASE DESIGN (factors not defined by generators)
  # ============================================================================

  # Number of base factors (not generated)
  num_base_factors <- num_factors - num_generators

  # Create full factorial for base factors
  base_levels_list <- lapply(
    vector("list", num_base_factors),
    function(x) c(-1, 1)
  )

  design_matrix <- do.call(expand.grid, base_levels_list)

  # Assign names to base factors (A, B, C, ...)
  base_factor_names <- LETTERS[seq_len(num_base_factors)]
  names(design_matrix) <- base_factor_names

  # ============================================================================
  # APPLY GENERATORS TO CREATE ADDITIONAL FACTORS
  # ============================================================================

  # Check for negative generators (indicated by "-" in original)
  is_negative_generator <- grepl("-", generators, fixed = TRUE)

  # Store original generator definitions for output
  generator_definitions <- generators_split

  for (i in seq_along(generators_split)) {
    gen_letters <- generators_split[[i]]
    gen_letters <- gen_letters[gen_letters != "-"]

    # Separate left-hand side (new factor) from right-hand side (existing factors)
    # Letters already in design_matrix are on RHS
    # Letters NOT in design_matrix are on LHS (new factor to create)
    rhs_letters <- gen_letters[gen_letters %in% names(design_matrix)]
    lhs_letters <- gen_letters[!(gen_letters %in% names(design_matrix))]

    # Validate generator structure
    if (length(lhs_letters) != 1 || length(rhs_letters) < 2) {
      stop(
        sprintf(
          "Generator %d ('%s') is not properly defined. ",
          i, generators[i]
        ),
        "The left-hand side must contain exactly one new factor letter, ",
        "and the right-hand side must contain 2+ existing factor letters from {",
        paste(names(design_matrix), collapse = ""), "}.",
        call. = FALSE
      )
    }

    # Apply generator: new factor = product of RHS factors
    # Multiply factor columns element-wise
    multiplier <- if (is_negative_generator[i]) -1 else 1
    design_matrix[[lhs_letters]] <- multiplier * apply(
      design_matrix[, rhs_letters, drop = FALSE],
      MARGIN = 1,
      FUN = prod
    )
  }

  # ============================================================================
  # DETECT CONFOUNDING STRUCTURE
  # ============================================================================

  # Compute generalized interactions (products of generators)
  # This reveals the complete confounding (aliasing) structure
  all_generators <- generator_definitions

  if (num_generators > 1) {
    # Generate all pairwise products of generators
    confounding_index <- num_generators + 1

    for (i in seq_len(num_generators - 1)) {
      for (j in (i + 1):num_generators) {
        # Find letters in generator i but not in generator j
        letters_only_i <- generator_definitions[[i]][
          !(generator_definitions[[i]] %in% generator_definitions[[j]])
        ]
        letters_only_i <- letters_only_i[letters_only_i != "-"]

        # Find letters in generator j but not in generator i
        letters_only_j <- generator_definitions[[j]][
          !(generator_definitions[[j]] %in% generator_definitions[[i]])
        ]
        letters_only_j <- letters_only_j[letters_only_j != "-"]

        # Product of generators gives confounded interaction
        confounded_letters <- sort(c(letters_only_i, letters_only_j))

        if (xor(is_negative_generator[i],is_negative_generator[j])) {
          confounded_letters[1] <- paste0("-", confounded_letters[1])
        }

        # Store for output
        all_generators[[confounding_index]] <- confounded_letters

        # Issue warnings about confounding
        if (length(confounded_letters) == 2) {
          warning(
            "Confounding detected: Main effect ", confounded_letters[1],
            " is confounded with main effect ", confounded_letters[2],
            call. = FALSE
          )
        } else if (length(confounded_letters) == 3) {
          warning(
            "Confounding detected: Main effect ", confounded_letters[1],
            " is confounded with interaction ", paste(confounded_letters[2:3], collapse = ":"),
            call. = FALSE
          )
        }
        # else if (length(confounded_letters) > 3) {
        #   warning(
        #     "Confounding detected: Interaction ",
        #     paste(confounded_letters[1:(length(confounded_letters)-1)], collapse = ":"),
        #     " is confounded with effect ", confounded_letters[length(confounded_letters)],
        #     call. = FALSE
        #   )
        # }

        confounding_index <- confounding_index + 1
      }
    }
  }

  # ============================================================================
  # DISPLAY DEFINING RELATION
  # ============================================================================

  # Print defining relation: I = generator1 = generator2 = ...
  defining_relation_terms <- sapply(all_generators,
                                    function(x) { paste(sort(x), collapse = '')})
  message(
    "Defining relation: I = ",
    paste(defining_relation_terms, collapse = " = ")
  )

  # ============================================================================
  # RETURN DESIGN MATRIX
  # ============================================================================

  return(design_matrix)
}
