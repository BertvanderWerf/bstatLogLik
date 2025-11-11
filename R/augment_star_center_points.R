#' Augment Design with Star and Center Points for Response Surface Methodology
#'
#' Adds axial (star) points and center points to a two-level factorial or fractional
#' factorial design to create a Central Composite Design (CCD) suitable for fitting
#' second-order (quadratic) response surface models.
#'
#' @details
#' **Central Composite Designs (CCD)** are the most popular response surface methodology
#' (RSM) designs for fitting quadratic models. A CCD consists of three components:
#' \itemize{
#'   \item **Factorial points**: 2^k or 2^(k-p) design from \code{fractional_factorial_design()}
#'   \item **Axial (star) points**: 2k points at distance ±α on each factor axis
#'   \item **Center points**: Replicated points at the origin (0, 0, ..., 0)
#' }
#'
#' **Three Types of CCD:**
#' \itemize{
#'   \item **CCC (Circumscribed)**: α > 1, star points outside factorial cube (default)
#'   \item **CCF (Face-centered)**: α = 1, star points on face of factorial cube
#'   \item **CCI (Inscribed)**: α < 1, star points inside factorial cube
#' }
#'
#' **Rotatability**: A design is rotatable when the variance of predicted response
#' depends only on distance from center, not direction. For rotatability:
#' \deqn{\alpha = (n_f)^{1/4}}
#' where \eqn{n_f} is the number of factorial points.
#'
#' **Orthogonality**: When blocks are used, orthogonal designs allow independent
#' estimation of model terms and block effects.
#'
#' @param factorial_design Data frame. Output from \code{fractional_factorial_design()}
#'   containing the factorial points coded as -1 and +1.
#' @param alpha Numeric. Distance of axial points from center. Default is "rotatable"
#'   which computes \eqn{\alpha = (n_f)^{1/4}}. Can also be:
#'   \itemize{
#'     \item Numeric value (e.g., 1.414 for 2 factors)
#'     \item "rotatable" (default) - ensures rotatability
#'     \item "orthogonal" - for orthogonal blocking
#'     \item "faces" or 1 - face-centered design (α=1)
#'   }
#' @param center_points Integer. Number of center point replicates. Default is
#'   calculated to balance design, typically 4-6 points. More replicates improve
#'   pure error estimation but increase experimental cost.
#' @param add_blocks Logical. If TRUE, adds a "block" column to identify factorial,
#'   axial, and center point blocks. Default FALSE.
#'
#' @return Data frame with the augmented design containing:
#'   \itemize{
#'     \item All original factorial points (coded -1, +1)
#'     \item 2k axial points (coded 0 except one factor at ±α)
#'     \item n center points (all factors coded 0)
#'     \item Optional "block" column if \code{add_blocks = TRUE}
#'     \item "point_type" column identifying: "factorial", "axial", or "center"
#'   }
#'
#' @note
#' **Design sizes:**
#' \itemize{
#'   \item Total runs = \eqn{n_f + 2k + n_c}
#'   \item \eqn{n_f} = factorial points (8 for 2^3, 16 for 2^4, etc.)
#'   \item \eqn{2k} = axial points (always twice number of factors)
#'   \item \eqn{n_c} = center points (user specified)
#' }
#'
#' **Alpha selection guide:**
#' \itemize{
#'   \item For k=2: α = 1.414 (rotatable)
#'   \item For k=3: α = 1.682 (rotatable)
#'   \item For k=4: α = 2.000 (rotatable)
#'   \item For k=5: α = 2.378 (rotatable)
#' }
#'
#' **Rotatability vs Face-centered:**
#' \itemize{
#'   \item Rotatable designs (α > 1) explore larger region, need 5 factor levels
#'   \item Face-centered (α = 1) need only 3 levels, easier to implement
#' }
#'
#' @references
#' Box, G.E.P. and Wilson, K.B. (1951) On the Experimental Attainment of
#' Optimum Conditions. Journal of the Royal Statistical Society B, 13, 1-45.
#'
#' Montgomery, D.C. (2017) Design and Analysis of Experiments, 9th ed. Wiley.
#'
#' Myers, R.H., Montgomery, D.C., and Anderson-Cook, C.M. (2016) Response
#' Surface Methodology: Process and Product Optimization Using Designed
#' Experiments, 4th ed. Wiley.
#'
#' @examples
#' # Example 1: 2^3 full factorial + CCD (rotatable)
#' factorial <- fractional_factorial_design(num_factors = 3)
#' ccd <- augment_with_star_center_points(factorial)
#' print(ccd)
#' # Total: 8 factorial + 6 axial + 6 center = 20 runs
#'
#' # Example 2: Face-centered design (alpha = 1)
#' factorial <- fractional_factorial_design(num_factors = 2)
#' ccf <- augment_with_star_center_points(
#'   factorial,
#'   alpha = "faces",
#'   center_points = 5
#' )
#' print(ccf)
#'
#' # Example 3: 2^(4-1) fractional + CCD with blocks
#' factorial <- fractional_factorial_design(
#'   num_factors = 4,
#'   generators = "D=ABC"
#' )
#' ccd_blocked <- augment_with_star_center_points(
#'   factorial,
#'   alpha = "rotatable",
#'   add_blocks = TRUE
#' )
#' table(ccd_blocked$block, ccd_blocked$point_type)
#'
#' # Example 4: Custom alpha value
#' factorial <- fractional_factorial_design(num_factors = 3)
#' ccd_custom <- augment_with_star_center_points(
#'   factorial,
#'   alpha = 1.5,
#'   center_points = 4
#' )
#'
#' @seealso
#' \code{\link{fractional_factorial_design}} for creating the base factorial design
#'
#' @export
augment_with_star_center_points <- function(factorial_design,
                                             alpha = "rotatable",
                                             center_points = NULL,
                                             add_blocks = FALSE) {

  # ============================================================================
  # INPUT VALIDATION
  # ============================================================================

  # Validate factorial_design
  if (!is.data.frame(factorial_design)) {
    stop("'factorial_design' must be a data frame", call. = FALSE)
  }

  if (nrow(factorial_design) == 0) {
    stop("'factorial_design' is empty", call. = FALSE)
  }

  # Ensure all columns are numeric
  if (!all(sapply(factorial_design, is.numeric))) {
    stop(
      "'factorial_design' must contain only numeric columns ",
      "(coded factor levels)",
      call. = FALSE
    )
  }

  # ============================================================================
  # EXTRACT DESIGN PARAMETERS
  # ============================================================================

  num_factors <- ncol(factorial_design)
  num_factorial_points <- nrow(factorial_design)
  factor_names <- names(factorial_design)

  # ============================================================================
  # DETERMINE ALPHA VALUE
  # ============================================================================

  if (is.character(alpha)) {
    alpha_input <- tolower(alpha)

    if (alpha_input == "rotatable") {
      # For rotatability: alpha = (n_factorial)^(1/4)
      alpha_value <- num_factorial_points^(1/4)
      message(
        "Using rotatable alpha = ", round(alpha_value, 4),
        " (computed as n_factorial^0.25)"
      )

    } else if (alpha_input %in% c("faces", "face", "face-centered")) {
      # Face-centered: alpha = 1
      alpha_value <- 1
      message("Using face-centered design with alpha = 1")

    } else if (alpha_input == "orthogonal") {
      # Orthogonal blocking - requires more complex calculation
      # Simplified: use alpha that maintains orthogonality
      # Formula depends on number of center points
      # For now, use rotatable as approximation
      alpha_value <- num_factorial_points^(1/4)
      warning(
        "Orthogonal blocking alpha calculation is complex. ",
        "Using rotatable alpha = ", round(alpha_value, 4),
        " as approximation. ",
        "For precise orthogonality, specify numeric alpha.",
        call. = FALSE
      )

    } else {
      stop(
        "'alpha' must be 'rotatable', 'orthogonal', 'faces', or a numeric value",
        call. = FALSE
      )
    }

  } else if (is.numeric(alpha)) {
    # User-specified alpha
    if (alpha <= 0) {
      stop("'alpha' must be positive", call. = FALSE)
    }
    alpha_value <- alpha

  } else {
    stop(
      "'alpha' must be 'rotatable', 'orthogonal', 'faces', or a numeric value",
      call. = FALSE
    )
  }

  # ============================================================================
  # DETERMINE NUMBER OF CENTER POINTS
  # ============================================================================

  if (is.null(center_points)) {
    # Default center points based on design size
    # Standard practice: 4-6 center points for good pure error estimation
    # More factors -> more center points
    center_points <- max(4, min(6, ceiling(num_factors * 1.5)))

    message(
      "Using ", center_points, " center points ",
      "(default based on ", num_factors, " factors)"
    )

  } else {
    # Validate user-specified center points
    if (!is.numeric(center_points) ||
        center_points < 1 ||
        center_points != round(center_points)) {
      stop("'center_points' must be a positive integer", call. = FALSE)
    }
  }

  # ============================================================================
  # CREATE AXIAL (STAR) POINTS
  # ============================================================================
  # Axial points: 2k points at (±alpha, 0, 0, ...) for each factor
  # Total: 2 * num_factors points

  num_axial_points <- 2 * num_factors

  # Initialize matrix for axial points
  axial_points <- matrix(
    0,
    nrow = num_axial_points,
    ncol = num_factors
  )

  # Fill in axial points
  # For factor i:
  #   - Row 2*i-1: alpha in position i, 0 elsewhere
  #   - Row 2*i:   -alpha in position i, 0 elsewhere
  for (i in seq_len(num_factors)) {
    axial_points[2*i - 1, i] <- alpha_value   # +alpha
    axial_points[2*i, i] <- -alpha_value      # -alpha
  }

  # Convert to data frame with proper column names
  axial_df <- as.data.frame(axial_points)
  names(axial_df) <- factor_names

  # ============================================================================
  # CREATE CENTER POINTS
  # ============================================================================
  # Center points: all factors set to 0 (center of design space)
  # Replicated for pure error estimation

  center_df <- as.data.frame(
    matrix(
      0,
      nrow = center_points,
      ncol = num_factors
    )
  )
  names(center_df) <- factor_names

  # ============================================================================
  # COMBINE ALL DESIGN POINTS
  # ============================================================================

  # Add point type identifier
  factorial_design$point_type <- "factorial"
  axial_df$point_type <- "axial"
  center_df$point_type <- "center"

  # Combine all points
  augmented_design <- rbind(
    factorial_design,
    axial_df,
    center_df
  )

  # Convert point_type to factor with ordered levels
  augmented_design$point_type <- factor(
    augmented_design$point_type,
    levels = c("factorial", "axial", "center")
  )

  # ============================================================================
  # ADD BLOCK INFORMATION (if requested)
  # ============================================================================

  if (add_blocks) {
    # Standard blocking scheme for CCD:
    # - Block 1: Factorial points
    # - Block 2: Axial points
    # - Block 3: Center points (can be split across blocks)
    #
    # Note: For orthogonal blocking, center points may be distributed
    # across blocks. Here we use simple scheme.

    augmented_design$block <- c(
      rep(1, num_factorial_points),          # Factorial block
      rep(2, num_axial_points),               # Axial block
      rep(3, center_points)                   # Center block
    )

    augmented_design$block <- factor(augmented_design$block)
  }

  # ============================================================================
  # ADD DESIGN ATTRIBUTES
  # ============================================================================

  attr(augmented_design, "design_type") <- "central_composite"
  attr(augmented_design, "alpha") <- alpha_value
  attr(augmented_design, "num_factors") <- num_factors
  attr(augmented_design, "num_factorial") <- num_factorial_points
  attr(augmented_design, "num_axial") <- num_axial_points
  attr(augmented_design, "num_center") <- center_points
  attr(augmented_design, "rotatable") <- (abs(alpha_value - num_factorial_points^0.25) < 1e-6)

  # ============================================================================
  # DISPLAY DESIGN SUMMARY
  # ============================================================================

  message("\n=== Central Composite Design Summary ===")
  message("Factors: ", num_factors)
  message("Factorial points: ", num_factorial_points)
  message("Axial points: ", num_axial_points)
  message("Center points: ", center_points)
  message("Total runs: ", nrow(augmented_design))
  message("Alpha: ", round(alpha_value, 4))

  if (attr(augmented_design, "rotatable")) {
    message("Design is ROTATABLE")
  } else {
    message("Design is NOT rotatable (rotatable alpha = ",
            round(num_factorial_points^0.25, 4), ")")
  }

  # ============================================================================
  # RETURN AUGMENTED DESIGN
  # ============================================================================

  return(augmented_design)
}


#' Print method for central composite designs
#'
#' @param x Central composite design from augment_with_star_center_points()
#' @param ... Additional arguments
#'
#' @export
print.central_composite <- function(x, ...) {
  cat("\nCentral Composite Design\n")
  cat("========================\n\n")

  # Extract attributes
  num_factors <- attr(x, "num_factors")
  alpha <- attr(x, "alpha")
  rotatable <- attr(x, "rotatable")

  cat("Number of factors:", num_factors, "\n")
  cat("Alpha value:", round(alpha, 4), "\n")
  cat("Rotatable:", ifelse(rotatable, "Yes", "No"), "\n")
  cat("\nDesign points:\n")
  cat("  Factorial:", attr(x, "num_factorial"), "\n")
  cat("  Axial:    ", attr(x, "num_axial"), "\n")
  cat("  Center:   ", attr(x, "num_center"), "\n")
  cat("  Total:    ", nrow(x), "\n\n")

  # Show point type distribution
  cat("Point type distribution:\n")
  print(table(x$point_type))

  cat("\nFirst few rows:\n")
  print(head(x, 10))

  invisible(x)
}
