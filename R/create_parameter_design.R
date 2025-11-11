#' Create Design of Experiments for Parameter Space Exploration
#'
#' Generates a Central Composite Design (CCD) in parameter space based on
#' \code{parameter_info} bounds. The design can be used to explore the
#' log-likelihood surface and find good initial parameter estimates through
#' systematic evaluation of the parameter space.
#'
#' @details
#' This function creates a rotatable or face-centered CCD for exploring the
#' parameter space defined by a \code{parameter_info} object. The design:
#' \itemize{
#'   \item Uses lower and upper parameter bounds to define the design region
#'   \item Creates a fractional factorial design in coded units (-1, +1)
#'   \item Adds axial (star) points at ±α distance from center
#'   \item Adds center point(s) for curvature estimation
#'   \item Transforms coded values back to original parameter scale
#' }
#'
#' **Use Cases:**
#' \itemize{
#'   \item **Initial value search**: Evaluate log-likelihood at design points to
#'     find good starting values for optimization
#'   \item **Likelihood surface visualization**: Map the response surface of
#'     log-likelihood
#'   \item **Parameter sensitivity**: Identify which parameters most affect
#'     log-likelihood
#'   \item **Multi-start optimization**: Use design points as multiple starting
#'     values for optimization
#' }
#'
#' **Coding and Scaling:**
#' The design is created in coded units (-1 to +1), then transformed to the
#' original parameter scale using:
#' \deqn{x_{real} = \frac{x_{coded} \times (upper - lower)}{2} + \frac{upper + lower}{2}}
#'
#' **Handling Infinite Bounds:**
#' Parameters with infinite bounds (±Inf) are handled by:
#' \itemize{
#'   \item Using current parameter value as center
#'   \item Creating symmetric range based on \code{infinite_range_factor}
#'   \item Default factor of 2 means range = \code{[value-delta, value+delta]} for positive values
#' }
#' delta logic
#' \itemize{
#'  \item Both ±Inf:    delta = \code{|value|/infinite_range_factor}
#'  \item Lower -Inf:   delta = \code{min(upper-value, |value|/infinite_range_factor)}
#'  \item Upper +Inf:   delta = \code{min(value-lower, |value|/infinite_range_factor)}
#' }
#'
#' @param parameter_info Object of class \code{parameter_info} from
#'   \code{create_parameter_info()}. Must have reasonable (non-infinite) bounds
#'   or current parameter values will be used to define ranges.
#' @param generators Character vector or NULL. Generator definitions for fractional
#'   factorial design in terms of expanded parameter names (e.g., "k.B=k.A*lambda.A").
#'   If NULL, creates full factorial. See \code{\link{fractional_factorial_design}}
#'   for generator syntax.
#' @param alpha Character or numeric. Distance of axial points from center:
#'   \itemize{
#'     \item "rotatable" (default) - ensures rotatability
#'     \item "faces" - face-centered design (3 levels only)
#'     \item Numeric - custom alpha value
#'   }
#' @param center_points Integer. Number of center point replicates. Default is 1
#'   (typically sufficient for initial value search). Use more (4-6) for response
#'   surface fitting.
#' @param infinite_range_factor Numeric. For parameters with infinite bounds,
#'   defines range as \code{[value/factor, value*factor]} for positive values, or
#'   \code{[value*factor, value/factor]} for negative values. Default 2.
#'
#' @return Data frame with the parameter space design containing:
#'   \itemize{
#'     \item One column per expanded parameter (e.g., "kA", "kB", "lambdaA", "lambdaB")
#'     \item Values in original parameter scale (not coded units)
#'     \item "point_type" column: "factorial", "axial", or "center"
#'     \item Attributes: design_type, num_parameters, parameter_names, bounds_used
#'   }
#'
#' @note
#' **Parameter expansion**: If parameter_info contains formulas creating multiple
#' coefficients (e.g., k with formula ~0+group creates kA and kB), the design
#' includes all expanded parameters.
#'
#' **Bounds validation**: The function will warn if:
#' \itemize{
#'   \item Any bound is infinite and must be estimated
#'   \item Bounds are very wide (range > 1000)
#'   \item Bounds are very narrow (range < 0.01)
#' }
#'
#' **Design size**: Total runs = \eqn{n_f + 2k + n_c} where:
#' \itemize{
#'   \item \eqn{n_f} = factorial points (2^k or 2^(k-p) if generators used)
#'   \item \eqn{k} = number of parameters
#'   \item \eqn{n_c} = center points
#' }
#'
#' @references
#' Box, G.E.P. and Draper, N.R. (2007) Response Surfaces, Mixtures, and Ridge
#' Analyses, 2nd ed. Wiley.
#'
#' Jones, B. and Nachtsheim, C.J. (2011) Efficient Designs With Minimal
#' Aliasing. Technometrics, 53(1), 62-71.
#'
#' @examples
#' \dontrun{
#' # Example 1: Simple 2-parameter design
#' param_info <- create_parameter_info(
#'   parameter_values = c(k = 2, lambda = 10),
#'   lower_bounds = c(k = 0.5, lambda = 5),
#'   upper_bounds = c(k = 5, lambda = 20),
#'   estimate_flag = c(k = 1, lambda = 1)
#' )
#'
#' design <- create_parameter_design(param_info)
#' print(design)
#' # Returns CCD with k and lambda in original scale
#'
#' # Example 2: With parameter expansion (formulas)
#' data <- data.frame(group = factor(c("A", "B", "A", "B")))
#' param_info <- create_parameter_info(
#'   parameter_values = c(k = 2, lambda = 10),
#'   formulas = list(k = ~0 + group, lambda = ~0 + group),
#'   lower_bounds = c(k = 0.5, lambda = 5),
#'   upper_bounds = c(k = 5, lambda = 20),
#'   model_data = data
#' )
#'
#' # Now have 4 parameters: kA, kB, lambdaA, lambdaB
#' design <- create_parameter_design(
#'   param_info,
#'   generators = "lambdaB=kA*kB*lambdaA"  # 2^(4-1) design
#' )
#'
#' # Example 3: Face-centered for limited levels
#' design_fccd <- create_parameter_design(
#'   param_info,
#'   alpha = "faces",
#'   center_points = 3
#' )
#'
#' # Example 4: Use design for initial value search
#' # Evaluate log-likelihood at each design point
#' design$loglik <- apply(design[, 1:4], 1, function(pars) {
#'   # Update parameter_info with design values
#'   param_info$parameter_values <- as.list(pars)
#'   # Compute log-likelihood
#'   ll <- compute_loglik(param_info, data)
#'   return(ll)
#' })
#'
#' # Find best starting point
#' best_idx <- which.max(design$loglik)
#' initial_values <- design[best_idx, 1:4]
#' }
#'
#' @seealso
#' \code{\link{fractional_factorial_design}} for creating the base factorial design
#'
#' \code{\link{augment_with_star_center_points}} for adding star and center points
#'
#' \code{\link{create_parameter_info}} for parameter structure
#'
#' @export
create_parameter_design <- function(parameter_info,
                                     generators = NULL,
                                     alpha = "rotatable",
                                     center_points = 1,
                                     infinite_range_factor = 2) {

  # ============================================================================
  # INPUT VALIDATION
  # ============================================================================

  if (!inherits(parameter_info, "parameter_info")) {
    stop("'parameter_info' must be an object of class 'parameter_info'",
         call. = FALSE)
  }

  # ============================================================================
  # EXTRACT PARAMETER INFORMATION
  # ============================================================================

  # Get expanded parameter values (accounts for design matrices)
  param_values <- unlist(parameter_info$parameter_values)
  param_names <- names(param_values)
  num_params <- length(param_values)

  # Get bounds
  lower_bounds <- unlist(parameter_info$lower_bounds)
  upper_bounds <- unlist(parameter_info$upper_bounds)

  # Only include parameters that should be estimated
  estimate_flags <- unlist(parameter_info$estimate_flag)
  params_to_design <- estimate_flags == 1

  if (sum(params_to_design) == 0) {
    stop(
      "No parameters marked for estimation (all estimate_flag = 0). ",
      "Cannot create design.",
      call. = FALSE
    )
  }

  # Subset to parameters being estimated
  param_values_design <- param_values[params_to_design]
  param_names_design <- param_names[params_to_design]
  lower_bounds_design <- lower_bounds[params_to_design]
  upper_bounds_design <- upper_bounds[params_to_design]
  num_params_design <- length(param_values_design)

  message("Creating design for ", num_params_design, " parameters:")
  message("  ", paste(param_names_design, collapse = ", "))

  # ============================================================================
  # HANDLE INFINITE BOUNDS (CENTER AT CURRENT VALUE)
  # ============================================================================
  # For parameters with infinite bounds, create reasonable ranges
  # centered at the current parameter value (not [value/2, value*2])
  #
  # Logic:
  # - Both infinite: delta = value/2, range = [value-delta, value+delta]
  # - Lower infinite: delta = min(upperbound-value, value/2)
  # - Upper infinite: delta = min(value-lowerbound, value/2)

  has_inf_lower <- is.infinite(lower_bounds_design)
  has_inf_upper <- is.infinite(upper_bounds_design)

  if (any(has_inf_lower) || any(has_inf_upper)) {
    warning(
      "Some parameters have infinite bounds. ",
      "Creating symmetric ranges centered at current values.",
      call. = FALSE
    )

    for (i in seq_along(param_values_design)) {
      current_val <- param_values_design[i]
      lower_val <- lower_bounds_design[i]
      upper_val <- upper_bounds_design[i]

      if (has_inf_lower[i] && has_inf_upper[i]) {
        # Both bounds infinite
        # Use delta = value/2 to create symmetric range
        delta <- abs(current_val) / 2

        # Handle case where current_val is 0
        if (delta < 1e-6) {
          delta <- infinite_range_factor  # Use factor as default delta
        }

        lower_bounds_design[i] <- current_val - delta
        upper_bounds_design[i] <- current_val + delta

      } else if (has_inf_lower[i]) {
        # Only lower bound infinite, upper bound is finite
        # delta = min(upperbound - value, value/2)
        delta_from_upper <- upper_val - current_val
        delta_from_value <- abs(current_val) / 2

        delta <- min(delta_from_upper, delta_from_value)

        # Ensure positive delta
        if (delta <= 0) {
          delta <- (upper_val - current_val) * 0.5
        }

        lower_bounds_design[i] <- current_val - delta

      } else if (has_inf_upper[i]) {
        # Only upper bound infinite, lower bound is finite
        # delta = min(value - lowerbound, value/2)
        delta_from_lower <- current_val - lower_val
        delta_from_value <- abs(current_val) / 2

        delta <- min(delta_from_lower, delta_from_value)

        # Ensure positive delta
        if (delta <= 0) {
          delta <- (current_val - lower_val) * 0.5
        }

        upper_bounds_design[i] <- current_val + delta
      }

      message(
        "  ", param_names_design[i], ": [",
        round(lower_bounds_design[i], 4), ", ",
        round(upper_bounds_design[i], 4), "] centered at ",
        round(current_val, 4)
      )
    }
  }

  # ============================================================================
  # VALIDATE BOUNDS
  # ============================================================================

  # Check that lower < upper
  if (any(lower_bounds_design >= upper_bounds_design)) {
    bad_params <- param_names_design[lower_bounds_design >= upper_bounds_design]
    stop(
      "Lower bounds must be less than upper bounds for: ",
      paste(bad_params, collapse = ", "),
      call. = FALSE
    )
  }

  # Warn about very wide or narrow ranges
  param_ranges <- upper_bounds_design - lower_bounds_design

  if (any(param_ranges > 1000)) {
    wide_params <- param_names_design[param_ranges > 1000]
    warning(
      "Very wide parameter ranges for: ",
      paste(wide_params, collapse = ", "),
      ". Consider narrowing bounds for better design.",
      call. = FALSE
    )
  }

  if (any(param_ranges < 0.01)) {
    narrow_params <- param_names_design[param_ranges < 0.01]
    warning(
      "Very narrow parameter ranges for: ",
      paste(narrow_params, collapse = ", "),
      ". Consider fixing these parameters.",
      call. = FALSE
    )
  }

  # ============================================================================
  # CREATE CODED DESIGN (in -1 to +1 scale)
  # ============================================================================

  # Create fractional factorial design in coded units
  factorial_coded <- fractional_factorial_design(
    num_factors = num_params_design,
    generators = generators
  )

  # Augment with star and center points
  ccd_coded <- augment_with_star_center_points(
    factorial_coded,
    alpha = alpha,
    center_points = center_points,
    add_blocks = FALSE  # Don't need blocks for parameter exploration
  )

  # ============================================================================
  # TRANSFORM FROM CODED (-1,+1) TO REAL PARAMETER SCALE
  # ============================================================================
  # Transformation: x_real = (x_coded * range/2) + center
  # where center = (upper + lower)/2
  #       range = upper - lower

  param_centers <- (upper_bounds_design + lower_bounds_design) / 2
  param_half_ranges <- (upper_bounds_design - lower_bounds_design) / 2

  # Extract just the parameter columns (exclude point_type)
  param_cols <- names(ccd_coded)[names(ccd_coded) != "point_type"]

  design_real <- ccd_coded
  for (i in seq_along(param_cols)) {
    # Transform coded values to real scale
    design_real[[param_cols[i]]] <- (
      ccd_coded[[param_cols[i]]] * param_half_ranges[i] + param_centers[i]
    )
  }

  # Rename columns to actual parameter names
  names(design_real)[match(param_cols, names(design_real))] <- param_names_design

  # ============================================================================
  # ADD FIXED PARAMETERS (if any)
  # ============================================================================
  # Parameters not being estimated are set to their current values

  if (any(!params_to_design)) {
    fixed_params <- param_names[!params_to_design]
    fixed_values <- param_values[!params_to_design]

    for (i in seq_along(fixed_params)) {
      design_real[[fixed_params[i]]] <- fixed_values[i]
    }

    message("\nFixed parameters (not in design):")
    message("  ", paste(paste0(fixed_params, " = ", round(fixed_values, 4)),
                        collapse = ", "))
  }

  # ============================================================================
  # REORDER COLUMNS
  # ============================================================================
  # Put all parameters first (in original order), then point_type

  design_real <- design_real[, c(param_names, "point_type")]

  # ============================================================================
  # ADD DESIGN ATTRIBUTES
  # ============================================================================

  attr(design_real, "design_type") <- "parameter_space_ccd"
  attr(design_real, "num_parameters") <- num_params_design
  attr(design_real, "parameter_names") <- param_names_design
  attr(design_real, "fixed_parameters") <- param_names[!params_to_design]
  attr(design_real, "bounds_used") <- data.frame(
    parameter = param_names_design,
    lower = lower_bounds_design,
    upper = upper_bounds_design,
    center = param_centers,
    range = param_ranges
  )
  attr(design_real, "alpha") <- attr(ccd_coded, "alpha")
  attr(design_real, "rotatable") <- attr(ccd_coded, "rotatable")
  attr(design_real, 'coded_design') <- ccd_coded

  # ============================================================================
  # DISPLAY DESIGN SUMMARY
  # ============================================================================

  message("\n=== Parameter Space Design Summary ===")
  message("Design points in parameter space")
  message("Total runs: ", nrow(design_real))
  message("  Factorial: ", sum(design_real$point_type == "factorial"))
  message("  Axial:     ", sum(design_real$point_type == "axial"))
  message("  Center:    ", sum(design_real$point_type == "center"))
  message("\nParameter ranges:")
  for (i in seq_along(param_names_design)) {
    message(sprintf("  %-15s [%8.4f, %8.4f]",
                    param_names_design[i],
                    lower_bounds_design[i],
                    upper_bounds_design[i]))
  }

  if (attr(design_real, "rotatable")) {
    message("\nDesign is ROTATABLE (alpha = ", round(attr(design_real, "alpha"), 4), ")")
  }

  # ============================================================================
  # RETURN DESIGN
  # ============================================================================

  class(design_real) <- c("parameter_design", "data.frame")
  return(design_real)
}


#' Print method for parameter_design objects
#'
#' @param x Object of class parameter_design
#' @param ... Additional arguments
#'
#' @export
print.parameter_design <- function(x, ...) {
  cat("\nParameter Space Design (CCD)\n")
  cat("============================\n\n")

  # Extract attributes
  num_params <- attr(x, "num_parameters")
  param_names <- attr(x, "parameter_names")
  bounds <- attr(x, "bounds_used")

  cat("Parameters in design:", num_params, "\n")
  cat("  ", paste(param_names, collapse = ", "), "\n\n")

  cat("Design points:\n")
  print(table(x$point_type))

  cat("\nParameter ranges:\n")
  print(bounds, row.names = FALSE)

  if (attr(x, "rotatable")) {
    cat("\nRotatable design (alpha =", round(attr(x, "alpha"), 4), ")\n")
  }

  cat("\nFirst few design points:\n")
  print(head(as.data.frame(x), 10))

  invisible(x)
}
