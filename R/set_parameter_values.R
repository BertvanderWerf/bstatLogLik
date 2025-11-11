#' Set Parameter Values in parameter_info Object
#'
#' Updates parameter values in a \code{parameter_info} object using a named vector
#' of new values. Handles hierarchical parameter structure where parameters may
#' have sub-components (e.g., "k.A", "k.B" for parameter k with formula ~0 + group).
#'
#' @details
#' This function provides a safe way to update parameter values in a
#' \code{parameter_info} object. It:
#' \itemize{
#'   \item Validates that \code{param_info} is correct class
#'   \item Checks that all parameter names in \code{newvals} exist in \code{param_info}
#'   \item Parses hierarchical parameter names using separator (default: \code{param_info$sep})
#'   \item Updates only specified parameters, leaving others unchanged
#'   \item Returns modified \code{parameter_info} object
#' }
#'
#' **Parameter Naming Convention:**
#' Parameters in \code{parameter_info} can have sub-components when formulas
#' create design matrices. For example:
#' \itemize{
#'   \item Base parameter: "k"
#'   \item With formula ~0 + group: creates "k.A", "k.B" (for groups A, B)
#'   \item Separator (default "."): "k.A" = "k" + "." + "A"
#' }
#'
#' The function splits each name in \code{newvals} using the separator to find
#' the base parameter and sub-component.
#'
#' @param param_info Object of class \code{parameter_info} from
#'   \code{\link{create_parameter_info}}.
#' @param newvals Named numeric vector. Names must match expanded parameter names
#'   from \code{param_info} (e.g., "k.A", "lambda.B"). If NULL, returns
#'   \code{param_info} unchanged.
#' @param sep Character. Separator used to split parameter names into base parameter
#'   and sub-component. Default is \code{param_info$sep}. Typically ".".
#'
#' @return Modified \code{parameter_info} object with updated parameter values.
#'
#' @note
#' **Error Handling:**
#' The function will stop with an error if:
#' \itemize{
#'   \item \code{param_info} is not of class \code{parameter_info}
#'   \item \code{newvals} is not a numeric vector (or NULL)
#'   \item \code{newvals} has unnamed elements
#'   \item Any name in \code{newvals} does not exist in \code{param_info}
#' }
#'
#' **Use Cases:**
#' \itemize{
#'   \item Update parameters after optimization
#'   \item Set parameters from design of experiments
#'   \item Initialize parameters from previous fit
#'   \item Batch update multiple parameters
#' }
#'
#' @examples
#' \dontrun{
#' # Example 1: Simple parameter update
#' param_info <- create_parameter_info(
#'   parameter_values = c(k = 2, lambda = 10),
#'   estimate_flag = c(k = 1, lambda = 1)
#' )
#'
#' # Update both parameters
#' param_info <- set_parameter_values(
#'   param_info,
#'   newvals = c(k = 3.5, lambda = 12.5)
#' )
#'
#' # Example 2: With formula expansion
#' data <- data.frame(group = factor(c("A", "B", "A", "B")))
#' param_info <- create_parameter_info(
#'   parameter_values = c(k = 2, lambda = 10),
#'   formulas = list(k = ~0 + group, lambda = ~1),
#'   model_data = data
#' )
#'
#' # Now have kA, kB, lambda
#' # Update only kA and kB
#' param_info <- set_parameter_values(
#'   param_info,
#'   newvals = c(k.A = 2.5, k.B = 3.0)
#' )
#'
#' # Example 3: Use with design of experiments
#' design <- create_parameter_design(param_info)
#' best_point <- design[which.max(design$loglik), ]
#'
#' # Extract parameter values from best design point
#' param_names <- names(param_info$parameter_values)
#' new_values <- as.numeric(best_point[param_names])
#' names(new_values) <- param_names
#'
#' # Update param_info
#' param_info <- set_parameter_values(param_info, new_values)
#'
#' # Example 4: Partial update
#' # Only update k, leave lambda unchanged
#' param_info <- set_parameter_values(
#'   param_info,
#'   newvals = c(k = 5.0)
#' )
#' }
#'
#' @seealso
#' \code{\link{create_parameter_info}} for creating parameter_info objects
#'
#' @export
set_parameter_values <- function(param_info, newvals, sep = param_info$sep) {

  # ============================================================================
  # INPUT VALIDATION
  # ============================================================================

  # Validate param_info class
  if (!inherits(param_info, 'parameter_info')) {
    stop(
      "The parameter 'param_info' must be of class 'parameter_info'.",
      call. = FALSE
    )
  }

  # Handle NULL case - return unchanged
  if (is.null(newvals)) {
    return(param_info)
  }

  # Validate newvals is numeric
  bstatErr::check_numeric_vector(newvals, allow_null = TRUE)

  # Validate that newvals has names
  if (is.null(names(newvals)) || any(names(newvals) == '')) {
    stop(
      "'newvals' must be a named vector with names corresponding to ",
      "full parameter names (e.g., 'k.A', 'lambda.B')",
      call. = FALSE
    )
  }

  # ============================================================================
  # UPDATE PARAMETER VALUES
  # ============================================================================

  # Loop through each new value
  for (i in seq_along(newvals)) {
    full_name <- names(newvals)[i]

    # Split full name into base parameter and sub-component
    # e.g., "k.A" -> parname = "k", subname = "A"
    name_parts <- strsplit(full_name, sep, fixed = TRUE)[[1]]

    # First part is base parameter name
    parname <- name_parts[1]

    # Remaining part is sub-component name
    # Use sub() to remove base parameter and separator from full name
    subname <- sub(paste0('^', parname, sep), '', full_name)

    # Flag to track whether parameter was successfully set
    set_successfully <- FALSE

    # Check if base parameter exists in param_info
    if (parname %in% names(param_info$parameter_values)) {

      # Check if sub-component exists within base parameter
      if (subname %in% names(param_info$parameter_values[[parname]])) {

        # Update the value
        param_info$parameter_values[[parname]][[subname]] <- newvals[i]
        set_successfully <- TRUE

      }
    }

    # Error if parameter could not be set
    if (!set_successfully) {
      stop(
        sprintf(
          "Parameter 'param_info$parameter_values$%s$%s' does not exist. ",
          parname, subname
        ),
        "Available parameters: ",
        paste(names(unlist(param_info$parameter_values)), collapse = ", "),
        call. = FALSE
      )
    }
  }

  # ============================================================================
  # RETURN UPDATED PARAMETER_INFO
  # ============================================================================

  return(param_info)
}
