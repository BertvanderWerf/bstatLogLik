#' Parameter Information Structure
#'
#' S3 class defining the structure to hold model parameter details for likelihood-based fitting.
#'
#' Objects of this class contain:
#'   \itemize{
#'     \item{\code{parameter_values}: Named list of initial parameter values used in modeling.}
#'     \item{\code{lower_bounds}: Lower bounds for parameter estimates.}
#'     \item{\code{upper_bounds}: Upper bounds for parameter estimates.}
#'   \item{\code{estimate_flag}: Logical or numeric vector (1/0) indicating whether each
#'     parameter is estimated (1) or fixed (0).}
#'     \item{\code{link_functions}: List of link functions for each parameter.}
#'     \item{\code{formulas}: List of formulas describing model structure for each parameter.}
#'     \item{\code{design_matrices}: List of design matrices generated from formulas and data.}
#'     \item{\code{model_data}: Data frame used to create design matrices.}
#'     \item{\code{sep}: String defining separator for compound parameter names.}
#'   }
#'
#' @name parameter_info
#' @aliases parameter_info-class
#' @docType class
#'
#' @seealso \code{\link{create_parameter_info}}, \code{\link{summary.parameter_info}}
#' @exportClass parameter_info
NULL


#' Create parameter information list for model fitting
#'
#' Constructs a \code{parameter_info} S3 object for likelihood-based model fitting,
#' supporting starting values, bounds, estimate flags, links, formulas, and compositional updates.
#'
#' @details
#' This function is compositional: you can use previous output as input, modifying specific slots
#' (e.g., formulas) while preserving others. The inclusion of \code{model_data} ensures reproducibility.
#'
#' The function performs the following steps:
#' \enumerate{
#'   \item Validates parameter_values as a named numeric vector or existing parameter_info object
#'   \item Initializes formulas (default ~1), link functions (default 'identity'), and bounds (default ?Inf)
#'   \item Expands bounds, flags, and other arguments to match parameter structure
#'   \item Generates design matrices from formulas and model_data
#'   \item Returns S3 object with class parameter_info
#' }
#'
#' For compositional updates, pass an existing parameter_info object as \code{parameter_values}
#' and specify only the arguments you wish to change.
#'
#' @param parameter_values Named numeric vector of initial parameter values, or an existing
#'   \code{parameter_info} object to update. Must have unique, non-empty names.
#' @param lower_bounds List or named vector of lower bounds for parameter estimates.
#'   If NULL, defaults to -Inf. Can be recycled to match parameter structure.
#' @param upper_bounds List or named vector of upper bounds for parameter estimates.
#'   If NULL, defaults to Inf. Can be recycled to match parameter structure.
#' @param estimate_flag List, vector, or logical of values (1/TRUE for estimated, 0/FALSE for fixed).
#'   If NULL, defaults to 1 (estimate all). Can be recycled to match parameter structure.
#' @param formulas List of formulas for each parameter, or a single formula to apply to all.
#'   If NULL, defaults to ~1 (intercept only).
#' @param model_data Data frame used to generate design matrices. If NULL, design matrices are NULL
#'   and a warning is issued.
#' @param link_functions List or character vector of link function names.
#'   If NULL, defaults to 'identity'. Supports character strings or function objects.
#' @param sep Separator for parameter names in composite parameter lists (default ".").
#'
#' @return An S3 object of class \code{parameter_info} with the following structure:
#'   \item{parameter_values}{Named list of parameter values (expanded by design matrix columns)}
#'   \item{lower_bounds}{Named list of lower bounds}
#'   \item{upper_bounds}{Named list of upper bounds}
#'   \item{estimate_flag}{Named list of 1/0 flags for estimation}
#'   \item{link_functions}{Named list of link functions}
#'   \item{formulas}{Named list of formulas}
#'   \item{design_matrices}{Named list of design matrices (one per parameter)}
#'   \item{model_data}{Data frame used for design matrix generation}
#'   \item{sep}{Separator for parameter names}
#'
#' @examples
#' # Example 1: Basic parameter setup with no formulas
#' df <- expand.grid(Capital = c("A", "B"), Lower = c("a", "b"))
#' info1 <- create_parameter_info(
#'   parameter_values = c(k = 2, lambda = 10, c = 0),
#'   estimate_flag = c(k = 1, lambda = 1, c = 0),
#'   lower_bounds = c(k = 0, lambda = 0, c = -1),
#'   upper_bounds = c(k = Inf, lambda = Inf, c = 1)
#' )
#'
#' # Example 2: With formulas and model data
#' info2 <- create_parameter_info(
#'   parameter_values = c(k = 2, lambda = 10, c = 0),
#'   estimate_flag = c(k = 1, lambda = 1, c = 0),
#'   formulas = list(k = ~0+Capital, lambda = ~1, c = ~Capital*Lower),
#'   model_data = df
#' )
#'
#' # Example 3: Compositional update - modify only the k formula
#' info3 <- create_parameter_info(
#'   parameter_values = info2,
#'   formulas = list(k = ~Capital + Lower)
#' )
#'
#' # Example 4: With link functions
#' info4 <- create_parameter_info(
#'   parameter_values = c(k = 2, lambda = 10, c = 0),
#'   formulas = list(k = ~Capital, lambda = ~1, c = ~1),
#'   link_functions = c("log", "log", "identity"),
#'   model_data = df
#' )
#'
#' @importFrom bstatUtils update_by_name assign_list_circular
#'
#' @export
#'
create_parameter_info <- function(
    parameter_values,
    lower_bounds = NULL,
    upper_bounds = NULL,
    estimate_flag = NULL,
    formulas = NULL,
    model_data = NULL,
    link_functions = NULL,
    sep = "."
) {

  # ============================================================================
  # Helper function: Expands arguments to match parameter list structure
  # ============================================================================
  # This function takes a scalar, vector, or list argument and expands it to
  # match the length and structure of param_list (which is the parameter values
  # list, potentially with multiple elements per parameter after design matrix
  # expansion). Returns the expanded argument as a list with the same structure.
  #
  expand_arg <- function(arg, param_list, default) {
    if (!is.null(arg)) {
      # Convert to list if not already
      arg <- as.list(arg)

      # Verify that lengths match at the parameter level
      if (length(param_list) != length(arg)) {
        stop("length of argument is not equal to the number of parameters. ",
          "Expected ", length(param_list), " but got ", length(arg),
          call. = FALSE
        )
      }

      # For each parameter, expand the argument value to match parameter length
      for (i in seq_along(param_list)) {
        if (!is.null(arg[[i]])) {
          arg[[i]] <- rep_len(arg[[i]], length(param_list[[i]]))
        } else {
          arg[[i]] <- rep_len(default, length(param_list[[i]]))
        }
      }
    } else {
      # If arg is NULL, fill with defaults for all parameters
      arg <- param_list
      for (i in seq_along(param_list)) {
        arg[[i]][] <- default
      }
    }

    return(arg)
  }

  # ============================================================================
  # Handle compositional updates: if parameter_values is a parameter_info object,
  # extract its slots and update with new arguments
  # ============================================================================
  #
  if (inherits(parameter_values, "parameter_info")) {
    # Extract parameter names from existing object
    parnames <- names(parameter_values$parameter_values)

    # Use bstatUtils::update_by_name to selectively update each component
    lower_bounds <- bstatUtils::update_by_name(
      parameter_values$lower_bounds, lower_bounds, parnames
    )
    upper_bounds <- bstatUtils::update_by_name(
      parameter_values$upper_bounds, upper_bounds, parnames
    )
    estimate_flag <- bstatUtils::update_by_name(
      parameter_values$estimate_flag, estimate_flag, parnames
    )
    formulas <- bstatUtils::update_by_name(
      parameter_values$formulas, formulas, parnames
    )
    link_functions <- bstatUtils::update_by_name(
      parameter_values$link_functions, link_functions, parnames
    )

    # Preserve model_data if not provided
    if (is.null(model_data)) {
      model_data <- parameter_values$model_data
    }

    # Preserve separator if not provided
    if (is.null(sep)) {
      sep <- parameter_values$sep
    }

    # Extract raw parameter values for further processing
      parameter_values <- parameter_values$parameter_values
  }

  # ============================================================================
  # Validate and standardize parameter_values
  # ============================================================================
  #
  param_list <- as.list(parameter_values)
  n_params <- length(param_list)
  param_names <- names(param_list)

  # Check that parameter_values is properly named
  if (is.null(param_names) || any(param_names == "") || any(duplicated(param_names))) {
    stop(
      "parameter_values must be a named vector with unique, non-empty names",
      call. = FALSE
    )
  }

  # Check that all parameter values are numeric
  if (!all(sapply(param_list, is.numeric))) {
    stop(
      "All parameter values must be numeric",
      call. = FALSE
    )
  }

  # ============================================================================
  # Handle formulas: standardize to list of formulas
  # ============================================================================
  #
  if (is.null(formulas)) {
    # Default: intercept-only formula for each parameter
    formulas <- setNames(rep(list(~1), n_params), param_names)
  } else if (inherits(formulas, "formula")) {
    # Single formula: apply to all parameters
    formulas <- setNames(rep(list(formulas), n_params), param_names)
  } else {
    # List of formulas: standardize and fill in any missing with ~1
    formulas <- as.list(formulas)
    formulas <- bstatUtils::update_by_name(
      setNames(rep(list(~1), n_params), param_names), formulas, param_names
    )
  }

  # Verify all formulas are actually formula objects
  stopifnot(all(sapply(formulas, inherits, "formula")))

  # ============================================================================
  # Handle link functions: standardize to list
  # ============================================================================
  #
  if (is.null(link_functions)) {
    # Default: identity link for all parameters
    link_functions <- setNames(rep(list("identity"), n_params), param_names)
  } else if (is.character(link_functions) && is.null(names(link_functions))) {
    # Character vector without names: recycle across parameters
    link_functions <- bstatUtils::assign_list_circular(
      link_functions, n_params, param_names
    )
  } else {
    # List or named vector: standardize and fill missing with 'identity'
    link_functions <- as.list(link_functions)
    link_functions <- bstatUtils::update_by_name(
      setNames(rep(list("identity"), n_params), param_names),
      link_functions, param_names
    )
  }

  # ============================================================================
  # Build design matrices for each parameter and expand parameter values
  # ============================================================================
  #
  if (is.null(model_data)) {
    # No design matrices if no data provided
    design_matrices <- setNames(rep(list(NULL), n_params), param_names)
    warning(
      "model_data is NULL; no design matrices are created. ",
      "Pass model_data to expand parameters via formulas.",
      call. = FALSE
    )
  } else {
    # Create design matrix for each parameter's formula
    design_matrices <- vector("list", n_params)
    for (i in seq_len(n_params)) {
      # Generate design matrix from formula and data
      design_matrices[[i]] <- model.matrix(formulas[[i]], data = model_data)
      # Clean column names: remove parentheses for consistency
      colnames(design_matrices[[i]]) <- gsub(
      "(", "", gsub(")", "",
      colnames(design_matrices[[i]]), fixed=TRUE), fixed=TRUE)

      # Expand parameter values to match number of design matrix columns
      param_list[[i]] <- rep_len(param_list[[i]], ncol(design_matrices[[i]]))
      names(param_list[[i]]) <- colnames(design_matrices[[i]])
    }

    # Assign names to design_matrices list
    names(design_matrices) <- param_names
  }

  # ============================================================================
  # Expand auxiliary arguments (bounds, flags) to parameter structure
  # ============================================================================
  #
  estimate_flag <- expand_arg(estimate_flag, param_list, 1)
  lower_bounds  <- expand_arg(lower_bounds, param_list, -Inf)
  upper_bounds  <- expand_arg(upper_bounds, param_list, Inf)

  # ============================================================================
  # Construct and return the parameter_info S3 object
  # ============================================================================
  #
  result <- list(
    parameter_values = param_list,
    lower_bounds = lower_bounds,
    upper_bounds = upper_bounds,
    estimate_flag = estimate_flag,
    link_functions = link_functions,
    formulas = formulas,
    design_matrices = design_matrices,
    model_data = model_data,
    sep = sep
  )

  # Assign parameter names to the first 7 components (all except model_data and sep)
  for (j in 1:7) {
    names(result[[j]]) <- param_names
  }

  # Assign S3 class
  class(result) <- c("parameter_info", "list")

  return(result)
}

#' Summary method for parameter_info objects
#'
#' Prints a human-readable summary of parameter information including
#' formulas with link functions and a table of parameter values, bounds,
#' and estimate flags.
#'
#' @param object A parameter_info object as returned by \code{\link{create_parameter_info}}.
#' @param ... Additional arguments (unused).
#'
#' @return Invisibly returns a list with two components:
#'   \item{info}{Data frame with parameter names and formula representations}
#'   \item{parameters}{Data frame with all parameter values, bounds, and estimate flags}
#'
#' @details
#' The returned list is structured as follows:
#' - \code{info}: Shows the formula for each parameter with its link function applied
#' - \code{parameters}: Tabular summary of all parameter-level details including
#'   base parameter name, full parameter name (after design matrix expansion),
#'   parameter values, bounds, and estimate flags
#'
#' @examples
#' df <- expand.grid(Group = c(\"A\", \"B\"))
#' info <- create_parameter_info(
#'   parameter_values = c(k = 2, lambda = 10),
#'   formulas = list(k = ~Group, lambda = ~1),
#'   model_data = df
#' )
#' summary(info)
#'
#' @export
#'
summary.parameter_info <- function(object, ...) {
  # Convert formulas to character strings for display
  formula_strs <- lapply(
    object$formulas,
    function(x) paste(as.character(x), collapse = " ")
  )

  # Create info data frame showing link functions and formulas
  info <- data.frame(
    parameter = names(object$formulas),
    formula = paste0(
      unlist(object$link_functions),
      "(", names(object$formulas), ") ",
      unlist(formula_strs)
    ),
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  # Create comprehensive parameter table by unlisting all components
  params <- as.data.frame(
    sapply(object[1:4], function(x) unlist(x), simplify = FALSE),
    stringsAsFactors = FALSE
  )

  # Add grouping columns: base parameter name and full parameter name
  params <- cbind(
    parameter = sapply(
      strsplit(row.names(params), object$sep, fixed = TRUE),
      function(x) x[1],
      simplify = TRUE
    ),
    parname = row.names(params),
    params,
    row.names = NULL
  )
  # Return as invisible list (for chaining) but print is also called
  invisible(list(info = info, parameters = params))
}
