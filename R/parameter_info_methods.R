#' Convert parameter_info to data frame
#'
#' Converts a \code{parameter_info} object to a data frame by computing
#' the design matrix multiplied by parameter values for each parameter.
#'
#' @details
#' This method expands parameter information using design matrices and returns
#' predicted parameter values. If the original \code{parameter_info} does not
#' have design matrices (when created without \code{model_data}), you must
#' provide a \code{data} argument.
#'
#' The result contains one column per parameter with all rows corresponding
#' to the observations in the design matrix.
#'
#' @param object A \code{parameter_info} object.
#' @param row.names Ignored (for compatibility with \code{as.data.frame} generic).
#' @param optional Ignored (for compatibility with \code{as.data.frame} generic).
#' @param data Optional data frame. If provided, design matrices are regenerated
#'   from the formulas in \code{parameter_info} using this data.
#'   If NULL (default), the \code{model_data} stored in \code{parameter_info}
#'   is used.
#' @param ... Additional arguments passed to \code{\link{as.data.frame}}.
#'
#' @return A data frame with one column per parameter. Number of rows equals
#'   the number of observations in the design matrix or provided data.
#'
#' @examples
#' df <- expand.grid(Group = c("A", "B"), Rep = 1:2)
#' info <- create_parameter_info(
#'   parameter_values = c(k = 2, lambda = 10),
#'   formulas = list(k = ~Group, lambda = ~1),
#'   model_data = df
#' )
#'
#' # Convert to data frame with stored model data
#' as.data.frame(info)
#'
#' # Convert with new data
#' new_df <- expand.grid(Group = c("A", "B", "C"), Rep = 1:2)
#' as.data.frame(info, data = new_df)
#'
#' @export
#'
as.data.frame.parameter_info <- function(object, row.names = NULL, optional = FALSE,
                                         data = NULL, ...) {

  # ============================================================================
  # Validate inputs
  # ============================================================================
  #
  if (!inherits(object, "parameter_info")) {
    stop(
      "Argument 'object' must be an object of class 'parameter_info'",
      call. = FALSE
    )
  }

  # ============================================================================
  # Handle data input: if provided, regenerate parameter_info with new data
  # ============================================================================
  #
  if (!is.null(data)) {
    # Validate data frame
    if (!is.data.frame(data)) {
      stop("'data' must be a data frame", call. = FALSE)
    }

    # Create new parameter_info with provided data
    object <- create_parameter_info(object, model_data = data)
  }

  # ============================================================================
  # Check that design matrices exist
  # ============================================================================
  #
  if (is.null(object$design_matrices[[1]])) {
    stop(
      "Design matrices not found. Provide 'data' argument or include ",
      "'model_data' when creating parameter_info",
      call. = FALSE
    )
  }

  # ============================================================================
  # Compute predicted parameter values: design matrix %*% parameter values
  # ============================================================================
  #
  parnames <- names(object$parameter_values)
  npars <- length(parnames)
  res <- NULL

  for (i in seq_len(npars)) {
    # Matrix multiplication: design matrix columns × parameter values
    pred <- object$design_matrices[[i]] %*% object$parameter_values[[i]]

    # Combine results column-wise
    res <- cbind(res, pred)
  }

  # Set column names to parameter names
  colnames(res) <- parnames

  # Convert to data frame
  as.data.frame(res)
}


#' Sample from parameter_info distribution
#'
#' Generates samples from a multivariate normal distribution centered at the
#' parameter values in a \code{parameter_info} object, using a provided
#' variance-covariance matrix.
#'
#' @details
#' This function is useful for bootstrapping and uncertainty quantification:
#' - Samples parameter values from a multivariate normal distribution
#' - Respects the \code{estimate_flag} (only samples estimated parameters)
#' - Returns a list of \code{parameter_info}-like structures for each sample
#' - Includes the original fitted parameters as the first element
#'
#' The variance-covariance matrix should correspond to the estimated parameters.
#' If the matrix is smaller than the number of parameters, it is padded with zeros
#' for fixed parameters.
#'
#' @param parameter_info A \code{parameter_info} object.
#' @param n Integer, number of samples to draw from the parameter distribution.
#'   Must be >= 1.
#' @param vcov Matrix, variance-covariance matrix for parameter estimates.
#'   Should be square and symmetric. Can be computed from a fitted model
#'   using \code{vcov()} or bootstrap methods.
#'
#' @return A list of length \code{n + 1}. Each element is a list (similar to
#'   \code{parameter_info$parameter_values}) containing sampled parameter values.
#'   The first element (named "fitted") contains the original parameter values.
#'   Subsequent elements are named "simulation_1", "simulation_2", etc.
#'
#' @examples
#' # Create parameter info
#' df <- expand.grid(Group = c("A", "B"))
#' info <- create_parameter_info(
#'   parameter_values = c(k = 2, lambda = 10, c = 0),
#'   estimate_flag = c(k = 1, lambda = 1, c = 0),
#'   model_data = df
#' )
#'
#' # Create a simple variance-covariance matrix
#' vcov_mat <- matrix(c(0.1, 0.02, 0.02, 0.15), nrow = 2)
#'
#' # Sample parameters
#' samples <- sample_parameter_info(info, n = 100, vcov = vcov_mat)
#'
#' # Access fitted parameters
#' samples$fitted
#'
#' # Access a specific simulation
#' samples$simulation_1
#'
#' @importFrom MASS mvrnorm
#'
#' @export
#'
sample_parameter_info <- function(parameter_info, n, vcov) {

  # ============================================================================
  # Validate parameter_info
  # ============================================================================
  #
  if (!inherits(parameter_info, "parameter_info")) {
    stop(
      "Argument 'parameter_info' must be an object of class 'parameter_info'",
      call. = FALSE
    )
  }

  # ============================================================================
  # Validate n
  # ============================================================================
  #
  if (!is.numeric(n) || n < 1 || n != round(n)) {
    stop("Argument 'n' must be a positive integer", call. = FALSE)
  }

  # ============================================================================
  # Build full parameter names and identify selected (estimated) parameters
  # ============================================================================
  #
  parnames <- names(parameter_info$parameter_values)
  npars <- length(parnames)

  # Full parameter names: base_name.subname format
  parnames_full <- paste(
    rep(parnames, sapply(parameter_info$parameter_values, length)),
    unlist(sapply(parameter_info$parameter_values, names)),
    sep = parameter_info$sep
  )
  npars_full <- length(parnames_full)

  # Determine which parameters are estimated (not fixed)
  cond_selected <- as.logical(unlist(parameter_info$estimate_flag))
  parnames_selected <- parnames_full[cond_selected]
  npars_selected <- sum(cond_selected)

  # ============================================================================
  # Validate variance-covariance matrix
  # ============================================================================
  #
  if (is.null(vcov) || !is.matrix(vcov) || ncol(vcov) != nrow(vcov)) {
    stop(
      "Argument 'vcov' must be a square matrix",
      call. = FALSE
    )
  }

  # Check if vcov is symmetric (or nearly symmetric for numerical stability)
  if (!isTRUE(all.equal(vcov, t(vcov)))) {
    warning(
      "vcov matrix is not symmetric; converting to symmetric form",
      call. = FALSE
    )
    vcov <- (vcov + t(vcov)) / 2
  }

  # ============================================================================
  # Expand vcov to full parameter space if necessary
  # ============================================================================
  #
  # If vcov only covers estimated parameters, pad it with zeros for fixed params
  if (ncol(vcov) < npars_full) {
    if (ncol(vcov) != npars_selected) {
      warning(
        "vcov matrix size (",  ncol(vcov), ") does not match number of ",
        "estimated parameters (", npars_selected, ")",
        call. = FALSE
      )
    }

    # Create full vcov matrix with zeros for fixed parameters
    v <- matrix(
      0, nrow = npars_full, ncol = npars_full,
      dimnames = list(parnames_full, parnames_full)
    )

    # Insert vcov values for estimated parameters
    v[cond_selected, cond_selected] <- vcov

    vcov <- v
  }

  # ============================================================================
  # Generate multivariate normal samples
  # ============================================================================
  #
  parvalues <- unlist(parameter_info$parameter_values)

  # mvrnorm returns matrix: rows are samples, columns are variables
  samples_mvn <- MASS::mvrnorm(
    n = n,
    mu = parvalues,
    Sigma = vcov,
    empirical = FALSE
  )

  # Add fitted values as first row
  parvalues <- rbind(parvalues, samples_mvn)
  rownames(parvalues) <- c("fitted", paste("simulation", seq_len(n), sep = "_"))

  # ============================================================================
  # Restructure samples back into parameter_info-like lists
  # ============================================================================
  #
  # Create index map: which columns belong to which parameter
  subnames <- unlist(sapply(parameter_info$parameter_values, names))
  groups <- factor(
    rep(parnames, sapply(parameter_info$parameter_values, length)),
    levels = parnames
  )
  indices <- setNames(seq_len(npars_full), subnames)
  indices <- split(indices, groups)

  # Convert each sample to the parameter_info structure
  L <- vector("list", n + 1)

  for (i in seq_len(n + 1)) {
    L_i <- indices

    # Assign sampled values to each parameter group
    for (j in seq_along(indices)) {
      L_i[[j]][] <- parvalues[i, indices[[j]]]
    }

    L[[i]] <- L_i
  }

  # Assign names to the list
  names(L) <- rownames(parvalues)

  return(L)
}


#' Extract parameter values
#'
#' Extracts the parameter values from a \code{parameter_info} object
#' as a named vector or list.
#'
#' @details
#' \code{coef.parameter_info} returns parameter values in flattened form
#' (as a named vector), similar to \code{\link{coef}} for model objects.
#'
#' Optionally filters to only return estimated parameters (those with
#' \code{estimate_flag = 1}).
#'
#' @param object A \code{parameter_info} object.
#' @param estimated_only Logical. If TRUE, only return parameters with
#'   \code{estimate_flag = 1}. Default is FALSE (return all parameters).
#' @param ... Additional arguments (unused).
#'
#' @return A named numeric vector of parameter values.
#'   Names follow the format "base_param.subparam" when design matrices
#'   are present.
#'
#' @examples
#' info <- create_parameter_info(
#'   parameter_values = c(k = 2, lambda = 10, c = 0),
#'   estimate_flag = c(k = 1, lambda = 1, c = 0)
#' )
#'
#' # All parameters
#' coef(info)
#'
#' # Only estimated parameters
#' coef(info, estimated_only = TRUE)
#'
#' @export
#'
coef.parameter_info <- function(object, estimated_only = FALSE, ...) {

  if (!inherits(object, "parameter_info")) {
    stop(
      "Argument 'object' must be an object of class 'parameter_info'",
      call. = FALSE
    )
  }

  # Flatten parameter values to a named vector
  pars <- unlist(object$parameter_values)

  # Filter to estimated parameters if requested
  if (estimated_only) {
    estimate_flags <- unlist(object$estimate_flag)
    pars <- pars[as.logical(estimate_flags)]
  }

  return(pars)
}


#' Extract variance-covariance information
#'
#' Extracts bounds and uncertainty information from a \code{parameter_info} object.
#'
#' @details
#' Returns a list containing lower bounds, upper bounds, and estimate flags
#' for all parameters. Useful for optimization or constraint checking.
#'
#' @param object A \code{parameter_info} object.
#' @param ... Additional arguments (unused).
#'
#' @return A list with components:
#'   \item{lower}{Lower bounds as named vector}
#'   \item{upper}{Upper bounds as named vector}
#'   \item{fixed}{Logical vector indicating fixed parameters}
#'   \item{estimated}{Logical vector indicating estimated parameters}
#'
#' @examples
#' info <- create_parameter_info(
#'   parameter_values = c(k = 2, lambda = 10, c = 0),
#'   estimate_flag = c(k = 1, lambda = 1, c = 0),
#'   lower_bounds = c(k = 0, lambda = 0, c = -1),
#'   upper_bounds = c(k = 100, lambda = 100, c = 1)
#' )
#'
#' vcov(info)
#'
#' @export
#'
vcov.parameter_info <- function(object, ...) {

  if (!inherits(object, "parameter_info")) {
    stop(
      "Argument 'object' must be an object of class 'parameter_info'",
      call. = FALSE
    )
  }

  estimate_flags <- unlist(object$estimate_flag)

  list(
    lower = unlist(object$lower_bounds),
    upper = unlist(object$upper_bounds),
    fixed = as.logical(1 - estimate_flags),
    estimated = as.logical(estimate_flags)
  )
}


#' Extract link functions
#'
#' Retrieves the link functions applied to each parameter.
#'
#' @details
#' Link functions describe transformations applied to linear predictors
#' to obtain parameter values. Common links include "identity", "log", "logit".
#'
#' @param object A \code{parameter_info} object.
#' @param ... Additional arguments (unused).
#'
#' @return A named list where each element is the link function
#'   for the corresponding parameter.
#'
#' @examples
#' info <- create_parameter_info(
#'   parameter_values = c(k = 2, lambda = 10, c = 0),
#'   link_functions = c("log", "log", "identity")
#' )
#'
#' links(info)
#'
#' @export
#'
links <- function(object, ...) {
  UseMethod("links")
}

#' @export
links.parameter_info <- function(object, ...) {

  if (!inherits(object, "parameter_info")) {
    stop(
      "Argument 'object' must be an object of class 'parameter_info'",
      call. = FALSE
    )
  }

  return(object$link_functions)
}


#' Extract formulas
#'
#' Retrieves the model formulas used to generate design matrices
#' for each parameter.
#'
#' @details
#' Formulas describe the structure of how parameters vary across
#' observations. For example, ~1 means intercept-only (constant parameter),
#' while ~Group means the parameter varies by levels of Group variable.
#'
#' @param object A \code{parameter_info} object.
#' @param ... Additional arguments (unused).
#'
#' @return A named list where each element is a formula for the
#'   corresponding parameter.
#'
#' @examples
#' df <- expand.grid(Group = c("A", "B"))
#' info <- create_parameter_info(
#'   parameter_values = c(k = 2, lambda = 10),
#'   formulas = list(k = ~Group, lambda = ~1),
#'   model_data = df
#' )
#'
#' formulas(info)
#'
#' @export
#'
formulas <- function(object, ...) {
  UseMethod("formulas")
}

#' @export
formulas.parameter_info <- function(object, ...) {

  if (!inherits(object, "parameter_info")) {
    stop(
      "Argument 'object' must be an object of class 'parameter_info'",
      call. = FALSE
    )
  }

  return(object$formulas)
}


#' Extract design matrices
#'
#' Retrieves the design matrices generated from model formulas
#' for each parameter.
#'
#' @details
#' Design matrices are used to expand parameter values across observations
#' based on the model formulas. Particularly useful when you want to inspect
#' how formulas translate to the data structure.
#'
#' @param object A \code{parameter_info} object.
#' @param which Optional parameter name(s) to extract specific design matrices.
#'   If NULL (default), returns design matrices for all parameters.
#' @param ... Additional arguments (unused).
#'
#' @return If \code{which} is NULL, a named list of design matrices.
#'   Otherwise, a design matrix or list of design matrices for the
#'   specified parameter(s).
#'
#' @examples
#' df <- expand.grid(Group = c("A", "B"), Rep = 1:2)
#' info <- create_parameter_info(
#'   parameter_values = c(k = 2, lambda = 10),
#'   formulas = list(k = ~Group, lambda = ~1),
#'   model_data = df
#' )
#'
#' # All design matrices
#' design_matrices(info)
#'
#' # Specific parameter
#' design_matrices(info, which = "k")
#'
#' @export
#'
design_matrices <- function(object, which = NULL, ...) {
  UseMethod("design_matrices")
}

#' @export
design_matrices.parameter_info <- function(object, which = NULL, ...) {

  if (!inherits(object, "parameter_info")) {
    stop(
      "Argument 'object' must be an object of class 'parameter_info'",
      call. = FALSE
    )
  }

  if (is.null(which)) {
    # Return all design matrices
    return(object$design_matrices)
  } else {
    # Validate which parameter names
    valid_params <- names(object$design_matrices)
    if (!all(which %in% valid_params)) {
      stop(
        "Invalid parameter names in 'which'. Valid parameters are: ",
        paste(valid_params, collapse = ", "),
        call. = FALSE
      )
    }

    # Return requested design matrices
    if (length(which) == 1) {
      return(object$design_matrices[[which]])
    } else {
      return(object$design_matrices[which])
    }
  }
}


#' Print method for parameter_info
#'
#' Displays a concise summary of a \code{parameter_info} object.
#'
#' @details
#' Shows parameter names, their values, estimate flags, and bounds.
#' Provides a quick overview of the parameter configuration.
#'
#' @param object A \code{parameter_info} object.
#' @param ... Additional arguments (unused).
#'
#' @return Invisibly returns \code{object}.
#'
#' @examples
#' info <- create_parameter_info(
#'   parameter_values = c(k = 2, lambda = 10, c = 0),
#'   estimate_flag = c(k = 1, lambda = 1, c = 0),
#'   lower_bounds = c(k = 0, lambda = 0, c = -1)
#' )
#'
#' print(info)
#'
#' @export
#'
print.parameter_info <- function(object, ...) {

  if (!inherits(object, "parameter_info")) {
    stop(
      "Argument 'object' must be an object of class 'parameter_info'",
      call. = FALSE
    )
  }

  cat("\n=== parameter_info object ===\n\n")

  # Display parameter names and values
  cat("Parameters (", length(object$parameter_values), " total):\n", sep = "")

  pars <- unlist(object$parameter_values)
  estimate_flags <- unlist(object$estimate_flag)
  lower <- unlist(object$lower_bounds)
  upper <- unlist(object$upper_bounds)

  df_print <- data.frame(
    Parameter = names(pars),
    Value = as.numeric(pars),
    Lower = as.numeric(lower),
    Upper = as.numeric(upper),
    Estimate = ifelse(as.logical(estimate_flags), "Yes", "No"),
    row.names = NULL
  )

  print(df_print, digits = 4)

  cat("\nSeparator: \"", object$sep, "\"\n", sep = "")

  # Show number of estimated vs fixed
  n_est <- sum(as.logical(estimate_flags))
  n_fixed <- length(estimate_flags) - n_est

  cat(
    "  Estimated parameters: ", n_est, "\n",
    "  Fixed parameters: ", n_fixed, "\n",
    sep = ""
  )

  # Show design matrix info if present
  if (!is.null(object$design_matrices[[1]])) {
    cat("\nDesign matrices present for all parameters\n")
  } else {
    cat("\nNo design matrices (model_data not provided)\n")
  }

  invisible(object)
}
