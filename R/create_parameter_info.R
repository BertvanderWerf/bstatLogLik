#' Parameter Information Structure
#'
#' S3 class defining the structure to hold model parameter details for likelihood-based fitting.
#'
#' Objects of this class contain:
#'   \itemize{
#'     \item{\code{parameter_values}: Named list of initial parameter values used in modeling.}
#'     \item{\code{lower_bounds}: Lower bounds for parameter estimates.}
#'     \item{\code{upper_bounds}: Upper bounds for parameter estimates.}
#'     \item{\code{estimate_flag}: Logical or numeric vector (1/0) indicating whether each parameter is estimated (1) or fixed (0).}
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
#' @seealso \code{\link{create_parameter_info}}
#' @exportClass parameter_info

#' Create parameter information list for model fitting
#'
#' Constructs a \code{parameter_info} S3 object for likelihood-based model fitting,
#' supporting starting values, bounds, estimate flags, links, formulas, and compositional updates.
#'
#' This function is compositional: you can use previous output as input, modifying specific slots
#' (e.g., formulas) while preserving others. The inclusion of \code{model_data} ensures reproducibility.
#'
#' @param parameter_values Named numeric vector of initial parameter values, or an existing \code{parameter_info} object to update.
#' @param lower_bounds List of lower bounds for parameter estimates (or NULL for \code{-Inf}).
#' @param upper_bounds List of upper bounds for parameter estimates (or NULL for \code{Inf}).
#' @param estimate_flag List or vector of 1s and 0s; 1 for estimated, 0 for fixed.
#' @param formulas List of formulas for each parameter, or NULL.
#' @param model_data Data frame used to generate design matrices.
#' @param link_functions List of link functions, or NULL.
#' @param sep Separator for parameter names (default ".").
#' @return An S3 object of class \code{parameter_info}, see \code{\link{parameter_info}} for structure.
#' @examples
#' df <- expand.grid(Capital = c("A", "B"), Lower = c("a", "b"))
#' info1 <- create_parameter_info(
#'   parameter_values = c(k = 2, lambda = 10, c = 0),
#'   estimate_flag = c(k = 1, lambda = 1, c = 0),
#'   formulas = list(k = ~0+Capital, lambda = ~1, c = ~Capital*Lower),
#'   model_data = df
#' )
#' # Update only the formula for 'k', keeping other slots
#' info2 <- create_parameter_info(
#'   parameter_values = info1,
#'   formulas = list(k = ~Capital + Lower)
#' )
#' @export
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
  # Helper expands any list or vector (bounds, flags etc.) to required length/naming
  expand_arg <- function(arg, param_list, default) {
    if (!is.null(arg)) {
      arg <- as.list(arg)
      if (length(param_list) != length(arg)) {
        stop('length of ', deparse(substitute(arg)), ' is not equal to the number of parameters')
      }
      for (i in seq_along(param_list)) {
        if (!is.null(arg[[i]])) {
          arg[[i]] <- rep_len(arg[[i]], length(param_list[[i]]))
        } else {
          arg[[i]] <- rep_len(default, length(param_list[[i]]))
        }
      }
    } else {
      arg <- param_list
      for (i in seq_along(param_list)) {
        arg[[i]][] <- default
      }
    }
    arg
  }

  # If input is a parameter_info object, extract underlying info and update with new arguments
  if (inherits(parameter_values, "parameter_info")) {
    parnames <- names(parameter_values$parameter_values)
    lower_bounds   <- bstatUtils::update_by_name(parameter_values$lower_bounds,   lower_bounds,   parnames)
    upper_bounds   <- bstatUtils::update_by_name(parameter_values$upper_bounds,   upper_bounds,   parnames)
    estimate_flag  <- bstatUtils::update_by_name(parameter_values$estimate_flag,  estimate_flag,  parnames)
    formulas       <- bstatUtils::update_by_name(parameter_values$formulas,       formulas,       parnames)
    link_functions <- bstatUtils::update_by_name(parameter_values$link_functions, link_functions, parnames)
    if (is.null(model_data)) model_data <- parameter_values$model_data
    if (is.null(sep)) sep <- parameter_values$sep
    parameter_values <- parameter_values$parameter_values
  }

  # Standardize parameter values
  param_list <- as.list(parameter_values)
  n_params <- length(param_list)
  param_names <- names(param_list)
  if (is.null(param_names) || any(param_names == '') || any(duplicated(param_names))) {
    stop("parameter_values must be a named vector with unique, non-empty names", call. = FALSE)
  }

  # Handle formulas
  if (is.null(formulas)) {
    formulas <- setNames(rep(list(~1), n_params), param_names)
  } else if (inherits(formulas, "formula")) {
    formulas <- setNames(rep(list(formulas), n_params), param_names)
  } else {
    formulas <- as.list(formulas)
    formulas <- bstatUtils::update_by_name(setNames(rep(list(~1), n_params), param_names), formulas, param_names)
  }
  stopifnot(all(sapply(formulas, inherits, "formula")))

  # Standardize link functions
  if (is.null(link_functions)) {
    link_functions <- rep(list('identity'), n_params)
  } else if (is.character(link_functions)) {
    link_functions <- bstatUtils::assign_list_circular(link_functions, n_params, param_names)
  } else {
    link_functions <- as.list(link_functions)
    link_functions <- bstatUtils::update_by_name(setNames(rep(list('identity'), n_params), param_names), link_functions, param_names)
  }

  # Build design matrices for each parameter
  if (is.null(model_data)) {
    design_matrices <- rep(list(NULL), n_params)
    warning('model_data in parameter_info is NULL; no design matrices are created.')
  } else {
    design_matrices <- vector("list", n_params)
    for (i in seq_len(n_params)) {
      design_matrices[[i]] <- model.matrix(formulas[[i]], data = model_data)
      # Clean colnames for consistency
      colnames(design_matrices[[i]]) <- gsub("(", "", gsub(")", "", colnames(design_matrices[[i]]), fixed=TRUE), fixed=TRUE)
      param_list[[i]] <- rep_len(param_list[[i]], ncol(design_matrices[[i]]))
      names(param_list[[i]]) <- colnames(design_matrices[[i]])
    }
  }

  # Expand auxiliary arguments to parameter lengths
  estimate_flag <- expand_arg(estimate_flag, param_list, 1)
  lower_bounds  <- expand_arg(lower_bounds, param_list, -Inf)
  upper_bounds  <- expand_arg(upper_bounds, param_list, Inf)

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

  # Assign parameter names to lists
  for (j in seq_along(result)[1:7]) names(result[[j]]) <- param_names

  class(result) <- c('parameter_info', class(result))
  result
}

#' Summary for "parameter_info" object
#'
#' Prints formula representation and tabulated parameter values.
#'
#' @param object A parameter_info object as returned by \code{create_parameter_info()}.
#' @return A list with info and parameter details for inspection or printing.
#' @export
summary.parameter_info <- function(object) {
  # Convert formulas to string for easier display
  formula_strs <- lapply(object$formulas, function(x) paste(as.character(x), collapse = ' '))
  info <- data.frame(parameter = names(object$formulas),
                     formula = paste0(unlist(object$link_functions), "(", names(object$formulas), ") ", unlist(formula_strs)))
  # Collate parameter table
  params <- as.data.frame(sapply(object[1:4], function(x) unlist(x), simplify = FALSE))
  params <- cbind(
    parameter = sapply(strsplit(row.names(params), ".", fixed = TRUE), function(x) x[1], simplify = TRUE),
    parname = row.names(params),
    params
  )
  row.names(params) <- NULL
  list(info = info, parameters = params)
}
