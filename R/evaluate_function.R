#' Evaluate a Model Function Dynamically
#'
#' Validates inputs, builds dynamic environments, applies link functions, differentiates if required, and evaluates for a given data frame.
#'
#' @param fun The function to evaluate (e.g. model or cdf).
#' @param vars Optional named character vector of variable mapping ("model_var" = "data_column").
#' @param parameter_info Metadata created by `create_parameter_info()`, containing parameter values, design matrices, links, formulas.
#' @param data Data frame containing evaluation data.
#'
#' @return Vector or list of results for the supplied data.
#' @export
evaluate_function <- function(
    fun,
    vars = NULL,
    parameter_info,
    data=NULL
) {
  # Error and type checks for robustness
  bstatErr::check_function(fun)
  if (!inherits(parameter_info, "parameter_info")) {
    stop("parameter 'parameter_info' is not of class 'parameter_info'.", call. = FALSE)
  }
  bstatErr::check_data_frame(data, allow_null = TRUE)
  if (is.null(data)) {
    if (is.null(parameter_info$model_data)) {
      stop("Both 'data' and 'parameter_info$model_data' cannot be NULL.", call. = FALSE)
    }
    data <- parameter_info$model_data
  } else {
    # Recreate parameter_info with data if needed
    parameter_info <- create_parameter_info(parameter_info, model_data = data)
  }
  bstatErr::check_string_vector(vars, allow_null = TRUE, must_have_names = TRUE)

  if (!is.null(vars)) {
    not_in_data <- vars[!(vars %in% colnames(data))]
    if (length(not_in_data) > 0) {
      stop(sprintf("Variables not found in data: %s", paste(not_in_data, collapse = ", ")), call. = FALSE)
    }
    vars <- setNames(vars, names(vars))
  } else {
    vars <- character(0)
  }

  param_names <- names(parameter_info$parameter_values)
  n_params <- length(param_names)

  # Extract and process formula terms, fixed vars
  formulas <- parameter_info$formulas
  fixed_vars <- c()
  for (i in seq_len(n_params)) {
    fml <- formulas[[i]]
    if (!is.null(fml)) {
      t <- stats::terms(fml)
      fixed_vars <- c(fixed_vars, labels(t))
    }
  }
  fixed_vars <- unique(unlist(strsplit(fixed_vars, ":")))
  fixed_vars <- fixed_vars[!(fixed_vars %in% c("0", "1"))]
  vars <- c(vars, fixed_vars[!(fixed_vars %in% vars)])

  # Symbolic differentiation and function construction
  func_body <- paste(deparse(Deriv::Deriv(body(fun), nderiv = 0, cache.exp = FALSE)), collapse = "\n")
  func_text <- sprintf("function() { local({ %s}) }", func_body)

  # Parameter usage check
  unused <- param_names[!sapply(param_names, function(x) grepl(sprintf("\\b%s\\b", x), func_text))]
  if (length(unused) > 0) {
    stop(sprintf("Model/parameter mismatch: %d of %d parameters not used: %s", length(unused), n_params, paste(unused, collapse = ",")), call. = FALSE)
  }

  # Insert link/inverse functions
  for (i in seq_len(n_params)) {
    if (parameter_info$link_functions[[i]] != "identity") {
      inv_fun <- bstatUtils::get_inverse_function(parameter_info$link_functions[[i]])
      if (is.null(inv_fun)) {
        stop(sprintf('%s is not found in the inverse_function_table, please use bstatUtils::add_inverse_function() to add it.',
                     parameter_info$link_functions[[i]]), call. = FALSE)
      }
      func_text <- gsub(sprintf("\\b%s\\b", param_names[i]), sprintf("%s(%s)", inv_fun, param_names[i]), func_text)
    }
  }

  # Build evaluable function
  func <- eval(parse(text = func_text))
  environment(func) <- new.env()

  # add default values for parameters added 11-11-2025
  fun_args <- formals(args(fun))
  has_value <- !sapply(fun_args, rlang::is_missing)
  if (any(has_value)) {
    for (i in seq_len(length(has_value))) {
      if (has_value[i]) {
        assign(names(has_value)[i],
               eval(fun_args[[i]], environment(func)),
               environment(func))
      }
    }
  }


  # Bind parameters (multiply through design matrix if present)
  params <- parameter_info$parameter_values
  for (i in seq_len(n_params)) {
    if (!is.null(parameter_info$design_matrices[[i]])) {
      params[[i]] <- as.vector(params[[i]] %*% t(parameter_info$design_matrices[[i]]))
    }
    assign(param_names[i], params[[i]], envir = environment(func))
    if (param_names[i] %in% names(has_value)) {
      has_value[param_names[i]] <- TRUE
    }
  }

  # Bind variables from data
  for (i in 1:length(vars)) {
    if(is.null(names(vars))) {
      name=vars[i]
    } else if (names(vars)[i]=='') {
      name=vars[i]
    } else {
      name=names(vars)[i]
    }
    assign(name, data[,vars[i]], envir=environment(func))
    if (name %in% names(has_value)) {
      has_value[name] <- TRUE
    }
  }

  if (!(all(has_value))) {
    stop(sprintf("the variable or parameter %s needed in the function has not been defined yet.",
                 paste0("'", paste(names(has_value)[!has_value], collapse="', '"), "'")),
         call. = FALSE)
  }

  # Evaluate and return result
  func()
}
