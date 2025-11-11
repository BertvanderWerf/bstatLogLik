#' Fitted Maximum Likelihood Model
#'
#' S3 class containing results from maximum likelihood estimation via \code{\link{loglik}}.
#'
#' Objects of this class contain:
#'   \itemize{
#'     \item{\code{fun}: The original likelihood function provided by the user.}
#'     \item{\code{control}: The \code{loglik_control} object used for optimization.}
#'     \item{\code{vars}: Character vector of variable names used in the model.}
#'     \item{\code{n}: Name of frequency/count variable (if applicable).}
#'     \item{\code{formula}: List of formulas from \code{parameter_info}.}
#'     \item{\code{terms}: List of parsed formula terms objects.}
#'     \item{\code{fixed_vars}: Character vector of fixed effect variables extracted from formulas.}
#'     \item{\code{data}: Data frame subset used for estimation.}
#'     \item{\code{nobs}: Total number of observations (accounting for frequency weights if present).}
#'     \item{\code{converged}: Logical indicating successful convergence.}
#'     \item{\code{lnlik}: Final maximized log-likelihood value.}
#'     \item{\code{npar_estimated}: Number of parameters actually estimated (non-fixed).}
#'     \item{\code{trace}: Matrix of iteration history (step, criterion, parameters, gradients, log-likelihood).}
#'     \item{\code{vcov}: Variance-covariance matrix (model-based) for all parameters.}
#'     \item{\code{summary}: Data frame with parameter estimates, variances, bounds, confidence intervals, and p-values.}
#'     \item{\code{coefficients}: Named vector of final parameter estimates.}
#'     \item{\code{parameter_info}: The \code{parameter_info} object with final parameter values.}
#'     \item{\code{gradient_matrix}: Matrix of score contributions (n_obs x n_params) for sandwich estimation.}
#'     \item{\code{hessian}: Final Hessian matrix at the optimum.}
#'     \item{\code{fitted_values}: Vector of fitted values from the model.}
#'     \item{\code{derivs}: List containing derivative evaluations from the final iteration.}
#' \item{\code{boundary_parameters}: (if any) Data frame identifying parameters that converged to bounds.}
#'   }
#'
#' @name loglik_fit
#' @aliases loglik-class loglik_fit-class
#' @docType class
#'
#' @seealso \code{\link{loglik}}, \code{\link{check_boundary_parameters}}, \code{\link{estfun.loglik}}, \code{\link{vcov_sandwich.loglik}}
#' @exportClass loglik
NULL

#' Maximum Likelihood Estimation with Sandwich Covariance
#'
#' Fits parameters via Newton-Raphson, Levenberg, or Marquardt optimization,
#' with automatic differentiation, parameter constraints, and robust standard errors.
#'
#' @details
#' The optimization algorithm solves the likelihood equations iteratively using
#' one of three methods (Newton-Raphson, Levenberg-Marquardt, Marquardt) with
#' convergence assessed using specified criteria (parameter change, gradient norm,
#' log-likelihood change, or combinations thereof).
#'
#' **Parameter bounds**: Parameters are automatically constrained to their specified
#' bounds during optimization. If a parameter converges to its bound, this indicates
#' potential issues with model specification or data identification - a warning will
#' be issued. Use \code{\link{check_boundary_parameters}} for detailed inspection.
#'
#' **Link functions**: Parameters can be transformed using link functions (e.g., "log",
#' "logit") specified in \code{parameter_info}. This allows fitting parameters on
#' different scales.
#'
#' **Design matrices**: If \code{parameter_info} contains formulas, design matrices
#' expand parameters across factor levels or covariate values automatically.
#'
#' **Sandwich covariance**: Robust (sandwich) standard errors are computed from the
#' gradient (score) matrix. Use \code{\link{vcov_sandwich.loglik}} for robust inference.
#'
#' @param fun Function. Likelihood function (not log-likelihood).
#' @param fun_vars Named character vector. Name(s) of primary variables for function;
#'   names should be column names in data.
#' @param parameter_info Object of class \code{parameter_info} from \code{create_parameter_info()}.
#' @param data Data frame or matrix with all analysis variables, or NULL to use
#'   parameter_info$model_data.
#' @param trace Logical. If TRUE, print iteration progress (default: TRUE).
#' @param control List from \code{\link{loglik_control}} with optimization settings.
#' @param n Character or NULL. Name of frequency/count variable (if data are aggregated).
#'
#' @return Object of class "loglik" with the fitted model results and diagnostics.
#'   See \code{\link{loglik_fit}} for detailed component descriptions.
#'
#' @note
#' **Boundary warnings**: If any parameter converges to its lower or upper bound,
#' a warning is issued. This typically indicates:
#' \itemize{
#'   \item Bounds may be too restrictive
#'   \item Parameters may not be well-identified
#'   \item Model formulation may need review
#' }
#' Inspect using \code{\link{check_boundary_parameters}} for details.
#'
#' **Missing fit information**: If convergence fails or numerical issues occur,
#' some results (e.g., Hessian) may be invalid. Check \code{$converged} status.
#'
#' @references
#' For model selection: Burnham, K.P. and Anderson, D.R. (2002) Model selection
#' and multimodel inference: a practical information-theoretic approach. 2nd ed.
#' Springer-Verlag, New York. DOI: 10.1007/b97636
#'
#' @examples
#' \dontrun{
#' # Fit a Weibull model with group-specific shape parameters
#' data <- data.frame(
#'   x = rweibull(100, shape = 2, scale = 5),
#'   group = rep(c("A", "B"), each = 50)
#' )
#'
#' param_info <- create_parameter_info(
#'   parameter_values = c(k = 2, lambda = 5),
#'   estimate_flag = c(k = 1, lambda = 1),
#'   lower_bounds = c(k = 0.01, lambda = 0.01),
#'   upper_bounds = c(k = 100, lambda = 100),
#'   formulas = list(k = ~0 + group, lambda = ~1),
#'   model_data = data
#' )
#'
#' fit <- loglik(
#'   fun = dweibull,
#'   fun_vars = c(x = "x"),
#'   parameter_info = param_info,
#'   data = data,
#'   trace = TRUE
#' )
#'
#' summary(fit)
#' check_boundary_parameters(fit)
#' }
#'
#' @export
loglik <- function(fun, fun_vars=NULL, parameter_info, data, trace = TRUE,
                   control = loglik_control(), n = NULL) {

  result <- list()

  # ============================================================================
  # INPUT VALIDATION
  # ============================================================================

  # Validate fun
  if (!is.function(fun)) {
    stop("Argument 'fun' must be a function", call. = FALSE)
  }
  result$fun <- fun

  # Validate parameter_info
  if (!inherits(parameter_info, 'parameter_info')) {
    stop("Argument 'parameter_info' must be an object of class 'parameter_info'",
         call. = FALSE)
  }

  # Handle data input
  if (is.null(data)) {
    if (is.null(parameter_info$model_data)) {
      stop('data must be set or model_data in parameter_info must be set',
           call. = FALSE)
    } else {
      data <- parameter_info$model_data
    }
  }
  bstatErr::check_data_frame(data)

  # Validate fun_vars
  bstatErr::check_string_vector(fun_vars, allow_null=TRUE, must_have_names=TRUE)
  if (!is.null(fun_vars)) {
    if (!all(fun_vars %in% names(data))) {
      stop("All elements of 'fun_vars' must be valid column names in 'data'.",
           call. = FALSE)
    }
  }
  result$vars <- fun_vars

  # Validate trace
  bstatErr::check_logical(trace)

  # Validate and set control
  if (!is.null(control)) {
    if (!inherits(control, 'loglik_control')) {
      stop('control must be a loglik_control', call. = FALSE)
    }
  } else {
    control <- loglik_control()
  }
  result$control <- control

  # Check frequency variable
  if (!is.null(n)) {
    n <- bstatErr::check_string(n)
    if (!(n %in% names(data))) {
      stop(sprintf("Argument '%s' must be a valid column name of argument 'data'.", n),
           call. = FALSE)
    }
    bstatErr::check_numeric_vector(data[,n])
    result$n <- n
    result$vars <- c(result$vars, setNames(n, 'n'))
  }

  # ============================================================================
  # RECREATE parameter_info WITH CURRENT DATA
  # ============================================================================

  parameter_info <- create_parameter_info(parameter_info, model_data = data)

  npars <- length(parameter_info$parameter_values)
  parnames <- names(parameter_info$parameter_values)

  # Extract fixed variables from formulas
  result$formula <- parameter_info$formulas
  result$terms <- result$formula
  result$fixed_vars <- c()

  for (i in seq_len(npars)) {
    if (!is.null(result$terms[[i]])) {
      result$terms[[i]] <- stats::terms(result$terms[[i]])
      result$fixed_vars <- c(result$fixed_vars, labels(result$terms[[i]]))
    }
  }
  result$fixed_vars <- unique(unlist(strsplit(result$fixed_vars, ":")))
  result$fixed_vars <- result$fixed_vars[!(result$fixed_vars %in% c("0", "1"))]

  result$vars <- c(result$vars, setNames(result$fixed_vars,result$fixed_vars))
  result$data <- data[, unique(result$vars), drop = FALSE]

  # ============================================================================
  # PARAMETER SETUP
  # ============================================================================

  parnames_full <- names(unlist(parameter_info$parameter_values))
  npars_full <- length(parnames_full)
  pars_cond <- as.logical(unlist(parameter_info$estimate_flag))
  parnames_selected <- parnames_full[pars_cond]
  npars_selected <- length(parnames_selected)

  par <- unlist(parameter_info$parameter_values)[pars_cond]
  par_min <- unlist(parameter_info$lower_bounds)[pars_cond]
  par_max <- unlist(parameter_info$upper_bounds)[pars_cond]

  # Enforce bounds on initial parameters
  parameter_info$parameter_values <- lapply(seq_len(npars), function(i) {
    pmin(pmax(parameter_info$parameter_values[[i]],
              parameter_info$lower_bounds[[i]]),
         parameter_info$upper_bounds[[i]])
  })
  names(parameter_info$parameter_values) <- parnames

  # ============================================================================
  # DATA REDUCTION (optional)
  # ============================================================================

  if (control$reduce) {
    if (!all(pars_cond)) {
      mat <- parameter_info$design_matrices
      mat[sapply(mat, is.null)] <- matrix(1, 1, 1)
      rcond <- rowSums(as.data.frame(mat)[, pars_cond, drop = FALSE]) > 0

      if (!all(rcond)) {
        if (!any(rcond)) {
          stop("No data available for estimating selected parameters")
        }
        result$data <- result$data[rcond, ]
        for (i in seq_len(npars)) {
          if (!is.null(parameter_info$design_matrices[[i]])) {
            parameter_info$design_matrices[[i]] <- parameter_info$design_matrices[[i]][rcond, , drop = FALSE]
          }
        }
      }
    }
  }

  # ============================================================================
  # COMPUTE NUMBER OF OBSERVATIONS
  # ============================================================================

  if (is.null(n) || all(data[[n]] == 1)) {
    n <- NULL
    result$nobs <- nrow(result$data)
  } else {
    result$nobs <- sum(result$data[[n]])
  }

# ============================================================================
  # CONSTRUCT LOG-LIKELIHOOD FUNCTION
  # ============================================================================

  if (is.null(n)) {
    ln_lik_i <- paste0("function() { log(",
                       paste(deparse(Deriv::Deriv(body(fun), "x", nderiv = 0, cache.exp = FALSE)), collapse = ""),
                       ")} " )
  } else {
    ln_lik_i <- paste0("function() { n * log(",
                       paste(deparse(Deriv::Deriv(body(fun), "x", nderiv = 0, cache.exp = FALSE)), collapse = ""),
                       ")}")
  }

  # Insert link functions if needed
  for (i in seq_len(npars)) {
    if (parameter_info$link_functions[[i]] != "identity") {
      inv_fun <- bstatUtils::get_inverse_function(parameter_info$link_functions[[i]])
      if (is.null(inv_fun)) {
        stop(sprintf('%s is not found in the inverse_function_table, please use bstatUtils::add_inverse_function() to add it.',
                     parameter_info$link_functions[[i]]), call. = FALSE)
      }
      ln_lik_i <- gsub(sprintf("\\b%s\\b", parnames[i]),
                       sprintf("%s(%s)", inv_fun, parnames[i]),
                       ln_lik_i)


    }
  }
  result$fun_text <- ln_lik_i
  ln_lik_i <- eval(parse(text = ln_lik_i))

  # --- Determine which parameters to estimate ---
  cond <- sapply(parameter_info$estimate_flag, sum) > 0
  nderiv <- 0:1
  if (control$hessian == "calculate") nderiv <- c(nderiv, 2)

  # Compute derivatives
  derivs <- Deriv::Deriv(f = ln_lik_i, x = parnames[cond],
                         cache.exp = TRUE, nderiv = nderiv, combine = "cbind")

  environment(derivs) <- new.env()

  # add default values for parameters added 11-11-2025
  fun_args <- formals(args(fun))
  has_value <- !sapply(fun_args, rlang::is_missing)
  if (any(has_default)) {
    for (i in seq_len(length(has_value))) {
      if (has_value[i]) {
        assign(names(has_value)[i],
               eval(fun_args[[i]], environment(derivs)),
               environment(derivs))
      }
    }
  }

  # Add parameters to derivative environment
  pars <- parameter_info$parameter_values
  for (i in seq_len(npars)) {
    if (!is.null(parameter_info$design_matrices[[i]])) {
      pars[[i]] <- c(pars[[i]] %*% t(parameter_info$design_matrices[[i]]))
    }
    assign(parnames[i], pars[[i]], envir = environment(derivs))
    if (parnames[i] %in% names(has_value)) {
      has_value[parnames[i]] <- TRUE
    }
  }

  # Add variables to derivative environment
  for (i in seq_along(result$vars)) {
    var_name <- if (is.null(names(result$vars)) || names(result$vars)[i] == "") {
      result$vars[i]
    } else {
      names(result$vars)[i]
    }
    assign(var_name, result$data[[result$vars[i]]], envir = environment(derivs))
    if (var_name %in% names(has_value)) {
      has_value[var_name] <- TRUE
    }
  }

  if (!(all(has_value))) {
    stop(sprintf("the variable or parameter %s needed in the function has not been defined yet.",
                 paste0("'", paste(names(has_value)[!has_value], collapse="', '"), "'")),
         call. = FALSE)
  }

  # --- Initial evaluation ---
  deriv_output <- derivs()
  result$derivs <- deriv_output
  deriv_output <- adjust_nan_inf(deriv_output, control)

  loglik <- sum(deriv_output$`0`)

  gradient_matrix <- compute_gradient(deriv_output, parameter_info, cond, npars, parnames_selected, control)
  jacobian_vec <- colSums(gradient_matrix, na.rm = TRUE)
  hessian <- compute_hessian(deriv_output, gradient_matrix, parameter_info, cond, npars, parnames_selected, control)

  # --- Optimization loop ---
  crit <- Inf
  step <- 0
  trace_data <- NULL
  loglik_old <- NA

  while (crit > control$eps && step < control$nsteps) {
    step <- step + 1

    # Solve for parameter update
    if (control$use_ginv) {
      delta <- solve_general(hessian, jacobian_vec)
    } else {
      delta <- solve(hessian, jacobian_vec)
    }
    delta <- adjust_nan_inf(delta, control)

    stopifnot(all(is.finite(delta)))

    # Update parameters with bounds enforcement
    par_old <- par
    par <- par - delta
    par <- pmax(pmin(par, par_max), par_min)
    delta <- par - par_old

    # Update parameter_info structure
    k <- 0
    for (i in seq_len(npars)) {
      for (j in seq_along(parameter_info$estimate_flag[[i]])) {
        if (parameter_info$estimate_flag[[i]][j] == 1) {
          k <- k + 1
          parameter_info$parameter_values[[i]][j] <- par[k]
        }
      }
    }

    # Update parameter vectors in derivative environment
    pars <- parameter_info$parameter_values
    for (i in seq_len(npars)) {
      if (!is.null(parameter_info$design_matrices[[i]])) {
        pars[[i]] <- c(pars[[i]] %*% t(parameter_info$design_matrices[[i]]))
      }
      assign(parnames[i], pars[[i]], envir = environment(derivs))
    }

    # Re-evaluate derivatives
    deriv_output <- derivs()
    result$derivs <- deriv_output
    deriv_output <- adjust_nan_inf(deriv_output, control)

    loglik <- sum(deriv_output$`0`)

    gradient_matrix <- compute_gradient(deriv_output, parameter_info, cond, npars, parnames_selected, control)
    jacobian_vec <- colSums(gradient_matrix, na.rm = TRUE)
    hessian <- compute_hessian(deriv_output, gradient_matrix, parameter_info, cond, npars, parnames_selected, control)

    # Compute convergence criterion
    crit_old <- crit
    crit <- Inf

    if ("delta" %in% control$criterium) {
      crit <- min(crit, sqrt(sum(delta^2)))
    }
    if ("jacobian" %in% control$criterium) {
      crit <- min(crit, sqrt(sum(jacobian_vec^2)))
    }
    if ("lnlik" %in% control$criterium) {
      if (step > 1) {
        crit <- min(crit, abs(loglik_old - loglik))
      }
      loglik_old <- loglik
    }

    # Trace output
    if (trace) {
      cat_par <- par
      cond_boundary <- (par == par_min) | (par == par_max)
      if (any(cond_boundary)) {
        cat_par[cond_boundary] <- paste0("\033[31m", cat_par[cond_boundary], "\033[0m")
      }

      if (step == 1) {
        trace_data <- t(c(step, crit, par, jacobian_vec, loglik))
        colnames(trace_data) <- c("step", "crit", names(par),
                                  paste("gradient", names(par), sep = "."), "loglik")
        cat(colnames(trace_data), "\n")
      } else {
        trace_data <- rbind(trace_data, c(step, crit, par, jacobian_vec, loglik))
      }
      cat(step, crit, cat_par, jacobian_vec, loglik, "\n")
    }

    # Lambda adjustment for Levenberg-Marquardt
    if (control$method != "NewtonRaphson") {
      if (crit > crit_old) {
        control$lambda <- control$lambda * control$lambda_multiplier
        par <- par_old
        crit <- crit_old
      } else {
        control$lambda <- control$lambda / control$lambda_multiplier
      }
    }
  }

  # --- Convergence check ---
  result$converged <- TRUE
  if (crit > control$eps && step == control$nsteps) {
    cat("No convergence within", control$nsteps, "steps\n")
    result$converged <- FALSE
  }

  # --- Covariance matrix ---
  cov <- matrix(0, nrow = npars_full, ncol = npars_full,
                dimnames = list(parnames_full, parnames_full))

  if (control$use_ginv) {
    cov[pars_cond, pars_cond] <- -solve_general(hessian)
  } else {
    cov[pars_cond, pars_cond] <- -solve(hessian)
  }

  var <- diag(cov)

  # --- Summary table ---
  result$parameters <- as.data.frame(pars)
  pars_full <- unlist(parameter_info$parameter_values)
  parname_factor <- factor(
    unlist(lapply(seq_len(npars), function(i) {
      rep(names(parameter_info$parameter_values)[i], length(parameter_info$parameter_values[[i]]))
    })),
    names(parameter_info$parameter_values)
  )

  summary_df <- data.frame(
    parname = parname_factor,
    estimate = pars_full,
    variance = var,
    is_constant = !pars_cond,
    min = unlist(parameter_info$lower_bounds),
    max = unlist(parameter_info$upper_bounds)
  )

  # Add confidence intervals and p-values if biostats package available
  if (requireNamespace("biostats", quietly = TRUE)) {
    summary_df <- biostats::addUpperLower(summary_df, df = result$nobs - sum(pars_cond))
    summary_df <- biostats::addPvalue(summary_df, df = result$nobs - sum(pars_cond))
  }

  # --- Fitted values ---
  if (is.null(n)) {
    result$fitted_values <- exp(deriv_output$`0`)
  } else {
    result$fitted_values <- exp(deriv_output$`0` / result$data[[n]])
  }

  if (any(is.infinite(result$fitted_values))) {
    warning("Iterations not optimal: expected values ± Inf, see fitted_values")
  }

  # ============================================================================
  # CHECK FOR PARAMETERS AT BOUNDARIES
  # ============================================================================
  # Detect parameters that converged to their specified bounds
  # This is important diagnostic information
  #
  pars_full <- unlist(parameter_info$parameter_values)
  par_min_full <- unlist(parameter_info$lower_bounds)
  par_max_full <- unlist(parameter_info$upper_bounds)

  # Numerical tolerance for boundary detection
  tolerance <- 1e-6

  # Identify which parameters are at boundaries
  at_lower <- abs(pars_full - par_min_full) < tolerance & !is.infinite(par_min_full)
  at_upper <- abs(pars_full - par_max_full) < tolerance & !is.infinite(par_max_full)
  at_boundary <- at_lower | at_upper

  # Issue warning if any parameters at boundaries
  if (any(at_boundary)) {
    boundary_names <- names(pars_full)[at_boundary]
    boundary_sides <- character(length(boundary_names))
    boundary_sides[at_lower[at_boundary]] <- "lower"
    boundary_sides[at_upper[at_boundary]] <- "upper"

    warning(
      length(boundary_names), " parameter(s) converged to bounds: ",
      paste(paste0(boundary_names, " (", boundary_sides, ")"), collapse=", "),
      ". This may indicate insufficient data or model over-specification.",
      call. = FALSE
    )

    # Store boundary information in result for later inspection
    result$boundary_parameters <- data.frame(
      parameter = boundary_names,
      bound_type = boundary_sides,
      estimate = pars_full[at_boundary],
      bound_value = ifelse(at_lower[at_boundary],
                           par_min_full[at_boundary],
                           par_max_full[at_boundary]),
      stringsAsFactors = FALSE
    )
  } else {
    # No boundary parameters - make this accessible
    result$boundary_parameters <- NULL
  }

  # ============================================================================
  # RETURN RESULTS
  # ============================================================================

  result$lnlik <- loglik
  result$npar_estimated <- npars_selected
  result$trace <- trace_data
  result$vcov <- cov
  result$summary <- summary_df
  result$coefficients <- unlist(parameter_info$parameter_values)
  result$parameter_info <- parameter_info
  result$gradient_matrix <- gradient_matrix
  result$hessian <- hessian

  class(result) <- c("loglik", class(result))
  return(result)
}
