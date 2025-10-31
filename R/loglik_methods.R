
#' Print method for loglik_control
#'
#' @param x Object of class \code{loglik_control}.
#' @param ... Additional arguments (unused).
#' @return Invisibly returns the input object.
#' @export
print.loglik_control <- function(x, ...) {
  cat("Likelihood Optimization Control\n")
  cat("===============================\n")
  cat("Method:               ", x$method, "\n")
  cat("Hessian:              ", x$hessian, "\n")
  cat("Convergence criteria: ", paste(x$criterium, collapse = ", "), "\n")
  cat("Tolerance (eps):      ", x$eps, "\n")
  cat("Max iterations:       ", x$nsteps, "\n")

  if (x$method %in% c("Levenberg", "Marquardt")) {
    cat("\nDamping parameters:\n")
    cat("  Lambda:             ", x$lambda, "\n")
    cat("  Multiplier:         ", x$lambda_multiplier, "\n")
  }

  cat("\nNumerical options:\n")
  cat("  Use ginv:           ", x$use_ginv, "\n")
  cat("  Replace NaN:        ", x$replace_nan, "\n")
  cat("  Replace Inf:        ", x$replace_inf, "\n")
  cat("  Reduce data:        ", x$reduce, "\n")

  invisible(x)
}

#' Print method for loglik fitted models
#'
#' @param x Object of class \code{loglik}.
#' @param ... Additional arguments (unused).
#' @return Invisibly returns the input object.
#' @export
print.loglik <- function(x, ...) {
  cat("Maximum Likelihood Estimation\n")
  cat("=============================\n")
  cat("Log-likelihood: ", round(x$lnlik, 4), "\n")
  cat("Converged:      ", x$converged, "\n")
  cat("Observations:   ", x$nobs, "\n")
  cat("Parameters:     ", x$npar_estimated, " estimated\n")

  cat("\nCoefficients:\n")
  print(round(x$coefficients, 4))

  invisible(x)
}

#' Summary method for loglik fitted models
#'
#' @param object Object of class \code{loglik}.
#' @param ... Additional arguments (unused).
#' @return Summary data frame with estimates, standard errors, and inference.
#' @export
summary.loglik <- function(object, ...) {
  cat("Maximum Likelihood Estimation Results\n")
  cat("======================================\n")
  cat("Log-likelihood: ", round(object$lnlik, 4), "\n")
  cat("Converged:      ", object$converged, "\n")
  cat("Observations:   ", object$nobs, "\n")
  cat("Parameters:     ", object$npar_estimated, " estimated,",
      length(object$coefficients) - object$npar_estimated, " fixed\n\n")

  cat("Parameter Estimates:\n")
  print(object$summary, digits = 4, row.names = FALSE)

  invisible(object$summary)
}

#' Extract Coefficients from loglik Object
#'
#' @param object An object of class \code{loglik}.
#' @param ... Additional arguments (unused).
#' @return Named numeric vector of parameter estimates.
#' @export
coef.loglik <- function(object, ...) {
  object$coefficients
}

#' Extract Variance-Covariance Matrix from loglik Object
#'
#' Returns the model-based variance-covariance matrix. For robust (sandwich)
#' covariance, use \code{\link{vcov_sandwich.loglik}}.
#'
#' @param object An object of class \code{loglik}.
#' @param ... Additional arguments (unused).
#' @return Variance-covariance matrix.
#' @export
vcov.loglik <- function(object, ...) {
  object$vcov
}

#' Extract Log-Likelihood from loglik Object
#'
#' @param object An object of class \code{loglik}.
#' @param ... Additional arguments (unused).
#' @return Object of class \code{logLik} with attributes for \code{df} and \code{nobs}.
#' @export
logLik.loglik <- function(object, ...) {
  val <- object$lnlik
  attr(val, "df") <- object$npar_estimated
  attr(val, "nobs") <- object$nobs
  class(val) <- "logLik"
  val
}

#' Extract Number of Observations
#'
#' @param object An object of class \code{loglik}.
#' @param ... Additional arguments (unused).
#' @return Integer, number of observations used in fitting.
#' @export
nobs.loglik <- function(object, ...) {
  object$nobs
}

#' Extract Residual Degrees of Freedom
#'
#' @param object An object of class \code{loglik}.
#' @param ... Additional arguments (unused).
#' @return Integer, residual degrees of freedom (n - p).
#' @export
df.residual.loglik <- function(object, ...) {
  object$nobs - object$npar_estimated
}

#' Extract Model Degrees of Freedom
#'
#' Returns the number of estimated parameters.
#'
#' @param object An object of class \code{loglik}.
#' @param ... Additional arguments (unused).
#' @return Integer, number of estimated parameters.
#' @export
df.model.loglik <- function(object, ...) {
  object$npar_estimated
}

#' Extract Fitted Values
#'
#' @param object An object of class \code{loglik}.
#' @param ... Additional arguments (unused).
#' @return Numeric vector of fitted values.
#' @export
fitted.loglik <- function(object, ...) {
  object$fitted_values
}

#' Extract Model Formula(s)
#'
#' @param x An object of class \code{loglik}.
#' @param ... Additional arguments (unused).
#' @return List of formulas for each parameter.
#' @export
formula.loglik <- function(x, ...) {
  x$formula
}

#' Extract Model Terms
#'
#' @param x An object of class \code{loglik}.
#' @param ... Additional arguments (unused).
#' @return List of terms objects for each parameter.
#' @export
terms.loglik <- function(x, ...) {
  x$terms
}

#' Extract Model Data
#'
#' Returns the data frame used for model fitting.
#'
#' @param object An object of class \code{loglik}.
#' @param ... Additional arguments (unused).
#' @return Data frame used in model fitting.
#' @export
model.frame.loglik <- function(object, ...) {
  object$data
}

#' Extract AIC from loglik Object
#'
#' @param object An object of class \code{loglik}.
#' @param ... Additional arguments (unused).
#' @param k Numeric, penalty per parameter (default 2 for AIC).
#' @return Numeric, AIC value.
#' @export
AIC.loglik <- function(object, ..., k = 2) {
  -2 * object$lnlik + k * object$npar_estimated
}

#' Extract BIC from loglik Object
#'
#' @param object An object of class \code{loglik}.
#' @param ... Additional arguments (unused).
#' @return Numeric, BIC value.
#' @export
BIC.loglik <- function(object, ...) {
  -2 * object$lnlik + log(object$nobs) * object$npar_estimated
}

#' Check Convergence Status
#'
#' @param object An object of class \code{loglik}.
#' @return Logical, TRUE if optimization converged.
#' @export
converged <- function(object) {
  UseMethod("converged")
}

#' @export
converged.loglik <- function(object) {
  object$converged
}

#' Extract Optimization Trace
#'
#' Returns iteration history from the optimization process.
#'
#' @param object An object of class \code{loglik}.
#' @return Matrix with columns: step, criterion, parameters, gradients, log-likelihood.
#' @export
trace_history <- function(object) {
  UseMethod("trace_history")
}

#' @export
trace_history.loglik <- function(object) {
  object$trace
}

#' Extract Hessian Matrix
#'
#' @param object An object of class \code{loglik}.
#' @return Numeric matrix, Hessian at the optimum.
#' @export
hessian <- function(object) {
  UseMethod("hessian")
}

#' @export
hessian.loglik <- function(object) {
  object$hessian
}

#' Extract Gradient/Score Matrix
#'
#' Returns the matrix of score contributions (gradient for each observation).
#'
#' @param object An object of class \code{loglik}.
#' @return Matrix (n_obs x n_params) of score contributions.
#' @export
gradient_matrix <- function(object) {
  UseMethod("gradient_matrix")
}

#' @export
gradient_matrix.loglik <- function(object) {
  object$gradient_matrix
}

#' Extract Parameter Information Structure
#'
#' @param object An object of class \code{loglik}.
#' @return Object of class \code{parameter_info}.
#' @export
parameter_info <- function(object) {
  UseMethod("parameter_info")
}

#' @export
parameter_info.loglik <- function(object) {
  object$parameter_info
}

#' Extract Control Settings
#'
#' @param object An object of class \code{loglik}.
#' @return Object of class \code{loglik_control}.
#' @export
control_settings <- function(object) {
  UseMethod("control_settings")
}

#' @export
control_settings.loglik <- function(object) {
  object$control
}

#' Predict Method for loglik Objects
#'
#' @param object An object of class \code{loglik}.
#' @param newdata Optional data frame for predictions. If NULL, uses training data.
#' @param ... Additional arguments (unused).
#' @return Numeric vector of predicted values.
#' @export
predict.loglik <- function(object, newdata = NULL, ...) {
  if (is.null(newdata)) {
    return(fitted(object))
  }

  # Reconstruct predictions using parameter_info and newdata
  # This requires recreating the linear predictors
  # Implementation depends on your model structure
  stop("predict() with newdata not yet implemented for loglik objects", call. = FALSE)
}
