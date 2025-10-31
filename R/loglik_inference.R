#' Likelihood Ratio Test for Nested Models
#'
#' Performs likelihood ratio test comparing two or more nested models fitted
#' with \code{\link{loglik}}. Models must be nested (one is a special case of another).
#'
#' The test statistic is:
#' \deqn{LR = -2(\ell_0 - \ell_1)}
#' where \eqn{\ell_0} is the log-likelihood of the restricted (simpler) model
#' and \eqn{\ell_1} is the log-likelihood of the unrestricted (more complex) model.
#' Under the null hypothesis (that the restrictions are true), LR follows a
#' chi-squared distribution with degrees of freedom equal to the difference
#' in number of parameters.
#'
#' @param object First fitted model (of class \code{loglik}).
#' @param ... Additional fitted models or arguments.
#' @param test Character. Always uses "LRT" for likelihood ratio test.
#' @return Data frame with test results including log-likelihoods, df, test statistic, and p-value.
#'
#' @details
#' Models are automatically ordered from simplest (fewest parameters) to most complex.
#' Each model is compared to the previous (simpler) one.
#'
#' **Important assumptions:**
#' \itemize{
#'   \item{Models must be nested (parameters in simpler model are subset of complex model)}
#'   \item{Models must be fit to the same data}
#'   \item{Models must have converged successfully}
#'   \item{Asymptotic chi-squared approximation must be valid (large sample, interior of parameter space)}
#' }
#'
#' @examples
#' \dontrun{
#' # Fit nested models
#' pars1 <- create_parameter_info(c(k = 2, lambda = 5), formulas = list(k = ~1, lambda = ~1))
#' pars2 <- create_parameter_info(c(k = 2, lambda = 5), formulas = list(k = ~group, lambda = ~1))
#'
#' fit1 <- loglik(weibull_pdf, "time", pars1, data)
#' fit2 <- loglik(weibull_pdf, "time", pars2, data)
#'
#' # Likelihood ratio test
#' lrt(fit1, fit2)
#' }
#'
#' @export
lrt <- function(object, ...) {
  UseMethod("lrt")
}

#' @export
lrt.loglik <- function(object, ..., test = "LRT") {
  # Collect all models
  models <- list(object, ...)
  models <- models[sapply(models, inherits, "loglik")]

  if (length(models) < 2) {
    stop("At least two models are required for likelihood ratio test", call. = FALSE)
  }

  # Extract key information
  get_model_info <- function(model) {
    list(
      loglik = model$lnlik,
      df = model$npar_estimated,
      nobs = model$nobs,
      converged = model$converged
    )
  }

  model_info <- lapply(models, get_model_info)

  # Check that all models converged
  converged <- sapply(model_info, `[[`, "converged")
  if (!all(converged)) {
    warning("Some models did not converge. Results may be unreliable.")
  }

  # Check same number of observations
  nobs_vec <- sapply(model_info, `[[`, "nobs")
  if (length(unique(nobs_vec)) > 1) {
    stop("All models must be fit to the same data (same number of observations)", call. = FALSE)
  }

  # Order models by degrees of freedom (simplest to most complex)
  df_vec <- sapply(model_info, `[[`, "df")
  order_idx <- order(df_vec)
  models <- models[order_idx]
  model_info <- model_info[order_idx]
  df_vec <- df_vec[order_idx]

  # Build results table
  n_models <- length(models)
  loglik_vec <- sapply(model_info, `[[`, "loglik")

  result <- data.frame(
    Model = paste0("Model ", seq_along(models)),
    Df = df_vec,
    LogLik = loglik_vec,
    Df_diff = c(NA, diff(df_vec)),
    LR_stat = c(NA, -2 * diff(loglik_vec)),
    Pr_Chi = rep(NA, n_models)
  )

  # Compute p-values
  for (i in 2:n_models) {
    df_diff <- result$Df_diff[i]
    lr_stat <- result$LR_stat[i]

    if (df_diff > 0 && lr_stat >= 0) {
      result$Pr_Chi[i] <- pchisq(lr_stat, df = df_diff, lower.tail = FALSE)
    } else if (df_diff == 0) {
      warning("Models ", i-1, " and ", i, " have same df; not nested?")
      result$Pr_Chi[i] <- NA
    } else if (lr_stat < 0) {
      warning("Negative LR statistic for comparison ", i-1, " vs ", i,
              "; simpler model fits better. Models may not be nested.")
      result$Pr_Chi[i] <- NA
    }
  }

  # Format output
  class(result) <- c("lrt", "data.frame")
  attr(result, "heading") <- "Likelihood Ratio Test for Nested Models"
  attr(result, "nobs") <- nobs_vec[1]

  result
}

#' Print Method for LRT Results
#'
#' @param x Object of class \code{lrt}.
#' @param digits Number of digits for printing (default 4).
#' @param ... Additional arguments (unused).
#' @export
print.lrt <- function(x, digits = 4, ...) {
  cat(attr(x, "heading"), "\n")
  cat("Number of observations:", attr(x, "nobs"), "\n\n")

  # Format output
  x_print <- x
  x_print$LogLik <- round(x_print$LogLik, digits)
  x_print$LR_stat <- round(x_print$LR_stat, digits)
  x_print$Pr_Chi <- format.pval(x_print$Pr_Chi, digits = digits)

  print.data.frame(x_print, row.names = FALSE)

  cat("\nNote: LR_stat = -2 * (LogLik_simple - LogLik_complex)\n")
  cat("      P-values from chi-squared distribution with Df_diff degrees of freedom\n")

  invisible(x)
}

#' ANOVA Method for loglik Objects (Likelihood Ratio Test)
#'
#' S3 method for \code{anova()} that performs likelihood ratio tests
#' on nested models fitted with \code{\link{loglik}}.
#'
#' This is a wrapper around \code{\link{lrt}} that provides familiar
#' \code{anova()} syntax while explicitly using likelihood ratio tests.
#'
#' @param object First fitted model (of class \code{loglik}).
#' @param ... Additional fitted models.
#' @param test Character. Must be "LRT" or "Chisq" (both perform LRT).
#' @return Object of class \code{lrt} containing test results.
#'
#' @examples
#' \dontrun{
#' anova(fit1, fit2, test = "LRT")
#' }
#'
#' @export
anova.loglik <- function(object, ..., test = "LRT") {
  if (!missing(test) && !(test %in% c("LRT", "Chisq"))) {
    stop("Only likelihood ratio test (test = 'LRT' or 'Chisq') is supported for loglik objects",
         call. = FALSE)
  }

  lrt.loglik(object, ..., test = test)
}

#' Compare Multiple Models with Information Criteria
#'
#' Compares models using AIC, BIC, and optionally likelihood ratio tests.
#'
#' @param ... Fitted models of class \code{loglik}.
#' @param lrt Logical. Perform likelihood ratio tests if models are nested (default TRUE).
#' @return Data frame with model comparison statistics.
#'
#' @examples
#' \dontrun{
#' compare_models(fit1, fit2, fit3)
#' }
#'
#' @export
compare_models <- function(..., lrt = TRUE) {
  models <- list(...)
  models <- models[sapply(models, inherits, "loglik")]

  if (length(models) < 2) {
    stop("At least two models required for comparison", call. = FALSE)
  }

  # Extract information
  model_names <- paste0("Model ", seq_along(models))
  if (!is.null(names(list(...)))) {
    named_idx <- names(list(...)) != ""
    model_names[named_idx] <- names(list(...))[named_idx]
  }

  result <- data.frame(
    Model = model_names,
    Df = sapply(models, function(x) x$npar_estimated),
    LogLik = sapply(models, function(x) x$lnlik),
    AIC = sapply(models, AIC),
    BIC = sapply(models, BIC),
    Converged = sapply(models, function(x) x$converged)
  )

  # Add delta AIC and BIC
  result$Delta_AIC <- result$AIC - min(result$AIC)
  result$Delta_BIC <- result$BIC - min(result$BIC)

  # Optionally add LRT
  if (lrt) {
    cat("Model Comparison Summary\n")
    cat("========================\n\n")
    print(result, row.names = FALSE)
    cat("\n")

    tryCatch({
      cat("Likelihood Ratio Tests:\n")
      cat("-----------------------\n")
      lrt_result <- do.call(lrt.loglik, models)
      print(lrt_result)
    }, error = function(e) {
      cat("Could not perform LRT (models may not be nested)\n")
    })
  }

  invisible(result)
}
