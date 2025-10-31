#' Distribution Family Functions for Likelihood Estimation
#'
#' This package provides a collection of distribution functions optimized
#' for use with \code{\link{loglik}} maximum likelihood estimation.
#'
#' Available distributions:
#' \itemize{
#'   \item{\code{\link{weibull_distribution}}: Three-parameter Weibull}
#'   \item{\code{\link{gamma_distribution}}: Two-parameter Gamma}
#'   \item{\code{\link{normal_distribution}}: Two-parameter Normal}
#'   \item{\code{\link{lognormal_distribution}}: Two-parameter Log-normal}
#' }
#'
#' Each distribution family includes:
#' \itemize{
#'   \item{PDF and CDF functions}
#'   \item{Hazard function (where applicable)}
#'   \item{Summary statistics (mean, median, mode)}
#'   \item{Initial parameter estimators}
#'   \item{Pre-computed derivative rules for automatic differentiation}
#' }
#'
#' @section Usage with loglik:
#'
#' Distribution functions are designed to work seamlessly with \code{\link{loglik}}.
#' The derivative rules stored in \code{\link{drule}} enable automatic differentiation:
#'
#' \preformatted{
#' # Example: Fit Weibull to survival data
#' pars <- create_parameter_info(
#'   parameter_values = weibull_init(data$time),
#'   formulas = list(k = ~1, lambda = ~group, c = ~0),
#'   model_data = data
#' )
#'
#' fit <- loglik(fun = weibull_pdf, x = "time", parameter_info = pars, data = data)
#' }
#'
#' @name distributions
#' @docType package
#' @keywords internal
"_PACKAGE"

#' @keywords internal
drule <- list()

# Note: Individual drule components are populated by each distribution module.
# See distributions_weibull.R, distributions_gamma.R, etc.
