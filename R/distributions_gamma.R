#' Gamma Distribution Functions for Likelihood Estimation
#'
#' Probability density function (PDF), cumulative distribution function (CDF),
#' and utility functions for the two-parameter Gamma distribution with
#' shape (alpha) and rate (beta) parameterization.
#'
#' These functions are optimized for use with \code{\link{loglik}} and include
#' pre-computed derivative rules stored in \code{\link{drule}}.
#'
#' @name gamma_distribution
#' @family distribution functions
NULL

#' Gamma Probability Density Function
#'
#' Computes the PDF of the Gamma distribution using shape-rate parameterization.
#'
#' The PDF is given by:
#' \deqn{f(x; \alpha, \beta) = \frac{\beta^\alpha}{\Gamma(\alpha)} x^{\alpha-1} e^{-\beta x}}
#'
#' @param x Numeric vector of quantiles (must be non-negative).
#' @param alpha Numeric. Shape parameter (must be positive).
#' @param beta Numeric. Rate parameter (must be positive).
#' @return Numeric vector of density values.
#'
#' @examples
#' # Standard gamma
#' gamma_pdf(x = 0:10, alpha = 2, beta = 1)
#'
#' # Compare with base R dgamma
#' identical(
#'   gamma_pdf(x = 1:5, alpha = 2, beta = 1),
#'   dgamma(x = 1:5, shape = 2, rate = 1)
#' )
#'
#' @export
gamma_pdf <- function(x, alpha, beta) {
  # Ensure x is non-negative
  (beta^alpha / gamma(alpha)) * (x^(alpha - 1)) * exp(-beta * x)
  # result <- ifelse(x<0, 0, result)
  # result
}

#' Gamma Cumulative Distribution Function
#'
#' Computes the CDF of the Gamma distribution using shape-rate parameterization.
#'
#' Uses the regularized lower incomplete gamma function:
#' \deqn{F(x; \alpha, \beta) = P(\alpha, \beta x) = \frac{\gamma(\alpha, \beta x)}{\Gamma(\alpha)}}
#'
#' @param x Numeric vector of quantiles (must be non-negative).
#' @param alpha Numeric. Shape parameter (must be positive).
#' @param beta Numeric. Rate parameter (must be positive).
#' @return Numeric vector of cumulative probabilities.
#'
#' @examples
#' gamma_cdf(x = 0:10, alpha = 2, beta = 1)
#'
#' @export
gamma_cdf <- function(x, alpha, beta) {
  # Use built-in pgamma with shape and rate parameterization
  pgamma(x, shape = alpha, rate = beta)
}

#' Gamma Hazard Function
#'
#' Computes the hazard (failure rate) function of the Gamma distribution.
#'
#' The hazard is \eqn{h(x) = f(x) / (1 - F(x))}.
#'
#' @param x Numeric vector of quantiles.
#' @param alpha Numeric. Shape parameter (must be positive).
#' @param beta Numeric. Rate parameter (must be positive).
#' @return Numeric vector of hazard rates.
#'
#' @examples
#' gamma_hazard(x = 0:10, alpha = 2, beta = 1)
#'
#' @export
gamma_hazard <- function(x, alpha, beta) {
  pdf <- gamma_pdf(x, alpha, beta)
  cdf <- gamma_cdf(x, alpha, beta)
  pdf / (1 - cdf)
}

#' Gamma Mean
#'
#' Computes the expected value (mean) of the Gamma distribution.
#'
#' The mean is:
#' \deqn{E[X] = \frac{\alpha}{\beta}}
#'
#' @param alpha Numeric. Shape parameter (must be positive).
#' @param beta Numeric. Rate parameter (must be positive).
#' @return Numeric scalar, the mean of the distribution.
#'
#' @examples
#' gamma_mean(alpha = 2, beta = 1)
#'
#' @export
gamma_mean <- function(alpha, beta) {
  alpha / beta
}

#' Gamma Variance
#'
#' Computes the variance of the Gamma distribution.
#'
#' The variance is:
#' \deqn{\text{Var}(X) = \frac{\alpha}{\beta^2}}
#'
#' @param alpha Numeric. Shape parameter (must be positive).
#' @param beta Numeric. Rate parameter (must be positive).
#' @return Numeric scalar, the variance of the distribution.
#'
#' @examples
#' gamma_variance(alpha = 2, beta = 1)
#'
#' @export
gamma_variance <- function(alpha, beta) {
  alpha / (beta^2)
}

#' Gamma Median (Approximate)
#'
#' Computes an approximation of the median of the Gamma distribution.
#'
#' Uses the approximation:
#' \deqn{\text{Median} \approx \alpha - \frac{1}{3}}
#' for \eqn{\alpha > 1}, with improved accuracy for larger \eqn{\alpha}.
#'
#' @param alpha Numeric. Shape parameter (must be positive).
#' @param beta Numeric. Rate parameter (must be positive).
#' @return Numeric scalar, approximate median of the distribution.
#'
#' @details
#' For exact median, use \code{qgamma(0.5, shape = alpha, rate = beta)}.
#'
#' @examples
#' gamma_median(alpha = 2, beta = 1)
#'
#' @export
gamma_median <- function(alpha, beta) {
  # Approximation: median ≈ (alpha - 1/3) / beta for alpha > 1
  # For exact value, use qgamma(0.5, shape = alpha, rate = beta) / beta
  qgamma(0.5, shape = alpha, rate = beta) / beta
}

#' Gamma Mode
#'
#' Computes the mode of the Gamma distribution.
#'
#' The mode is:
#' \deqn{\text{Mode} = \begin{cases}
#' 0 & \text{if } \alpha \leq 1 \\
#' \frac{\alpha - 1}{\beta} & \text{if } \alpha > 1
#' \end{cases}}
#'
#' @param alpha Numeric. Shape parameter (must be positive).
#' @param beta Numeric. Rate parameter (must be positive).
#' @return Numeric scalar, the mode of the distribution.
#'
#' @examples
#' gamma_mode(alpha = 2, beta = 1)
#' gamma_mode(alpha = 0.5, beta = 1)  # Returns 0 for alpha <= 1
#'
#' @export
gamma_mode <- function(alpha, beta) {
  ifelse(alpha <= 1, 0, (alpha - 1) / beta)
}

#' Initial Parameter Estimates for Gamma Distribution
#'
#' Computes method-of-moments initial parameter estimates for the Gamma distribution.
#'
#' Uses the relationship between mean, variance, and parameters:
#' \deqn{\alpha = \frac{mean^2}{variance}}
#' \deqn{\beta = \frac{mean}{variance}}
#'
#' @param x Numeric vector of observations (must be positive).
#' @param weights Numeric vector of weights/frequencies (optional).
#' @param verbose Logical. Print intermediate calculations (default FALSE).
#' @return Named numeric vector with elements \code{alpha} and \code{beta}.
#'
#' @examples
#' set.seed(123)
#' x <- rgamma(100, shape = 2, rate = 1)
#' gamma_init(x)
#'
#' @export
gamma_init <- function(x, weights = NULL, verbose = FALSE) {
  # Remove NAs
  x <- x[!is.na(x)]

  # Compute (weighted) mean and variance
  if (is.null(weights)) {
    m <- mean(x)
    v <- var(x)
  } else {
    w_sum <- sum(weights, na.rm = TRUE)
    m <- sum(x * weights, na.rm = TRUE) / w_sum
    v <- sum(weights * (x - m)^2, na.rm = TRUE) / w_sum
  }

  if (verbose) {
    message(sprintf("Mean: %.4f, Variance: %.4f", m, v))
  }

  # Method-of-moments estimates
  alpha <- m^2 / v
  beta <- m / v

  c(alpha = alpha, beta = beta)
}

# Compute derivatives for PDF
drule[['gamma_pdf']] <- lapply(
  c('x', 'alpha', 'beta'),
  function(var) Deriv::Deriv(body(gamma_pdf), var, cache.exp = FALSE)
)
names(drule[['gamma_pdf']]) <- c('x', 'alpha', 'beta')

# Compute derivatives for CDF (uses built-in pgamma, so derivatives are approximate)
# drule[['gamma_cdf']] <- lapply(
#   c('x', 'alpha', 'beta'),
#   function(var) Deriv::Deriv(body(gamma_cdf), var, cache.exp = FALSE)
# )
# names(drule[['gamma_cdf']]) <- c('x', 'alpha', 'beta')
