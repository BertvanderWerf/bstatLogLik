#' Log-Normal Distribution Functions for Likelihood Estimation
#'
#' Probability density function (PDF), cumulative distribution function (CDF),
#' and utility functions for the two-parameter Log-Normal distribution with
#' meanlog (mu) and sdlog (sigma) parameterization.
#'
#' These functions are optimized for use with \code{\link{loglik}} and include
#' pre-computed derivative rules stored in \code{\link{drule}}.
#'
#' @name lognormal_distribution
#' @family distribution functions
NULL

#' Log-Normal Probability Density Function
#'
#' Computes the PDF of the Log-Normal distribution.
#'
#' The PDF is given by:
#' \deqn{f(x; \mu, \sigma) = \frac{1}{x\sigma\sqrt{2\pi}} \exp\left(-\frac{(\log x - \mu)^2}{2\sigma^2}\right)}
#'
#' for \eqn{x > 0}.
#'
#' @param x Numeric vector of quantiles (must be positive).
#' @param mu Numeric. Mean of log-transformed data.
#' @param sigma Numeric. Standard deviation of log-transformed data (must be positive).
#' @return Numeric vector of density values.
#'
#' @examples
#' lognormal_pdf(x = 1:10, mu = 0, sigma = 1)
#'
#' @export
lognormal_pdf <- function(x, mu, sigma) {
  # Log-Normal PDF
  (1 / (x * sigma * sqrt(2 * pi))) * exp(-0.5 * ((log(x) - mu) / sigma)^2)
  #result <- ifelse(x <= 0, 0, result)
  #result
}

#' Log-Normal Cumulative Distribution Function
#'
#' Computes the CDF of the Log-Normal distribution.
#'
#' The CDF is given by:
#' \deqn{F(x; \mu, \sigma) = \Phi\left(\frac{\log x - \mu}{\sigma}\right)}
#'
#' where \eqn{\Phi} is the standard normal CDF.
#'
#' @param x Numeric vector of quantiles (must be positive).
#' @param mu Numeric. Mean of log-transformed data.
#' @param sigma Numeric. Standard deviation of log-transformed data (must be positive).
#' @return Numeric vector of cumulative probabilities.
#'
#' @examples
#' lognormal_cdf(x = 1:10, mu = 0, sigma = 1)
#'
#' @export
lognormal_cdf <- function(x, mu, sigma) {
  # Use built-in plnorm
  plnorm(x, meanlog = mu, sdlog = sigma)
}

#' Log-Normal Quantile Function
#'
#' Computes the quantile (inverse CDF) of the Log-Normal distribution.
#'
#' @param p Numeric vector of probabilities.
#' @param mu Numeric. Mean of log-transformed data.
#' @param sigma Numeric. Standard deviation of log-transformed data (must be positive).
#' @return Numeric vector of quantiles.
#'
#' @examples
#' lognormal_quantile(p = c(0.025, 0.5, 0.975), mu = 0, sigma = 1)
#'
#' @export
lognormal_quantile <- function(p, mu, sigma) {
  qlnorm(p, meanlog = mu, sdlog = sigma)
}

#' Log-Normal Hazard Function
#'
#' Computes the hazard (failure rate) function of the Log-Normal distribution.
#'
#' The hazard is \eqn{h(x) = f(x) / (1 - F(x))}.
#'
#' @param x Numeric vector of quantiles.
#' @param mu Numeric. Mean of log-transformed data.
#' @param sigma Numeric. Standard deviation of log-transformed data (must be positive).
#' @return Numeric vector of hazard rates.
#'
#' @examples
#' lognormal_hazard(x = 1:10, mu = 0, sigma = 1)
#'
#' @export
lognormal_hazard <- function(x, mu, sigma) {
  pdf <- lognormal_pdf(x, mu, sigma)
  cdf <- lognormal_cdf(x, mu, sigma)
  pdf / (1 - cdf)
}

#' Log-Normal Mean
#'
#' Computes the expected value (mean) of the Log-Normal distribution.
#'
#' The mean is:
#' \deqn{E[X] = \exp\left(\mu + \frac{\sigma^2}{2}\right)}
#'
#' @param mu Numeric. Mean of log-transformed data.
#' @param sigma Numeric. Standard deviation of log-transformed data (must be positive).
#' @return Numeric scalar, the mean of the distribution.
#'
#' @examples
#' lognormal_mean(mu = 0, sigma = 1)
#'
#' @export
lognormal_mean <- function(mu, sigma) {
  exp(mu + 0.5 * sigma^2)
}

#' Log-Normal Variance
#'
#' Computes the variance of the Log-Normal distribution.
#'
#' The variance is:
#' \deqn{\text{Var}(X) = \left(e^{\sigma^2} - 1\right) e^{2\mu + \sigma^2}}
#'
#' @param mu Numeric. Mean of log-transformed data.
#' @param sigma Numeric. Standard deviation of log-transformed data (must be positive).
#' @return Numeric scalar, the variance of the distribution.
#'
#' @examples
#' lognormal_variance(mu = 0, sigma = 1)
#'
#' @export
lognormal_variance <- function(mu, sigma) {
  exp(sigma^2) - 1
}

#' Log-Normal Median
#'
#' Computes the median of the Log-Normal distribution.
#'
#' The median is:
#' \deqn{\text{Median} = e^\mu}
#'
#' @param mu Numeric. Mean of log-transformed data.
#' @return Numeric scalar, the median of the distribution.
#'
#' @examples
#' lognormal_median(mu = 0)
#'
#' @export
lognormal_median <- function(mu) {
  exp(mu)
}

#' Log-Normal Mode
#'
#' Computes the mode of the Log-Normal distribution.
#'
#' The mode is:
#' \deqn{\text{Mode} = \exp(\mu - \sigma^2)}
#'
#' @param mu Numeric. Mean of log-transformed data.
#' @param sigma Numeric. Standard deviation of log-transformed data (must be positive).
#' @return Numeric scalar, the mode of the distribution.
#'
#' @examples
#' lognormal_mode(mu = 0, sigma = 1)
#'
#' @export
lognormal_mode <- function(mu, sigma) {
  exp(mu - sigma^2)
}

#' Initial Parameter Estimates for Log-Normal Distribution
#'
#' Computes initial parameter estimates for the Log-Normal distribution
#' using method-of-moments on log-transformed data.
#'
#' The transformation is:
#' \deqn{\mu = \text{mean}(\log x)}
#' \deqn{\sigma = \text{sd}(\log x)}
#'
#' @param x Numeric vector of observations (must be positive).
#' @param weights Numeric vector of weights/frequencies (optional).
#' @param verbose Logical. Print intermediate calculations (default FALSE).
#' @return Named numeric vector with elements \code{mu} and \code{sigma}.
#'
#' @examples
#' set.seed(123)
#' x <- rlnorm(100, meanlog = 0, sdlog = 1)
#' lognormal_init(x)
#'
#' @export
lognormal_init <- function(x, weights = NULL, verbose = FALSE) {
  # Remove NAs and ensure positive
  x <- x[!is.na(x) & x > 0]

  if (length(x) == 0) {
    stop("No valid positive observations", call. = FALSE)
  }

  # Compute (weighted) mean and SD of log-transformed data
  log_x <- log(x)

  if (is.null(weights)) {
    mu <- mean(log_x)
    sigma <- sd(log_x)
  } else {
    w_sum <- sum(weights, na.rm = TRUE)
    mu <- sum(log_x * weights, na.rm = TRUE) / w_sum
    sigma <- sqrt(sum(weights * (log_x - mu)^2, na.rm = TRUE) / w_sum)
  }

  if (verbose) {
    message(sprintf("Mean of log(x): %.4f, SD of log(x): %.4f", mu, sigma))
  }

  c(mu = mu, sigma = sigma)
}

# Compute derivatives for PDF
drule[['lognormal_pdf']] <- lapply(
  c('x', 'mu', 'sigma'),
  function(var) Deriv::Deriv(body(lognormal_pdf), var, cache.exp = FALSE)
)
names(drule[['lognormal_pdf']]) <- c('x', 'mu', 'sigma')

# Compute derivatives for CDF
# drule[['lognormal_cdf']] <- lapply(
#   c('x', 'mu', 'sigma'),
#   function(var) Deriv::Deriv(body(lognormal_cdf), var, cache.exp = FALSE)
# )
# names(drule[['lognormal_cdf']]) <- c('x', 'mu', 'sigma')
