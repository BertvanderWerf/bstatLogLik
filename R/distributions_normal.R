#' Normal Distribution Functions for Likelihood Estimation
#'
#' Probability density function (PDF), cumulative distribution function (CDF),
#' and utility functions for the two-parameter Normal (Gaussian) distribution
#' with mean (mu) and standard deviation (sigma) parameterization.
#'
#' These functions are optimized for use with \code{\link{loglik}} and include
#' pre-computed derivative rules stored in \code{\link{drule}}.
#'
#' @name normal_distribution
#' @family distribution functions
NULL

#' Normal Probability Density Function
#'
#' Computes the PDF of the Normal distribution.
#'
#' The PDF is given by:
#' \deqn{f(x; \mu, \sigma) = \frac{1}{\sigma\sqrt{2\pi}} \exp\left(-\frac{(x-\mu)^2}{2\sigma^2}\right)}
#'
#' @param x Numeric vector of quantiles.
#' @param mu Numeric. Mean parameter.
#' @param sigma Numeric. Standard deviation parameter (must be positive).
#' @return Numeric vector of density values.
#'
#' @examples
#' normal_pdf(x = -3:3, mu = 0, sigma = 1)
#'
#' @export
normal_pdf <- function(x, mu, sigma) {
  # Normal PDF with mean mu and standard deviation sigma
  (1 / (sigma * sqrt(2 * pi))) * exp(-0.5 * ((x - mu) / sigma)^2)
}

#' Normal Cumulative Distribution Function
#'
#' Computes the CDF of the Normal distribution.
#'
#' The CDF is given by:
#' \deqn{F(x; \mu, \sigma) = \Phi\left(\frac{x - \mu}{\sigma}\right)}
#'
#' where \eqn{\Phi} is the standard normal CDF.
#'
#' @param x Numeric vector of quantiles.
#' @param mu Numeric. Mean parameter.
#' @param sigma Numeric. Standard deviation parameter (must be positive).
#' @return Numeric vector of cumulative probabilities.
#'
#' @examples
#' normal_cdf(x = -3:3, mu = 0, sigma = 1)
#'
#' @export
normal_cdf <- function(x, mu, sigma) {
  # Use built-in pnorm
  pnorm(x, mean = mu, sd = sigma)
}

#' Normal Quantile Function
#'
#' Computes the quantile (inverse CDF) of the Normal distribution.
#'
#' @param p Numeric vector of probabilities.
#' @param mu Numeric. Mean parameter.
#' @param sigma Numeric. Standard deviation parameter (must be positive).
#' @return Numeric vector of quantiles.
#'
#' @examples
#' normal_quantile(p = c(0.025, 0.5, 0.975), mu = 0, sigma = 1)
#'
#' @export
normal_quantile <- function(p, mu, sigma) {
  qnorm(p, mean = mu, sd = sigma)
}

#' Normal Mean
#'
#' Returns the mean parameter of the Normal distribution.
#'
#' @param mu Numeric. Mean parameter.
#' @return Numeric scalar, the mean of the distribution.
#'
#' @examples
#' normal_mean(mu = 5)
#'
#' @export
normal_mean <- function(mu) {
  mu
}

#' Normal Variance
#'
#' Computes the variance of the Normal distribution.
#'
#' The variance is:
#' \deqn{\text{Var}(X) = \sigma^2}
#'
#' @param sigma Numeric. Standard deviation parameter (must be positive).
#' @return Numeric scalar, the variance of the distribution.
#'
#' @examples
#' normal_variance(sigma = 2)
#'
#' @export
normal_variance <- function(sigma) {
  sigma^2
}

#' Normal Median
#'
#' Returns the median of the Normal distribution.
#'
#' For Normal distributions, the median equals the mean.
#'
#' @param mu Numeric. Mean parameter.
#' @return Numeric scalar, the median of the distribution.
#'
#' @examples
#' normal_median(mu = 5)
#'
#' @export
normal_median <- function(mu) {
  mu
}

#' Normal Mode
#'
#' Returns the mode of the Normal distribution.
#'
#' For Normal distributions, the mode equals the mean.
#'
#' @param mu Numeric. Mean parameter.
#' @return Numeric scalar, the mode of the distribution.
#'
#' @examples
#' normal_mode(mu = 5)
#'
#' @export
normal_mode <- function(mu) {
  mu
}

#' Initial Parameter Estimates for Normal Distribution
#'
#' Computes initial parameter estimates for the Normal distribution.
#'
#' Uses simple sample estimates:
#' \deqn{\mu = \bar{x}}
#' \deqn{\sigma = s}
#'
#' where \eqn{\bar{x}} is the sample mean and \eqn{s} is the sample standard deviation.
#'
#' @param x Numeric vector of observations.
#' @param weights Numeric vector of weights/frequencies (optional).
#' @param verbose Logical. Print intermediate calculations (default FALSE).
#' @return Named numeric vector with elements \code{mu} and \code{sigma}.
#'
#' @examples
#' set.seed(123)
#' x <- rnorm(100, mean = 5, sd = 2)
#' normal_init(x)
#'
#' @export
normal_init <- function(x, weights = NULL, verbose = FALSE) {
  # Remove NAs
  x <- x[!is.na(x)]

  # Compute (weighted) mean and SD
  if (is.null(weights)) {
    mu <- mean(x)
    sigma <- sd(x)
  } else {
    w_sum <- sum(weights, na.rm = TRUE)
    mu <- sum(x * weights, na.rm = TRUE) / w_sum
    sigma <- sqrt(sum(weights * (x - mu)^2, na.rm = TRUE) / w_sum)
  }

  if (verbose) {
    message(sprintf("Mean: %.4f, SD: %.4f", mu, sigma))
  }

  c(mu = mu, sigma = sigma)
}

# Compute derivatives for PDF
drule[['normal_pdf']] <- lapply(
  c('x', 'mu', 'sigma'),
  function(var) Deriv::Deriv(body(normal_pdf), var, cache.exp = FALSE)
)
names(drule[['normal_pdf']]) <- c('x', 'mu', 'sigma')

# Compute derivatives for CDF
drule[['normal_cdf']] <- lapply(
  c('x', 'mu', 'sigma'),
  function(var) Deriv::Deriv(body(normal_cdf), var, cache.exp = FALSE)
)
names(drule[['normal_cdf']]) <- c('x', 'mu', 'sigma')
