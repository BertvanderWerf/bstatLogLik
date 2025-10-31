#' Weibull Distribution Functions for Likelihood Estimation
#'
#' Probability density function (PDF), cumulative distribution function (CDF),
#' and utility functions for the three-parameter Weibull distribution with
#' shape (k), scale (lambda), and location (c) parameters.
#'
#' These functions are optimized for use with \code{\link{loglik}} and include
#' pre-computed derivative rules stored in \code{\link{drule}}.
#'
#' @name weibull_distribution
#' @family distribution functions
NULL

#' Weibull Probability Density Function
#'
#' Computes the PDF of the three-parameter Weibull distribution.
#'
#' The PDF is given by:
#' \deqn{f(x; k, \lambda, c) = \frac{k}{\lambda}\left(\frac{x+c}{\lambda}\right)^{k-1} \exp\left(-\left(\frac{x+c}{\lambda}\right)^k\right)}
#'
#' @param x Numeric vector of quantiles.
#' @param k Numeric. Shape parameter (must be positive).
#' @param lambda Numeric. Scale parameter (must be positive).
#' @param c Numeric. Location/offset parameter (default 0).
#' @return Numeric vector of density values.
#'
#' @examples
#' # Standard two-parameter Weibull
#' weibull_pdf(x = 0:10, k = 2, lambda = 5, c = 0)
#'
#' # With location parameter
#' weibull_pdf(x = 0:10, k = 2, lambda = 5, c = 1)
#'
#' @export
weibull_pdf <- function(x, k, lambda, c = 0) {
  # Weibull PDF with shape k, scale lambda, and location c
  (k / lambda) * ((x + c) / lambda)^(k - 1) * exp(-((x + c) / lambda)^k)
}

#' Weibull Cumulative Distribution Function
#'
#' Computes the CDF of the three-parameter Weibull distribution.
#'
#' The CDF is given by:
#' \deqn{F(x; k, \lambda, c) = 1 - \exp\left(-\left(\frac{x+c}{\lambda}\right)^k\right)}
#'
#' @param x Numeric vector of quantiles.
#' @param k Numeric. Shape parameter (must be positive).
#' @param lambda Numeric. Scale parameter (must be positive).
#' @param c Numeric. Location/offset parameter (default 0).
#' @return Numeric vector of cumulative probabilities.
#'
#' @examples
#' weibull_cdf(x = 0:10, k = 2, lambda = 5, c = 0)
#'
#' @export
weibull_cdf <- function(x, k, lambda, c = 0) {
  # Weibull CDF with shape k, scale lambda, and location c
  1 - exp(-((x + c) / lambda)^k)
}

#' Weibull Hazard Function
#'
#' Computes the hazard (failure rate) function of the Weibull distribution.
#'
#' The hazard function is:
#' \deqn{h(x; k, \lambda, c) = \frac{k}{\lambda}\left(\frac{x+c}{\lambda}\right)^{k-1}}
#'
#' @param x Numeric vector of quantiles.
#' @param k Numeric. Shape parameter (must be positive).
#' @param lambda Numeric. Scale parameter (must be positive).
#' @param c Numeric. Location/offset parameter (default 0).
#' @return Numeric vector of hazard rates.
#'
#' @examples
#' weibull_hazard(x = 0:10, k = 2, lambda = 5, c = 0)
#'
#' @export
weibull_hazard <- function(x, k, lambda, c = 0) {
  # Hazard is 0 for x < 0 (below support)
  ifelse(x < 0, 0, (k / lambda) * ((x + c) / lambda)^(k - 1))
}

#' Weibull Mean
#'
#' Computes the expected value (mean) of the Weibull distribution.
#'
#' The mean is:
#' \deqn{E[X] = \lambda \Gamma\left(1 + \frac{1}{k}\right) - c}
#'
#' @param k Numeric. Shape parameter (must be positive).
#' @param lambda Numeric. Scale parameter (must be positive).
#' @param c Numeric. Location/offset parameter (default 0).
#' @return Numeric scalar, the mean of the distribution.
#'
#' @examples
#' weibull_mean(k = 2, lambda = 5, c = 0)
#'
#' @export
weibull_mean <- function(k, lambda, c = 0) {
  lambda * gamma(1 + 1 / k) - c
}

#' Weibull Median
#'
#' Computes the median of the Weibull distribution.
#'
#' The median is:
#' \deqn{\text{Median} = \lambda (\log 2)^{1/k} - c}
#'
#' @param k Numeric. Shape parameter (must be positive).
#' @param lambda Numeric. Scale parameter (must be positive).
#' @param c Numeric. Location/offset parameter (default 0).
#' @return Numeric scalar, the median of the distribution.
#'
#' @examples
#' weibull_median(k = 2, lambda = 5, c = 0)
#'
#' @export
weibull_median <- function(k, lambda, c = 0) {
  lambda * log(2)^(1 / k) - c
}

#' Weibull Mode
#'
#' Computes the mode of the Weibull distribution.
#'
#' The mode is:
#' \deqn{\text{Mode} = \begin{cases}
#' 0 & \text{if } k \leq 1 \\
#' \lambda \left(\frac{k-1}{k}\right)^{1/k} - c & \text{if } k > 1
#' \end{cases}}
#'
#' @param k Numeric. Shape parameter (must be positive).
#' @param lambda Numeric. Scale parameter (must be positive).
#' @param c Numeric. Location/offset parameter (default 0).
#' @return Numeric scalar, the mode of the distribution.
#'
#' @examples
#' weibull_mode(k = 2, lambda = 5, c = 0)
#' weibull_mode(k = 0.5, lambda = 5, c = 0)  # Returns 0 for k <= 1
#'
#' @export
weibull_mode <- function(k, lambda, c = 0) {
  ifelse(k <= 1, 0, lambda * ((k - 1) / k)^(1 / k)) - c
}

#' Initial Parameter Estimates for Weibull Distribution
#'
#' Computes method-of-moments initial parameter estimates for the Weibull distribution,
#' suitable as starting values for \code{\link{loglik}}.
#'
#' Uses empirical approximations based on the mean and standard deviation:
#' \deqn{k \approx (sd/mean)^{-1.086}}
#' \deqn{\lambda \approx mean \cdot k^{2.6674} / (0.184 + 0.816 \cdot k^{2.73855})}
#'
#' Reference: Method-of-moments approximations for Weibull parameters.
#'
#' @param x Numeric vector of observations (or difftime object).
#' @param weights Numeric vector of weights/frequencies (optional). If provided, weighted mean and variance are used.
#' @param verbose Logical. Print intermediate calculations (default FALSE).
#' @return Named numeric vector with elements \code{k}, \code{lambda}, and \code{c}.
#'
#' @examples
#' set.seed(123)
#' x <- rweibull(100, shape = 2, scale = 5)
#' weibull_init(x)
#'
#' # With weights
#' weibull_init(x, weights = rep(1:2, 50))
#'
#' @export
weibull_init <- function(x, weights = NULL, verbose = FALSE) {
  # Convert difftime objects to numeric
  if (inherits(x, "difftime")) {
    x <- as.numeric(x)
  }

  # Compute (weighted) mean and standard deviation
  if (is.null(weights)) {
    m <- mean(x, na.rm = TRUE)
    v <- sd(x, na.rm = TRUE)
  } else {
    # Weighted mean and SD
    w_sum <- sum(weights, na.rm = TRUE)
    m <- sum(x * weights, na.rm = TRUE) / w_sum
    v <- sqrt(sum(weights * x^2, na.rm = TRUE) / w_sum - m^2)
  }

  if (verbose) {
    message(sprintf("Mean: %.4f, SD: %.4f", m, v))
  }

  # Method-of-moments approximations
  cv <- v / m  # Coefficient of variation
  k <- cv^(-1.086)
  lambda <- m * k^2.6674 / (0.184 + 0.816 * k^2.73855)

  c(k = k, lambda = lambda, c = 0)
}

# Compute derivatives for PDF
drule[['weibull_pdf']] <- lapply(
  c('x', 'k', 'lambda', 'c'),
  function(var) Deriv::Deriv(body(weibull_pdf), var, cache.exp = FALSE)
)
names(drule[['weibull_pdf']]) <- c('x', 'k', 'lambda', 'c')

# Compute derivatives for CDF
drule[['weibull_cdf']] <- lapply(
  c('x', 'k', 'lambda', 'c'),
  function(var) Deriv::Deriv(body(weibull_cdf), var, cache.exp = FALSE)
)
names(drule[['weibull_cdf']]) <- c('x', 'k', 'lambda', 'c')
