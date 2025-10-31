#' Control Structure for Maximum Likelihood Optimization
#'
#' S3 class defining control parameters for likelihood-based optimization routines.
#'
#' Objects of this class contain:
#'   \itemize{
#'     \item{\code{method}: Optimization algorithm - "NewtonRaphson", "Levenberg", or "Marquardt".}
#'     \item{\code{criterium}: Character vector of convergence criteria to evaluate: "delta" (parameter change), "jacobian" (gradient norm), "lnlik" (log-likelihood change).}
#'     \item{\code{hessian}: Method for computing Hessian matrix - "calculate" (analytic second derivatives) or "expected" (outer product of scores).}
#'     \item{\code{eps}: Numeric tolerance threshold for convergence.}
#'     \item{\code{nsteps}: Maximum number of optimization iterations.}
#'     \item{\code{lambda}: Initial damping parameter for Levenberg-Marquardt methods.}
#'     \item{\code{lambda_multiplier}: Factor for adjusting lambda during optimization.}
#'     \item{\code{reduce}: Logical flag to reduce dataset to only rows needed for parameter estimation.}
#'     \item{\code{replace_nan}: Logical flag to replace NaN values with zero in derivative calculations.}
#'     \item{\code{use_ginv}: Logical flag to use generalized inverse (MASS::ginv) instead of solve() for singular matrices.}
#'     \item{\code{replace_inf}: Logical flag to replace infinite values with machine precision limits.}
#'   }
#'
#' @name loglik_control_class
#' @aliases loglik_control-class
#' @docType class
#'
#' @seealso \code{\link{loglik_control}}
#' @exportClass loglik_control

#' Create Control Parameters for Maximum Likelihood Optimization
#'
#' Constructs a \code{loglik_control} S3 object specifying algorithmic choices,
#' convergence criteria, and numerical stability settings for \code{\link{loglik}}.
#'
#' The Newton-Raphson method uses the standard update rule. Levenberg adds a
#' diagonal damping term \code{lambda * I}, while Marquardt uses \code{lambda * diag(H)}.
#' Multiple convergence criteria can be specified simultaneously; optimization stops
#' when any criterion falls below \code{eps}.
#'
#' @param method Character. Optimization algorithm: "NewtonRaphson" (default), "Levenberg", or "Marquardt".
#' @param criterium Character vector. Convergence criteria: "delta", "jacobian", "lnlik" (default all three).
#' @param eps Numeric. Convergence tolerance (default 1e-7).
#' @param nsteps Integer. Maximum iterations (default 200).
#' @param hessian Character. Hessian computation: "calculate" (analytic, default) or "expected" (outer product).
#' @param lambda Numeric. Initial damping parameter for Levenberg/Marquardt (default 1).
#' @param lambda_multiplier Numeric. Factor for lambda adjustment during optimization (default 10).
#' @param reduce Logical. Reduce data to relevant rows only (default TRUE).
#' @param replace_nan Logical. Replace NaN with 0 in derivatives (default FALSE).
#' @param use_ginv Logical. Use generalized inverse for singular matrices (default FALSE).
#' @param replace_inf Logical. Replace Inf with machine limits (default FALSE).
#'
#' @return An S3 object of class \code{loglik_control}, see \code{\link{loglik_control_class}} for structure.
#'
#' @examples
#' # Standard Newton-Raphson with default settings
#' ctrl1 <- loglik_control()
#'
#' # Levenberg-Marquardt with stricter convergence
#' ctrl2 <- loglik_control(method = "Marquardt", eps = 1e-10)
#'
#' # Robust settings for difficult problems
#' ctrl3 <- loglik_control(
#'   use_ginv = TRUE,
#'   replace_nan = TRUE,
#'   replace_inf = TRUE,
#'   nsteps = 500
#' )
#'
#' @export
loglik_control <- function(
    method = c("NewtonRaphson", "Levenberg", "Marquardt"),
    criterium = c("delta", "jacobian", "lnlik"),
    eps = 1e-7,
    nsteps = 200,
    hessian = c("calculate", "expected"),
    lambda = 1,
    lambda_multiplier = 10,
    reduce = TRUE,
    replace_nan = FALSE,
    use_ginv = FALSE,
    replace_inf = FALSE
) {
  result <- list(
    method = match.arg(method),
    criterium = criterium,
    hessian = match.arg(hessian),
    eps = eps,
    nsteps = nsteps,
    lambda = lambda,
    lambda_multiplier = lambda_multiplier,
    reduce = reduce,
    replace_nan = replace_nan,
    use_ginv = use_ginv,
    replace_inf = replace_inf
  )

  class(result) <- c("loglik_control", class(result))
  result
}
