#' Extract Score Contributions (estfun) for Sandwich Estimator
#'
#' Returns the matrix of score contributions (one row per observation).
#'
#' @param object Object of class "loglik".
#' @return Matrix of score contributions.
#' @export
estfun <- function(object) {
  UseMethod("estfun")
}

#' @export
estfun.loglik <- function(object) {
  # Return the stored gradient matrix from the fitted model
  object$gradient_matrix
}


#' Compute Meat Matrix for Sandwich Estimator
#'
#' Calculates the "meat" component: sum of outer products of scores.
#'
#' @param object Object of class "loglik".
#' @return Meat matrix.
#' @export
meat.loglik <- function(object) {
  scores <- estfun.loglik(object)
  crossprod(scores)
}


#' Compute Bread Matrix for Sandwich Estimator
#'
#' Returns the "bread" component: negative of the Hessian.
#'
#' @param object Object of class "loglik".
#' @return Bread matrix.
#' @export
bread.loglik <- function(object) {
  -object$hessian
}


#' Robust Sandwich Variance-Covariance Matrix
#'
#' Computes the sandwich (robust) covariance estimator.
#'
#' @param object Object of class "loglik".
#' @return Sandwich variance-covariance matrix.
#' @export
vcov_sandwich.loglik <- function(object) {
  bread <- bread.loglik(object)
  meat <- meat.loglik(object)

  bread_inv <- if (object$control$use_ginv) {
    solve_general(bread)
  } else {
    solve(bread)
  }

  bread_inv %*% meat %*% t(bread_inv)
}
