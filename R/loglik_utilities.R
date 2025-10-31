#' Adjust NaN and Inf Values in Derivative Objects
#'
#' Replaces NaN and infinite values according to control settings.
#'
#' @param x Numeric vector, matrix, or list.
#' @param control List from loglik_control().
#' @return Cleaned object.
#' @keywords internal
adjust_nan_inf <- function(x, control) {
  if (is.list(x)) {
    return(lapply(x, adjust_nan_inf, control = control))
  }

  if (isTRUE(control$replace_nan)) {
    x[is.na(x)] <- 0
  }

  if (isTRUE(control$replace_inf)) {
    x[is.infinite(x) & sign(x) < 0] <- .Machine$double.xmin
    x[is.infinite(x) & sign(x) > 0] <- .Machine$double.xmax
  }

  x
}

#' Generalized Matrix Solver with Pseudo-Inverse Fallback
#'
#' Solves Ax = b using standard solver or generalized inverse if singular.
#'
#' @param a Numeric matrix.
#' @param b Numeric vector or NULL.
#' @param tol Numeric. Tolerance for pseudo-inverse computation.
#' @return Inverse matrix (if b is NULL) or solution vector x.
#' @keywords internal
solve_general <- function(a, b = NULL, tol = sqrt(.Machine$double.eps)) {
  ginv <- MASS::ginv(a, tol)
  if (is.null(b)) return(ginv)
  res <- c(ginv %*% b)
  names(res) <- colnames(a)
  res
}
