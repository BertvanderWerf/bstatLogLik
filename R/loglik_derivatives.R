#' Compute Jacobian (Gradient) Matrix for Log-Likelihood
#'
#' Assembles the gradient matrix from automatic differentiation output,
#' accounting for parameter design matrices and conditional estimation.
#'
#' @param deriv_output List. Output from Deriv::Deriv containing derivative components.
#' @param parameter_info Object of class \code{parameter_info}.
#' @param cond_vec Logical vector. Which parameters are being estimated.
#' @param npars Integer. Total number of parameter groups.
#' @param parnames_selected Character vector. Names of parameters being estimated.
#' @param control List from loglik_control().
#' @return Numeric matrix of gradients (n_obs x n_params).
#' @export
compute_gradient <- function(deriv_output, parameter_info, cond_vec, npars, parnames_selected, control) {
  d1 <- deriv_output$`1`

  # Ensure d1 is a matrix
  if (is.null(dim(d1))) {
    dim(d1) <- c(length(d1), 1)
  }

  k <- 0
  jacobian <- NULL

  # Loop over parameter groups
  for (j in seq_len(npars)) {
    if (cond_vec[j]) {
      k <- k + 1
      mat_j <- parameter_info$design_matrices[[j]]

      if (!is.null(mat_j)) {
        # Multiply design matrix columns by derivative
        jacobian <- cbind(jacobian, mat_j[, parameter_info$estimate_flag[[j]] == TRUE, drop = FALSE] * d1[, k])
      } else {
        jacobian <- cbind(jacobian, d1[, k, drop = FALSE])
      }
    }
  }

  colnames(jacobian) <- parnames_selected
  adjust_nan_inf(jacobian, control)
}



#' Compute Hessian (Second Derivative) Matrix for Log-Likelihood
#'
#' Assembles the Hessian from automatic differentiation or outer product of scores,
#' with optional Levenberg-Marquardt damping.
#'
#' @param deriv_output List. Output from Deriv::Deriv.
#' @param gradient Matrix. Computed gradient/Jacobian.
#' @param parameter_info Object of class \code{parameter_info}.
#' @param cond_vec Logical vector. Which parameters are being estimated.
#' @param npars Integer. Total number of parameter groups.
#' @param parnames_selected Character vector. Names of parameters being estimated.
#' @param control List from loglik_control().
#' @return Symmetric Hessian matrix (n_params x n_params).
#' @export
compute_hessian <- function(deriv_output, gradient, parameter_info, cond_vec, npars, parnames_selected, control) {

  if (control$hessian == "calculate") {
    # Analytic Hessian from second derivatives
    d2 <- deriv_output$`2`

    if (is.null(dim(d2))) {
      dim(d2) <- c(length(d2), 1)
    }

    # Build nested list structure for Hessian blocks
    hess_list <- matrix(list(), nrow = sum(cond_vec), ncol = sum(cond_vec))

    k <- 0
    for (j in seq_len(npars)) {
      for (i in seq_len(npars)) {
        if (cond_vec[i] && cond_vec[j]) {
          k <- k + 1
          mat_i <- parameter_info$design_matrices[[i]]
          mat_j <- parameter_info$design_matrices[[j]]
          cond_i <- parameter_info$estimate_flag[[i]]
          cond_j <- parameter_info$estimate_flag[[j]]

          # if (is.null(mat_i) && is.null(mat_j)) {
          #   # Both parameters are scalars
          #   hess_list[[k]] <- sum(d2[, k, drop = FALSE])
          #   dim(hess_list[[k]]) <- c(1, 1)
          #
          # } else if (is.null(mat_i)) {
          #   # i is scalar, j has design matrix
          #   hess_list[[k]] <- t(colSums(d2[, k] * mat_j[, cond_j == TRUE, drop = FALSE]))
          #
          # } else if (is.null(mat_j)) {
          #   # i has design matrix, j is scalar
          #   hess_list[[k]] <- t(t(colSums(d2[, k] * mat_i[, cond_i == TRUE, drop = FALSE])))
          #
          # } else {
            # Both have design matrices
            n_i <- sum(cond_i)
            n_j <- sum(cond_j)
            h_block <- matrix(NA, nrow = n_i, ncol = n_j)

            kk <- 0
            for (jj in seq_len(ncol(mat_j))) {
              if (cond_j[jj]) {
                for (ii in seq_len(ncol(mat_i))) {
                  if (cond_i[ii]) {
                    kk <- kk + 1
                    h_block[[kk]] <- sum(d2[, k] * mat_i[, ii] * mat_j[, jj])
                  }
                }
              }
            }
            hess_list[[k]] <- h_block
          # }
        }
      }
    }

    # Assemble full Hessian from blocks
    hessian <- NULL
    k <- 0
    for (i in seq_len(nrow(hess_list))) {
      row_blocks <- NULL
      for (j in seq_len(ncol(hess_list))) {
        k <- k + 1
        row_blocks <- rbind(row_blocks, hess_list[[k]])
      }
      hessian <- cbind(hessian, row_blocks)
    }

    rownames(hessian) <- colnames(hessian) <- parnames_selected

    # Symmetrize
    hessian <- (hessian + t(hessian)) / 2

  } else {
    # Expected information: outer product of gradients
    hessian <- crossprod(gradient)
  }

  # Add damping for Levenberg-Marquardt methods
  n_param <- length(parnames_selected)
  if (control$method == "Levenberg") {
    hessian <- hessian + control$lambda * diag(nrow = n_param)
  } else if (control$method == "Marquardt") {
    hessian <- hessian + control$lambda * diag(diag(hessian))
  }

  adjust_nan_inf(hessian, control)
}
