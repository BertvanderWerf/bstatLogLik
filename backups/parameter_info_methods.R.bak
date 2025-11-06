
#' Title
#'
#' @param parameter_info
#' @param data
#' @param ...
#'
#' @returns
#' @export
#'
#' @examples
as.data.frame.parameter_info <- function(parameter_info, data=NULL, ...) {

  # --- Input validation ---
  # * parameter_info *
  if (!inherits(parameter_info, 'parameter_info')) {
    stop("Argument 'parameter_info' must be an object of class 'parameter_info'", call. = FALSE)
  }

  # * data *
  if (is.null(data)) {
    if (is.null(parameter_info$model_data)) {
      stop('data must be set or model_data in parameter_info must be set', call. = FALSE)
    }
  } else {
    bstatErr::check_data_frame(data)
    parameter_info <- create_parameter_info(parameter_info, model_data=data)
  }

  parnames <- names(parameter_info$parameter_values)
  npars <- length(parnames)
  res <- NULL
  for (i in seq_len(npars)) {
    res <- cbind(res, parameter_info$design_matrices[[i]] %*%  parameter_info$parameter_values[[i]])
  }
  colnames(res) <- parnames
  as.data.frame(res)
}

#' Title
#'
#' @param parameter_info
#' @param n
#' @param vcov
#'
#' @returns
#' @export
#'
#' @examples
sample_parameter_info <- function(parameter_info, n, vcov) {
  # --- Input validation ---
  # * parameter_info *
  if (!inherits(parameter_info, 'parameter_info')) {
    stop("Argument 'parameter_info' must be an object of class 'parameter_info'", call. = FALSE)
  }

  parnames <- names(parameter_info$parameter_values)
  npars <- length(parnames)
  parnames_full <- paste(rep(names(parameter_info$parameter_values), sapply(parameter_info$parameter_values, length)),
                         unlist(sapply(parameter_info$parameter_values, names)),
                         sep=parameter_info$sep)
  npars_full <- length(parnames_full)
  cond_selected <- as.logical(unlist(parameter_info$estimate_flag))
  parnames_selected <- parnames_full[cond_selected]
  npars_selected <- sum(cond_selected)

  vcov <- vcov(fit)
  # * data *
  if (is.null(vcov) || !is.matrix(vcov) || ncol(vcov) != nrow(vcov)) { # || !(all(vcov==t(vcov)))) {
      stop('vcov must be a symmetric variance-covariance matrix ', call. = FALSE)
  }
  if (ncol(vcov)< npars_full) {
    v <- matrix(0, nrow=npars_full, ncol=npars_full,
                dimnames=list(parnames_full, parnames_full))
    c <- match(row.names(v), row.names(vcov), nomatch = 0)
    v[c, c] <- vcov
    vcov <- v
  }

  parvalues <- unlist(parameter_info$parameter_values)
  parvalues <- rbind(parvalues,MASS::mvrnorm(n, parvalues, vcov))
  row.names(parvalues) <- c('fitted', paste('simulation', 1:n, sep="_"))

  subnames <- unlist(sapply(parameter_info$parameter_values, names))
  groups <- factor(rep(parnames, sapply(parameter_info$parameter_values, length)), parnames)
  indices <- setNames(1:npars_full, unlist(sapply(parameter_info$parameter_values, names)))
  indices <- split(indices, groups)

  L <- list()
  for (i in 1:(n+1)) {
    L_i <- indices
    for (j in 1:length(indices)) {
      L_i[[j]][] <- parvalues[i,indices[[j]]]
    }
    L[[i]] <- L_i
  }
  names(L) <- row.names(parvalues)
  L
}
